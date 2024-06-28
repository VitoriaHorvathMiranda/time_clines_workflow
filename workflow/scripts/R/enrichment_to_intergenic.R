#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(ggeffects)
library(MASS)

#different from enrichment_effects script, this script compares each genic class ("effect")
# to the intergenic class

#parse arguments 
parser <- ArgumentParser(
  description="makes enrichiment plots and enchriment table, based on intergenic regions")
parser$add_argument('--eff', '-eff', help = "table with qvalues and effects")
parser$add_argument('--enrPlot', '-pen', help= 'enrichment plot (jpeg)')
parser$add_argument('--enOut', '-out', 
                    help= 'table (tsv) with enrichment coeficient and confidance intervals, all these values ARE NOT exponentiated')
xargs<- parser$parse_args()

# get clinal year ---------------------------------------------------------

file_name <- xargs$eff
year <- 
  file_name |> str_split_i(pattern = "_", i = 14) |> str_remove(".tsv")

year <- ifelse(year == "97", "1997", "2009/2010")

# wrangle q_eff file ------------------------------------------------------

q_eff <- fread(xargs$eff,
               select = c("chrom", "POS", "qvalue", "effect"))

q_eff <- q_eff[effect != "-"]

effects <- unique(q_eff$effect)

q_eff[, effect := factor(effect,
                         levels = c("INTERGENIC","INTRON", "UPSTREAM","-",                  
                                    "SPLICE","UTR_3","NON_SYNONYMOUS_CODING",
                                    "SYNONYMOUS_CODING","UTR_5","DOWNSTREAM"))]


#cria uma coluna para cada limite de qvalue 
FDR <- c(0.1, 0.075, 0.05)
for (i in seq_along(FDR)) {
  
  q_eff[, paste0("clinal_at_", FDR[[i]]) := 
          ifelse(qvalue <= FDR[[i]], 1, 0)]
  
}

clinal_status <- colnames(q_eff)[(length(q_eff)-2):length(q_eff)]
models <- vector("list", length = length(FDR))
conf_int <- vector("list", length(FDR))
for (i in seq_along(clinal_status)) {
  
  test <- paste0(clinal_status[[i]], " ~ effect")
  print(test)
  models[[i]] <- glm(data = q_eff,
               formula = test,
               family = binomial)
  conf_int[[i]] <- confint(models[[i]])
  conf_int[[i]] <- as.data.table(conf_int[[i]])
  models[[i]] <- summary(models[[i]])
  models[[i]] <- coef(models[[i]])
  models[[i]] <- as.data.table(models[[i]])
  models[[i]][, c("lower_ci", "upper_ci") := conf_int[[i]]]
  models[[i]][, clinal_at := FDR[[i]]]
}

odds_ratio <- rbindlist(models)


# plot enrichment ---------------------------------------------------------

odds_ratio[, effect2 := rep(str_replace_all(effects, pattern = "_", " "),
                             times = 3)] 

fwrite(odds_ratio, file = xargs$enOut, sep = "\t")

# exponencia os valores para tudo ficar em razão de chance
exp_values <- c("exp_estimate", "exp_lower_ci", "exp_upper_ci")
cols <- c("Estimate", "lower_ci", "upper_ci")
odds_ratio[,  (exp_values) := lapply(.SD, exp),
           .SDcols = cols]

q_values_names <- c("0.05" = "5%",
                    "0.075" = "7.5%",
                    "0.1" = "10%")

enrichment_plot <- 
  odds_ratio[effect2 != "INTERGENIC"] |>
  ggplot(aes(x = exp_estimate, y = effect2)) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  geom_errorbar(aes(xmax = exp_upper_ci, xmin = exp_lower_ci), 
                linewidth = 0.5,
                width = 0.2,
                color = "grey50") +
  geom_point(color = "firebrick") +
  facet_wrap(vars(as.character(clinal_at)), labeller = as_labeller(q_values_names)) +
  theme_minimal() +
  labs(y = "", x = "Odds Ratio", 
       title = paste("Enrichment of clinal SNPs in", year)) 


jpeg(filename = xargs$enrPlot,
     width = 30,
     height = 15,
     units = "cm",
     res = 1200)

enrichment_plot

dev.off()


