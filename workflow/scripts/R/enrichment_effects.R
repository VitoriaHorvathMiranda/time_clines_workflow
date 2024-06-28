#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(ggeffects)
library(MASS)
library(RColorBrewer)
library(patchwork)

#parse arguments 
parser <- ArgumentParser(description= "makes enrichiment plots and enchriment table")
parser$add_argument('--eff', '-eff', help = "table with qvalues and effects")
parser$add_argument('--enrPlot', '-pen', help= 'enrichment plot (jpeg)')
parser$add_argument('--enOut', '-out', 
                    help= 'table (tsv) with enrichment coeficient and confidance intervals, all these values ARE NOT exponentiated')
parser$add_argument('--CountPlot', '-cp', help= 'bar and proportion plot (jpeg)')
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

dcasted_eff <- #cada efeito vira uma coluna
  dcast.data.table(q_eff, chrom + POS + qvalue ~ effect,
                   fill = 0,
                   value.var = "effect")

#se o snp for aquele efeito 1 se não for 0
dcasted_eff[, (effects) := lapply(.SD, function(x) ifelse(x != 0, 1, x)),
            .SDcols = effects]

#transforma os efeitos em numéricos
dcasted_eff[, (effects) := lapply(.SD, as.numeric),
            .SDcols = effects]

#cria uma coluna para cada limite de qvalue 
FDR <- c(0.1, 0.075, 0.05)
for (i in seq_along(FDR)) {
  
  dcasted_eff[, paste0("clinal_at_", FDR[[i]]) := 
                ifelse(qvalue <= FDR[[i]], 1, 0)]
    
}

clinal_status <- colnames(dcasted_eff)[(length(dcasted_eff)-2):length(dcasted_eff)]


models <- vector("list", length = length(effects))
conf_int <- vector("list", length = length(effects))
output <- vector("list", length = length(effects))
for (i in seq_along(effects)) { #para cada um dos efeitos
  for (j in seq_along(clinal_status)) { #e cada clinal_status
    
    test <- paste0(clinal_status[[j]], " ~ ", effects[[i]]) #roda um modelo
    models[[i]][[j]] <- glm(data = dcasted_eff, test, family = binomial)
    
    conf_int[[i]][[j]] <- confint(models[[i]][[j]], level = 0.95) #calcula os intervalos de confiança
    conf_int[[i]][[j]] <- as.data.table(conf_int[[i]][[j]])
    conf_int[[i]][[j]] <- conf_int[[i]][[j]][2] # drop os intervalos do intercepto
    
    output[[i]][[j]] <- summary(models[[i]][[j]]) #pega os coeficientes de cada modelo
    output[[i]][[j]] <- coef(output[[i]][[j]])
    output[[i]][[j]] <- as.data.table(output[[i]][[j]])
    output[[i]][[j]] <- output[[i]][[j]][2] #drop os coeficinebtes do intercepto
    
  }
  
  output[[i]] <- rbindlist(output[[i]])
  conf_int[[i]] <- rbindlist(conf_int[[i]])
}

#cria as colunas extras para a tabela
test_effect <- rep(effects, each = 3, length.out = (9*3))
test_clinal <- rep(FDR, times = 9)

odds_ratio <- rbindlist(output)
all_int <- rbindlist(conf_int)


odds_ratio[, effect := test_effect]
odds_ratio[, test_clinal := test_clinal]
odds_ratio[, c("lower_ci", "upper_ci") := all_int]

fwrite(odds_ratio, file = xargs$enOut, sep = "\t")


# plot enrichment ---------------------------------------------------------

odds_ratio[, effect2 := str_replace_all(effect, pattern = "_", " ")] 

# exponencia os valores para tudo ficar em razão de chance
exp_values <- c("exp_estimate", "exp_lower_ci", "exp_upper_ci")
cols <- c("Estimate", "lower_ci", "upper_ci")
odds_ratio[,  (exp_values) := lapply(.SD, exp),
           .SDcols = cols]

q_values_names <- c("0.05" = "5%",
                    "0.075" = "7.5%",
                    "0.1" = "10%")

enrichment_plot <- 
odds_ratio |>
  ggplot(aes(x = exp_estimate, y = effect2)) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  geom_errorbar(aes(xmax = exp_upper_ci, xmin = exp_lower_ci), 
                linewidth = 0.5,
                width = 0.2,
                color = "grey50") +
  geom_point(color = "firebrick") +
  facet_wrap(vars(test_clinal), labeller = as_labeller(q_values_names)) +
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


# summary plots -----------------------------------------------------------

q_eff[, clinal_status := fcase(qvalue <= 0.05, 0.05,
                               qvalue <= 0.1, 0.1,
                               qvalue <= 0.2, 0.2,
                               qvalue <= 0.3, 0.3,
                               qvalue <= 0.4, 0.4,
                               qvalue <= 0.5, 0.5,
                               default = 1)]

q_eff[, clinal_status := as.factor(clinal_status)]

q_eff[, effect2 := str_replace_all(effect, pattern = "_", " ")] 

COUNT <- 
q_eff |>
  ggplot() +
  geom_bar(aes(y = effect2,  fill = clinal_status), 
           position = position_stack(reverse = TRUE)) +
  # scale_fill_manual(labels = c("5%", "10%", "none"), 
  #                   values = c("firebrick", "darkolivegreen", "darkslategray")) +
  scale_fill_brewer(palette = "Reds") +
  theme_minimal() +
  labs(y = "") +
  #theme(text = element_text(family = "serif")) +
  guides(fill=guide_legend(title="q-value <=")) +
  theme(axis.text.y = element_blank())


PROP <- 
q_eff |>
  ggplot() +
  geom_bar(aes(y = effect2,  fill = clinal_status), 
           #position = position_stack(reverse = TRUE),
           position = "fill") +
  # scale_fill_manual(labels = c("5%", "10%", "none"), 
  #                   values = c("firebrick", "darkolivegreen", "darkslategray")) +
  scale_fill_brewer(palette = "Reds") +
  theme_minimal() +
  labs(y = "", title = year, 
       x = "proportion") +
  theme(legend.position = "none") #+
  #guides(fill=guide_legend(title="q-value <="))

jpeg(xargs$CountPlot,
     width = 27,
     height = 15,
     units = "cm",
     res = 1200)

PROP + COUNT

dev.off()
