#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(parallel)
library(patchwork)

#parse arguments 
parser <- ArgumentParser(description= "sample 100000 random snps and get the odds ratio of each genomic region for a critinal qvalue 10000 times")
#inputs
parser$add_argument('--incr', '-incr', help = "table with clinal SNPs which slopes are increasing")
parser$add_argument('--q97', '-q97', help= 'q-values_97 file')
parser$add_argument('--q0910', '-q0910', help= 'q-values_0910 file')
#outputs
parser$add_argument('--output', '-out', help= 'table with odds ratio for all 10000 trials')
parser$add_argument('--plot', '-plot', help= 'plot with odds ratio')
xargs<- parser$parse_args()

#read data --------------------------------------------------------------------

# q97 <- fread("/dados/time_clines/analysis/time_GLM_lat/effects_with_shortI_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv")
# q0910 <- fread("/dados/time_clines/analysis/time_GLM_lat/effects_with_shortI_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910all.tsv")
# 
# incr <- fread("/dados/time_clines/analysis/time_GLM_lat/incerasing_slope_snps_all.tsv",
#               select = c(1,2,12))

q97 <- fread(xargs$q97)
q0910 <- fread(xargs$q0910)
incr <- fread(xargs$incr,
              select = c(1,2,12))

#get non concordant SNPs ----------------------------------------------------

setkey(q97, chrom, POS)
setkey(q0910, chrom, POS)

setcolorder(q97, c("chrom", "POS"))
setcolorder(q0910, c("chrom", "POS"))

newnames97 <- paste0(colnames(q97)[-1:-2], "_97")
newnames0910 <- paste0(colnames(q0910)[-1:-2], "_0910")

setnames(q97, colnames(q97)[-1:-2], newnames97)
setnames(q0910, colnames(q0910)[-1:-2], newnames0910)

all_q <- 
  merge.data.table(q97, q0910)

only_clinal_0910 <- 
  all_q[qvalue_0910 <= 0.1 & qvalue_97 > 0.1]

#get test tables -----------------------------------------------------------

effect_clinal0910 <- 
  merge.data.table(q97[, .(chrom, POS, effect_97)], 
                   only_clinal_0910[, .(chrom, POS, qvalue_0910)], 
                   by.x = c("chrom", "POS"), 
                   by.y = c("chrom", "POS"), all = TRUE)

effect_incr <- 
merge.data.table(q97[, .(chrom, POS, effect_97)], incr, by.x = c("chrom", "POS"), 
                 by.y = c("CHROM", "POS"), all = TRUE)



#calcula a chance de um snp clinal ser daquela determindade classe genômica
chance_target <- function(data, column, effect_column) {
  snps <- data[!is.na(get(column)), ]  
  n_snps <- snps[, .(n_r = .N), by = effect_column]
  n_snps[, total_clinal := sum(n_r) - n_r]
  n_snps[, chance_r := n_r / total_clinal]
  
  return(n_snps)
}


n_snps_clinal0910 <- 
chance_target(effect_clinal0910, "qvalue_0910", "effect_97")

n_snps_incr <- 
  chance_target(effect_incr, "diff_slope", "effect_97")


#pega 100000 snps aleatórios e calcula a chance dele daquela determindade classe genômica
get_controls <- 
  function(dt, column, effect_column){
    control_snps <- dt[is.na(get(column))][sample(.N, 100000)]
    
    n_control <- control_snps[,.(n = nrow(.SD)), by = effect_column]
    n_control[, total_n_class := sum(n_control$n)-n]
    n_control[, chance_control := n/total_n_class]
    return(n_control)
  }


num_cores <- 30


#repete isso 10000 vezes
all_controls_incr <- 
  mclapply(1:10000, function(x) get_controls(effect_incr, "diff_slope", "effect_97"),
           mc.cores = num_cores)

all_controls_clinal0910 <- 
  mclapply(1:10000, function(x) get_controls(effect_clinal0910, "qvalue_0910", "effect_97"),
           mc.cores = num_cores)


all_controls_incr <- 
  mclapply(all_controls_incr, function(x) merge.data.table(n_snps_incr, x, 
                                                      by = "effect_97"),
           mc.cores = num_cores) |> rbindlist()

all_controls_clinal0910 <- 
  mclapply(all_controls_clinal0910, function(x) merge.data.table(n_snps_clinal0910, x, 
                                                           by = "effect_97"),
           mc.cores = num_cores) |> rbindlist()

all_controls_incr[, odds_ratio := chance_r/chance_control]
all_controls_clinal0910[, odds_ratio := chance_r/chance_control]


#saves odd ratio tables ----------------------------------------------------------
all_controls_incr[, data_from := "increasing_slopes"]
all_controls_clinal0910[, data_from := "only_clinal0910"]

out <- 
rbindlist(list(all_controls_incr,all_controls_clinal0910))

fwrite(out, file = xargs$output, sep = "\t")

out[, effect_plot := str_replace_all(effect_97, "_", " ")]
out[, effect_plot := fcase(effect_97 == "NON_SYNONYMOUS_CODING", "NSYN",
                                         effect_97 == "SYNONYMOUS_CODING", "SYN",
                                         effect_97 == "UTR_3", "3'UTR",
                                         effect_97 == "UTR_5", "5'UTR",
                                         effect_97 == "SHORT_INTRON", "SHORT INTRON",
                                         rep(TRUE, .N), effect_97)]

PLOT_incr <- 
  out[effect_97 != "-" & data_from == "increasing_slopes"] |>
  ggplot(aes(x= odds_ratio)) +
  geom_histogram(aes(),  binwidth = 0.005) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  facet_grid(cols = vars(effect_plot), #rows = vars(clinal_year),
             scales = "free_x") +
  theme_bw() +
  labs(x = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))


PLOT_clinal0910 <- 
  out[effect_97 != "-" & data_from == "only_clinal0910"] |>
  ggplot(aes(x= odds_ratio)) +
  geom_histogram(aes(),  binwidth = 0.005) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  facet_grid(cols = vars(effect_plot), #rows = vars(clinal_year),
             scales = "free_x") +
  theme_bw() +
  labs(x = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))



pdf(file = xargs$plot,
    width = 11.5,
    height = 5.5)

(PLOT_incr + labs(title = "Increased Slopes")) /
  (PLOT_clinal0910 + labs(title = "SNPs only clinal in 2009/20010"))

dev.off()







