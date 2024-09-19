#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)

#parse arguments 
parser <- ArgumentParser(description= "sample 100000 random snps and get the odds ratio of each genomic region for a critinal qvalue 10000 times")
parser$add_argument('--odr97', '-odr97', help= 'table with odds ratio for all 10000 trials in 1997')
parser$add_argument('--odr0910', '-odr0910', help= 'table with odds ratio for all 10000 trials in 2009/2010')
parser$add_argument('--out', '-out', help= 'plot - odds ratio distribution for each genomic region')
xargs<- parser$parse_args()

# odds_ratio97 <- fread("/dados/time_clines/analysis/time_GLM_lat/genomic_region_odr/odds_ratio_97_clinal_at_0.1.tsv")
# odds_ratio0910 <- fread("/dados/time_clines/analysis/time_GLM_lat/genomic_region_odr/odds_ratio_0910all_clinal_at_0.1.tsv")
odds_ratio97 <- fread(xargs$odr97)
odds_ratio0910 <- fread(xargs$odr0910)
odds_ratio97[, clinal_year := "1997"][, i := .I]
odds_ratio0910[, clinal_year := "2009/2010"][, i := .I]

all_or <- rbind(odds_ratio97, odds_ratio0910)
all_or[, effect_plot := str_replace_all(effect, "_", " ")]
all_or[, effect_plot := fcase(effect == "NON_SYNONYMOUS_CODING", "NSYN",
                              effect == "SYNONYMOUS_CODING", "SYN",
                              effect == "UTR_3", "3'UTR",
                              effect == "UTR_5", "5'UTR",
                              rep(TRUE, .N), effect)]

PLOT <- 
all_or[effect != "-"] |>
  ggplot(aes(x= odds_ratio)) +
  geom_histogram(aes(),  binwidth = 0.005) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  facet_grid(cols = vars(effect_plot), rows = vars(clinal_year),
             scales = "free_x") +
  theme_bw() +
  labs(x = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))

jpeg(filename = xargs$out,
     width = 25,
     height = 15,
     units = "cm",
     res = 500)

PLOT

dev.off()
