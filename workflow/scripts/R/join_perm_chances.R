#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

#parse arguments 
parser <- ArgumentParser(description= "joins all permuted effect chances")
parser$add_argument('--cyear', '-cyear', help = "clinal year permuted")
parser$add_argument('--chances', '-pc', help='path to chances files')
parser$add_argument('--out', '-out', help= 'table with all chances values')
xargs<- parser$parse_args()

year <- xargs$cyear

perm_chances_file_names <- list.files(xargs$chances,
                                      full.names = TRUE,
                                      pattern = paste0("chance_perm_n_[0-9]{1,4}_",year,".tsv"))


# perm_chances_file_names <- list.files("/dados/time_clines/analysis/time_GLM_lat/Perm",
#                                       full.names = TRUE,
#                                       pattern = "chance_perm_n_[0-9]{1,4}_97.tsv")

perm_chances <- lapply(perm_chances_file_names, fread) 

for (i in seq_along(perm_chances)) {
  perm_chances[[i]][, perm_n := i]
}

perm_chances <- rbindlist(perm_chances)

fwrite(perm_chances, file = xargs$out, sep = "\t")

