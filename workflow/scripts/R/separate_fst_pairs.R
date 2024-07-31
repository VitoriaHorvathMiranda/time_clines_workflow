#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)

#parse arguments 
parser <- ArgumentParser(description= "separates fst_cutoffs files into each pair")
parser$add_argument('--fst', '-fst', help = "fst file from outliers_*_FST_all_comp.R")
parser$add_argument('--output', '-o', help = "output path and prefix")
xargs<- parser$parse_args()


cutoff <- fread(xargs$fst)

pairs <- unique(cutoff$test) 

s_cutoffs <- lapply(pairs, function(x) cutoff[test == x])

full_o_names <- paste0(xargs$output, "_", pairs, ".tsv")

for (i in seq_along(pairs)) {
  fwrite(s_cutoffs[[i]], file = full_o_names[[i]],
         col.names = FALSE, sep = "\t")
  
}

