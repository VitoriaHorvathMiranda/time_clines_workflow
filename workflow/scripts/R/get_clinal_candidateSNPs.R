#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "gets candidate snps for GOWINDA")
parser$add_argument('--effect97', '-e97', help = 'output from get_effects.R')
parser$add_argument('--effect0910', '-e0910', help = 'output from get_effects.R')
parser$add_argument('--FDR', '-fdr', help = 'fdr to filter candidate snps')
parser$add_argument('--oldclinal', '-clinal97', help = 'output with candidade snps ony clinal in 97')
parser$add_argument('--newclinal', '-clinal0910', help = 'output with candidade snps ony clinal in 0910')
xargs<- parser$parse_args()

#--------------------------------------------------------------------------------

q_97 <- fread(xargs$effect97)
q_0910 <- fread(xargs$effect0910)
fdr <- xargs$FDR

setcolorder(q_97, c("chrom", "POS"))
setcolorder(q_0910, c("chrom", "POS"))

new_names_97 <- colnames(q_97)[3:6] |> paste0("_97")
new_names_0910 <- colnames(q_0910)[3:6] |> paste0("_0910")

setnames(q_97, colnames(q_97)[3:6], new_names_97)
setnames(q_0910, colnames(q_0910)[3:6], new_names_0910)


all_tests <- merge.data.table(q_97, q_0910, by = c("chrom", "POS", "region", "effect", "ANN", "ALT"))

decrease <- all_tests[qvalue_97 <=fdr & qvalue_0910 >= fdr]
increase <- all_tests[qvalue_97 >=fdr & qvalue_0910 <= fdr]


fwrite(decrease, file = xargs$oldclinal, sep = "\t", col.names = FALSE)
fwrite(increase, file = xargs$newclinal, sep = "\t", col.names = FALSE)

