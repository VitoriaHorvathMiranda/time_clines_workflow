#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "make a bed file for neutral regions in short introns")
parser$add_argument('--gtf', '-gtf', help = 'a gtf file')
parser$add_argument('--introns', '-i', help = 'a bed file with all introns')
parser$add_argument('--output', '-o', help = "a bed file with neutral regions of short introns")
xargs<- parser$parse_args()


#get genes and gene info to find intron strand info
# full_genes <- fread(cmd = paste0("awk '$3 == \"gene\"' ", '/dados/time_clines/data/reference/dmel-6.55.gtf'),
#                     drop = c(2,6,8))
full_genes <- fread(cmd = paste0("awk '$3 == \"gene\"' ", xargs$gtf),
                    drop = c(2,6,8))
setnames(full_genes, colnames(full_genes), c("chrom", "type", "start", "end", "strand", "info"))
setkey(full_genes, chrom, start)

# Short intron sites -> bases 8–30 of introns 65 bp -> Parsch et al. 2010
# Short intron sites -> < 86 bp, 16 bp away from the intron start and 6 bp away from the intron end -> Lawrie 2013

#introns <- fread("/dados/time_clines/data/reference/dmel-6.32_introns.bed")
introns <- fread(xargs$introns)
setnames(introns, colnames(introns), c("chrom", "start", "end"))

introns[, i_length := end - start]

short_introns <- introns[i_length <= 86 & i_length > 22]
short_introns[, i_start := start]
short_introns[, i_end := end]

short_i_with_info <- short_introns[full_genes, on = .(chrom, start >= start, end <= end),
                       nomatch = 0L]

#filtra os introns que estão em mais de um gene
short_i_with_info[, filter := nrow(.SD), by = c("chrom", "i_start", "i_end")]
short_i_with_info <- short_i_with_info[!(filter > 1)]

#pega as posições neutras dos introns curtos
short_i_with_info[, neutral_start := ifelse(strand == '+',
                                            i_start + 16,
                                            i_start + 6)]

short_i_with_info[, neutral_end := ifelse(strand == '+',
                                          i_end - 6,
                                          i_end - 16)]

# Parsch et al. 2010:
# short_i_with_info[, neutral_start := ifelse(strand == '+',
#                                             i_start + 8,
#                                             i_end - 30)]
# 
# short_i_with_info[, neutral_end := ifelse(strand == '+',
#                                           i_start + 30,
#                                           i_end - 8)]

#sanity check
short_i_with_info[, .(check = neutral_end - neutral_start,
                      i_length)][i_length - check != 22]


#final table

final <- short_i_with_info[, .(chrom, neutral_start, neutral_end)]

fwrite(final, xargs$output, sep = "\t")




