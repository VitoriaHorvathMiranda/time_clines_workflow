#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "gets latitude and effective Nº of chrom")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--vcf', '-vcf', help = 'vcf file')
parser$add_argument('--output', '-o', help= 'csv with population and mean number of chromosomes (for autosomes)')
xargs<- parser$parse_args()

#precisa colocar a ordem certa porque o arquivo sync não tem nome das colunas
small_vcf <- fread(xargs$vcf,
                   skip = "#CHROM",
                   nrows = 1)

pop_order <- colnames(small_vcf)[10:28]

meta <- fread(xargs$metadata)

meta[, n_chrom := fcase(flies_per_lin == 2, n_females*1.5,
                        flies_per_lin == 1, n_females*2)]

pool_sizes <- meta[, c("population", "n_chrom"), with = FALSE] 

pool_sizes[, population:= factor(population, levels = pop_order)]

setorder(pool_sizes, population)

pool_sizes <- na.omit(pool_sizes)

pool_sizes[, sample_name := paste0("PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.", c(1:19))]

pool_sizes <- data.table(pool_sizes$sample_name, pool_sizes$n_chrom)

fwrite(pool_sizes, file = xargs$output,
       col.names = FALSE, sep = ",")
