#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "gets latitude and effective Nº of chrom")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--vcf', '-vcf', help = 'vcf file')
parser$add_argument('--outputAuto', '-oa', help= 'csv with population and mean number of chromosomes (for autosomes)')
parser$add_argument('--outputX', '-ox', help= 'csv with population and mean number of chromosomes (for X chrom)')
xargs<- parser$parse_args()

#precisa colocar a ordem certa porque o arquivo sync não tem nome das colunas
small_vcf <- fread(xargs$vcf,
                   skip = "#CHROM",
                   nrows = 1)

pop_order <- colnames(small_vcf)[10:28]

meta <- fread(xargs$metadata)

meta[, n_chrom_auto := fcase(flies_per_lin == 2, n_females*1.5,
                        flies_per_lin == 1, n_females*2)]

meta[, n_chrom_X := fcase(flies_per_lin == 2 & fly_sex == "female", n_females*1.5,
                          flies_per_lin == 1 & fly_sex == "female", n_females*2,
                          flies_per_lin == 2 & fly_sex == "male", n_females*0.5,
                          flies_per_lin == 1 & fly_sex == "male", n_females*1)]


pool_sizes_auto <- meta[, c("population", "n_chrom_auto"), with = FALSE] 
pool_sizes_X <- meta[, c("population", "n_chrom_X"), with = FALSE]

pool_sizes_auto[, population:= factor(population, levels = pop_order)]
pool_sizes_X[, population := factor(population, levels = pop_order)]

setorder(pool_sizes_auto, population)
setorder(pool_sizes_X, population)

pool_sizes_auto <- na.omit(pool_sizes_auto)
pool_sizes_X <- na.omit(pool_sizes_X)

pool_sizes_auto[, sample_name := paste0("PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.", c(1:19))]
pool_sizes_X[, sample_name := paste0("PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.", c(1:19))]

pool_sizes_auto <- data.table(pool_sizes_auto$sample_name, pool_sizes_auto$n_chrom_auto)
pool_sizes_X <- data.table(pool_sizes_X$sample_name, pool_sizes_X$n_chrom_X)


fwrite(pool_sizes_auto, file = xargs$outputAuto,
       col.names = FALSE, sep = ",")

fwrite(pool_sizes_X, file = xargs$outputX,
       col.names = FALSE, sep = ",")