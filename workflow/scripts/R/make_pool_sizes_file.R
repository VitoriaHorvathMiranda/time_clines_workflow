#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "gets latitude and effective Nº of chrom")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--vcf', '-vcf', help = 'vcf file')
parser$add_argument('--syncString', '-string', help = 'name of the sync file')
parser$add_argument('--outputAuto', '-oa', help= 'csv with population and mean number of chromosomes (for autosomes)')
parser$add_argument('--outputX', '-ox', help= 'csv with population and mean number of chromosomes (for X chrom)')
parser$add_argument('--outputName', '-on', help= 'tsv with current name given by grenedalf and corresponding sample name')
xargs<- parser$parse_args()

#precisa colocar a ordem certa porque o arquivo sync não tem nome das colunas
small_vcf <- fread(xargs$vcf,
                   skip = "#CHROM",
                   nrows = 1)

pop_order <- colnames(small_vcf)[10:length(small_vcf)]

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

rename_samples <- 
data.table(sync_name = paste(xargs$syncString,
                             c(1:length(colnames(small_vcf)[10:length(small_vcf)])),
                             sep = "."),
           names = colnames(small_vcf)[10:length(small_vcf)])


fwrite(pool_sizes_auto, file = xargs$outputAuto,
       col.names = FALSE, sep = ",")

fwrite(pool_sizes_X, file = xargs$outputX,
       col.names = FALSE, sep = ",")

fwrite(rename_samples, file = xargs$outputName,
       col.names = FALSE, sep = "\t")


