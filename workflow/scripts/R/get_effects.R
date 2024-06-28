#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparse)


#parse arguments
parser <- ArgumentParser(description= "compute significative SNPs for a given FDR")
parser$add_argument('--vcf', '-vcf',
                    help = 'annotated vcf')
parser$add_argument('--effects', '-ef',
                    help = 'table with sequence ontology and effects')
parser$add_argument('--qValue', '-qv',
                    help = 'table with q-values for each tested snp')
parser$add_argument('--output', '-o',
                    help = "table (tsv) with snps positions, q-values
                    and it's annotation")
xargs<- parser$parse_args()


effects <- fread(xargs$effects)
VCF <- fread(xargs$vcf,
             select = c(1,2,5,8))

q_values <- fread(xargs$qValue)

setnames(VCF, "#CHROM", "chrom")

# filter for biallelic snps
VCF <- VCF[!(ALT %like% ","),]

VCF[, ANN := gsub(".*ANN=[A-Z].(.+?)\\|.*", "\\1", INFO)]
VCF[, impact := str_split_i(INFO, pattern = "\\|", i = 3)]
VCF <- VCF[, !"INFO", with = FALSE]

VCF_effects <- merge.data.table(VCF, effects, by = "ANN")

q_value_eff <- merge.data.table(q_values, VCF_effects, by = c("POS", "chrom"))


fwrite(q_value_eff, file = xargs$output, sep = "\t")



