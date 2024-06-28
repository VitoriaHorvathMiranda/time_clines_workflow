#!/usr/bin/env Rscript

#this script takes the vcf with african and european populations
#and transforms it in a pre_AHH input

library(argparse)
library(tidyverse)
library(data.table)

#parse arguments 
parser <- ArgumentParser(description= "makes a raw ref ancestry painel,
                         columns: chrom, pos(snp), alt, ref,
                         ref_count on european pops, alt_count on european pops,
                         ref_count on African pops, alt_count on African pops")
parser$add_argument('--vcf', '-vcf', help = "vcf with ref ancestral haploid genomes")
parser$add_argument('--output', '-out', help= 'raw ancestry painel')
xargs<- parser$parse_args()

vcf_painel <- fread(xargs$vcf)

setnames(vcf_painel, "#CHROM", "CHROM")

# os cromosomos estavam com esses nomes porque originalmente eu pretendia usar o \n
#plick para calcular a ancestralidade global dos meus genomas haploides
vcf_painel[, CHROM := fcase(CHROM == "1", "2L",
                            CHROM == "2", "2R",
                            CHROM == "3", "3L",
                            CHROM == "4", "3R",
                            rep(TRUE, .N), as.character(CHROM))]


vcf_painel_tidy <- melt.data.table(vcf_painel, id.vars = 1:5,
                                   variable.name = "sample",
                                   value.name = "genotype")

#dt[, c("PX", "PY") := tstrsplit(PREFIX, "_", fixed=TRUE)]
vcf_painel_tidy[, c("GT","AD","DP","GQ","PL") := 
                  tstrsplit(genotype, ":", fixed = TRUE)]

#GT = genotype
#DP = Depth
vcf_painel_tidy[, c("GT", "DP") := lapply(.SD, as.integer), .SDcols = c("GT","DP")]

vcf_painel_tidy <- vcf_painel_tidy[DP > 1]

vcf_painel_tidy[, ancestry := fcase(sample %like% "FR" | sample %like% "SU", "EU",
                                    default = "AFR")]

vcf_painel_tidy[, ref_gt := fcase(GT == 0, 1,
                                  GT == 1, 0)]

painel_count <- 
  vcf_painel_tidy[, .(alt_allele_count = sum(GT),
                      ref_allele_count = sum(ref_gt)),
                  by = c("CHROM", "POS", "REF", "ALT", "ancestry")]

painel_count_wide <- 
  dcast.data.table(painel_count,
                   CHROM + POS + REF + ALT ~ ancestry,
                   value.var = c("alt_allele_count", "ref_allele_count"))

#NAs tanto na coluna AFR quanto na coluna EU indicam que não temos
#infomação dessas posições nessas pops, então eu vou descartar
painel_count_wide <- painel_count_wide[!(is.na(ref_allele_count_AFR)) & !(is.na(ref_allele_count_EU))]

setcolorder(painel_count_wide,
            c("CHROM", "POS", "ALT", "REF",
              "ref_allele_count_EU", "alt_allele_count_EU",
              "ref_allele_count_AFR", "alt_allele_count_AFR"))

fwrite(painel_count_wide, file = xargs$output,
       sep = "\t")

