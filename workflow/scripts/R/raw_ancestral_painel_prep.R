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



#format- GT:AD:DP:GQ:PL

# vcf <- fread("/dados/time_clines/data/database_seqs/call/AFR_EU_chrom_2L.biallelic.clean.vcf",
#              skip = "#CHROM",
#              drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT"))

vcf <- fread(xargs$vcf,
             skip = "#CHROM",
             drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT"))


setnames(vcf, "#CHROM", "CHROM")


vcf_tidy <- melt.data.table(vcf, id.vars = 1:4,
                                   variable.name = "sample",
                                   value.name = "genotype")

#dt[, c("PX", "PY") := tstrsplit(PREFIX, "_", fixed=TRUE)]
vcf_tidy[, c("GT","AD","DP","GQ","PL") := 
                  tstrsplit(genotype, ":", fixed = TRUE)]

#GT = genotype
#DP = Depth
vcf_tidy[, c("GT", "DP") := lapply(.SD, as.integer), .SDcols = c("GT","DP")]

vcf_tidy <- vcf_tidy[DP > 2]

vcf_tidy[, ancestry := fcase(sample %like% "FR" | sample %like% "SU", "EU",
                                    default = "AFR")]

vcf_tidy[, ref_gt := fcase(GT == 0, 1,
                           GT == 1, 0)]

painel_count <- 
  vcf_tidy[, .(alt_allele_count = sum(GT),
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

