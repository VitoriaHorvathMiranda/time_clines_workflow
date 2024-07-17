#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments 
parser <- ArgumentParser(description= "makes a list of positions with higth fst snps between eu and afr")
parser$add_argument('--painel', '-painel', help= 'ancestry painel from raw_ancestral_painel_prep.R script')
parser$add_argument('--snpfst', '-fst', help= 'snp fst with eu and afr from grenedalf')
parser$add_argument('--output', '-o', help= 'a list of positions')
xargs<- parser$parse_args()

anc_painel <- fread(xargs$painel)
snp_fst <- fread(xargs$snpfst)

# anc_painel <- fread("/dados/time_clines/analysis/ancestry/raw_ancestry_painel.tsv")
# snp_fst <- fread("/dados/time_clines/analysis/ancestry/FST_ancestry/EU_AFR_SNPs_fst.csv")

setnames(snp_fst, 
         c("chrom","start","EU:AFR"),
         c("CHROM","POS","EU_AFR"))

snp_fst <- snp_fst[, !c("end", "snps"), with = FALSE]

anc_painel[, total_eu := ref_allele_count_EU + alt_allele_count_EU]
anc_painel[, total_afr := ref_allele_count_AFR + alt_allele_count_AFR]
anc_painel <- anc_painel[total_eu >= 15 & total_afr >= 15]

anc_painel[, freq_eu := alt_allele_count_EU/total_eu]
anc_painel[, freq_afr := alt_allele_count_AFR/total_afr]
anc_painel <- anc_painel[!(freq_eu == 0 & freq_afr == 0)]
anc_painel[, freq_diff := abs(freq_afr-freq_eu)]


fst_merged <- merge(snp_fst, anc_painel, all = FALSE)
setorder(fst_merged, -EU_AFR)

top_snp_fst <- fst_merged[1:ceiling(0.001 * .N)]

top_snp_fst[, position2 := paste0(CHROM, ":", POS)]

top_positions <- top_snp_fst[, .(position2)]


fwrite(top_positions, file = xargs$output, col.names = FALSE)

