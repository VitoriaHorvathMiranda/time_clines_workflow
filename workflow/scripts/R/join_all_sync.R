#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments 
parser <- ArgumentParser(description= "transformns ancestral pops vcf into sync")
parser$add_argument('--vcf', '-vcf', help= 'vcf for colnames')
parser$add_argument('--ancSync', '-ancSync', help = "sync file with ancestral pops")
parser$add_argument('--mySync', '-mySync', help= 'sync file with sample pops')
parser$add_argument('--output', '-o', help= 'sync with all pops')
parser$add_argument('--outPath', '-opath', help= 'path to output files, it must end with /')
xargs<- parser$parse_args()

# anc_sync <- fread("/dados/time_clines/data/database_seqs/call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean.sync")
# 
# my_sync <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync")
# 
# pop_names <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.ann.vcf",
#                    drop = c("#CHROM", "POS",  "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
#                    skip = "#CHROM",
#                    nrows = 1)

anc_sync <- fread(xargs$ancSync)

my_sync <- fread(xargs$mySync)

pop_names <- fread(xargs$vcf,
                   drop = c("#CHROM", "POS",  "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
                   skip = "#CHROM",
                   nrows = 1)


setnames(anc_sync, colnames(anc_sync)[4:7], c("EU", "westAFR", "ZI", "EG"))
setnames(my_sync, 
         colnames(my_sync)[4:length(my_sync)],
         colnames(pop_names))


setkey(anc_sync, V1,V2,V3)
setkey(my_sync, V1,V2,V3)

pool_sync <- 
  merge.data.table(
    anc_sync, 
    my_sync,
    all = TRUE)


pool_sync[, colnames(pool_sync)[4:length(pool_sync)] := 
       lapply(.SD, function(x) 
         ifelse(is.na(x), 
                "0:0:0:0:0:0",
                x)),
     .SDcols = c(4:length(pool_sync))]

chrom <- c("2L", "2R", "3R", "3L", "X")

lapply(chrom, function(x) {
  fwrite(pool_sync[V1 == x], paste0(xargs$outPath, x, ".sync"),
         col.names = F,
         sep = "\t")
})

fwrite(pool_sync, xargs$output,
       col.names = F,
       sep = "\t")

