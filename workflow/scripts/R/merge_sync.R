#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments 
parser <- ArgumentParser(description= "")
parser$add_argument('--painel', '-painel', help= 'path to sync files with ancestral pops')
parser$add_argument('--sync', '-sync', help= ' original sync file')
parser$add_argument('--output', '-o', help= 'sync file with all pops')
xargs<- parser$parse_args()

# sync_samples <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync")

# anc_files <- list.files(path = "/dados/time_clines/analysis/ancestry",
#                         pattern = "ancestral_chrom_(2|3|R|L|X){1,2}.sync",
#                         full.names = TRUE)

anc_files <- list.files(path = xargs$painel,
                        pattern = "ancestral_chrom_(2|3|R|L|X){1,2}.sync",
                        full.names = TRUE)

sync_anc <- lapply(anc_files, fread)
sync_anc <- rbindlist(sync_anc)

sync_samples <- fread(xargs$sync)

setnames(sync_samples, paste0("V", 1:3), colnames(sync_anc)[1:3])

all_sync <- merge.data.table(sync_samples, sync_anc,
                             by = colnames(sync_anc)[1:3],
                             all = TRUE)

fwrite(all_sync, xargs$output, sep = "\t", 
       col.names = TRUE, na = "NA", quote = FALSE)
