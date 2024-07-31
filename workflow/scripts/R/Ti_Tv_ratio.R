#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)

#parse arguments 
parser <- ArgumentParser(description= "makes Ti/Tv ratio table")
parser$add_argument('--freqs', '-freqs', help = "frequency per pop, position with alt and ref")
parser$add_argument('--out', '-o', help= 'Tv/Tv ratio table')
xargs<- parser$parse_args()


freqs <- fread(xargs$freqs)
#freqs <- fread("/dados/time_clines/data/seqs/calls/freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")

freqs[, freq := as.double(freq)]
freqs <- freqs[!is.na(freq)]

freqs[freq != 0, mutation_type := fcase(REF == "A" & ALT == "G", "Ti",
                                        REF == "C" & ALT == "T", "Ti",
                                        REF == "G" & ALT == "A", "Ti",
                                        REF == "T" & ALT == "C", "Ti",
                                        default = "Tv")]

dcasted_freqs <- 
  dcast.data.table(freqs,
                   formula = population + CHROM + POS + REF + ALT ~ mutation_type,
                   value.var = "freq")

freq_ti_tv <- 
  dcasted_freqs[, .(Ti = sum(Ti, na.rm = T), Tv = sum(Tv, na.rm = T)),
                by = population]

freq_ti_tv[, ratio := Ti/Tv]


fwrite(freq_ti_tv, file = xargs$out)
