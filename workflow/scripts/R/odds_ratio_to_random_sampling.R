#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(parallel)

#parse arguments 
parser <- ArgumentParser(description= "sample 100000 random snps and get the odds ratio of each genomic region for a critinal qvalue 10000 times")
parser$add_argument('--cutoff', '-cut', type="numeric", help = "q-value to consider a clinal snp")
parser$add_argument('--qvalues', '-qv', help= 'table with snps freqs of the populations in that specific clinal year')
parser$add_argument('--out', '-out', help= 'table with odds ratio for all 10000 trials')
xargs<- parser$parse_args()

critical_q <- xargs$cutoff
#critical_q <- 0.1
#read data

qvalues <- fread(xargs$qvalues,
                 select = c("chrom","POS","qvalue","effect"))

# qvalues <- fread("/dados/time_clines/analysis/time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv",
#                  select = c("chrom","POS","qvalue","effect"))



#calcula a chance de um snp clinal ser daquela determindade classe genômica
clinal_snps <- qvalues[qvalue <= critical_q]
n_snps <- clinal_snps[,.(n_r = nrow(.SD)), by = effect]
n_snps[, total_clinal := sum(n_snps$n_r)-n_r]
n_snps[, chance_r := n_r/total_clinal]

#pega 100000 snps aleatórios e calcula a chance dele daquela determindade classe genômica
get_controls <- 
function(dt, effects){
  control_snps <- dt[qvalue > 0.1][sample(.N, 100000)]
  
  n_control <- control_snps[,.(n = nrow(.SD)), by = "effect"]
  n_control[, total_n_class := sum(n_control$n)-n]
  n_control[, chance_control := n/total_n_class]
  return(n_control)
}

num_cores <- 30

#repete isso 10000 vezes
all_controls <- 
mclapply(1:10000, function(x) get_controls(dt = qvalues,
                                      effects = snp_effects),
         mc.cores = num_cores)

all_controls <- 
mclapply(all_controls, function(x) merge.data.table(n_snps, x, 
                                                  by = "effect"),
         mc.cores = num_cores)

all_controls <- 
rbindlist(all_controls)

all_controls[, odds_ratio := chance_r/chance_control]


fwrite(x = all_controls,
       file = xargs$out,
       sep = "\t")
