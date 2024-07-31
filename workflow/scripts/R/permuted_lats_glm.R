#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(combinat)
library(qvalue)
#library(parallel)

#parse arguments 
parser <- ArgumentParser(description= "permute latitudes perform glms and get p-values")
parser$add_argument('--cyear', '-cyear', help = "clinal year to be permuted")
parser$add_argument('--permN', '-permN', type="integer", help= 'permutation number')
parser$add_argument('--meta', '-meta', help= 'metadata table')
parser$add_argument('--freqs', '-freqs', help= 'table with snps freqs of the populations in that specific clinal year')
parser$add_argument('--effects', '-e', help='effect of each snp')
parser$add_argument('--out', '-out', help= 'p-values for each snp for that permuted latitude')
xargs<- parser$parse_args()

year <- xargs$cyear
perm_n <- xargs$permN

I_info <-  fread(xargs$meta,
                 select = c("population", "latitude"))
# I_info <-  fread("/dados/time_clines/data/meta/seq_metadata.tsv",
#                select = c("population", "latitude"))

#data_prep ---------------------------------------------------------------------

freqs_filt <- fread(xargs$freqs,
                    drop = c(3,6,7,8))

# freqs_filt <- fread("/dados/time_clines/data/seqs/calls/TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv",
#                     drop = c(3,6,7,8))

setnames(freqs_filt, colnames(freqs_filt), c("population", "latitude",
                                               "CHROM", "POS","freq","NE"))


#finds SNPs that aren't present in that period or that are present in only
#one  or two population
freqs_filt[, freq_status := 1*(freq!= 0)]
no_0freq <- freqs_filt[, .(freq_status_sum = sum(freq_status)),
                        by = c("POS", "CHROM")][freq_status_sum <= 2,]

#filter those SNPs out
freqs_filt <- freqs_filt[!no_0freq, on = c("POS", "CHROM")]
freqs_filt[, freq_status := NULL]

#data.table function for nesting: ---------------------------------------------
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}


#permute latitudes -------------------------------------------------------------
#get relevant columns 
I_info_filt <-
  I_info[population %like% year & !(population %in% c("ESC97", "HFL97_new","SNC10","dlFL10"))]


lat_perm <- permn(I_info_filt$latitude)[perm_n]
perm_meta <- I_info_filt[, latitude := lat_perm[[1]]]


permuted <- freqs_filt[perm_meta, on = "population"]

nested_snps <- permuted[,group_nest_dt(.SD, CHROM, POS)]

get_pvalue <- function(dt) {
  fit <- glm(freq ~ i.latitude, weights = NE, 
             data = dt,
             family = binomial()) |>
    summary()
  fit$coefficients[2, 4]
}

nested_snps[, p_value := sapply(data, get_pvalue)]
nested_snps[, data := NULL]


#get qvalue --------------------------------------------------------------------

p_values <- nested_snps$p_value
qobj <- qvalue(p = p_values)
nested_snps[, qvalue := qobj$qvalues]


#get effects -------------------------------------------------------------------

snp_effects <- fread(xargs$effects,
                     select = c(2,1,10))

# snp_effects <- fread("/dados/time_clines/analysis/time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv",
#                      select = c(2,1,10))
setnames(snp_effects, colnames(snp_effects), c("CHROM","POS", "effect"))


pvalues_perm_ann <- merge.data.table(nested_snps, snp_effects,
                                     by = c("CHROM","POS"), all = FALSE)


pvalues_perm_ann[, clinal_10 := ifelse(qvalue<=0.1, "clinal", "not_clinal")]

clinal_effect_chance <- 
  pvalues_perm_ann[, .(n = nrow(.SD)), by = c("effect", "clinal_10")] |>
  dcast.data.table(formula = effect ~ clinal_10,
                   value.var = "n")

clinal_effect_chance[, chance := clinal/not_clinal]

fwrite(x = clinal_effect_chance,
       file = xargs$out,
       sep = "\t")
