#!/usr/bin/env Rscript
#como o meu vcf dos paineis é individual é ruim usar um script pronto para criar o .sync
#então esse script tenta fazer isso a partir do arquivo do painel
#sync files are A:T:C:G:N:del format

library(argparse)
library(data.table)
library(tidyverse)


#parse arguments 
parser <- ArgumentParser(description= "")
parser$add_argument('--painel', '-painel', help= 'ancestry painel from raw_ancestral_painel_prep.R script')
parser$add_argument('--output', '-o', help= 'sync file with ancetral pops')
xargs<- parser$parse_args()

#anc_painel <- fread("/dados/time_clines/analysis/ancestry/raw_ancestry_painel.tsv")
anc_painel <- fread(xargs$painel)

bases <- c("A", "T", "C", "G")

ref_cols <- c("ref_allele_count_EU", "ref_allele_count_AFR")
alt_cols <- c("alt_allele_count_EU", "alt_allele_count_AFR")

refEU_refAFR <- lapply(ref_cols, function(x) dcast.data.table(anc_painel, CHROM + POS + REF ~ REF,
                                              value.var = x,
                                              fill = 0))

altEU_altAFR <- lapply(alt_cols, function(x) dcast.data.table(anc_painel, CHROM + POS + REF ~ ALT,
                                                              value.var = x,
                                                              fill = 0))

# Function to sum columns of same ancestry
sum_tables <- function(dt1, dt2) {
  merged_dt <- merge(dt1, dt2, by = c("CHROM", "POS", "REF"))
  cols_to_sum <- setdiff(names(dt1), c("CHROM", "POS", "REF")) 
  result <- merged_dt[, (cols_to_sum) := lapply(cols_to_sum, function(col) get(paste0(col, ".x")) + get(paste0(col, ".y")))]
  result <- result[, .SD, .SDcols = c("CHROM", "POS", "REF", cols_to_sum)]
  return(result)
}

pre_sync_EU <- sum_tables(refEU_refAFR[[1]], altEU_altAFR[[1]])
pre_sync_AFR <- sum_tables(refEU_refAFR[[2]], altEU_altAFR[[2]])

pre_sync_AFR[, sync_AFR := do.call(paste, c(.SD, "0", "0", sep = ":")), .SDcols = bases]
pre_sync_EU[, sync_EU := do.call(paste, c(.SD, "0", "0", sep = ":")), .SDcols = bases]

sync <- merge.data.table(pre_sync_EU[, !(bases), with = FALSE],
                         pre_sync_AFR[, !(bases), with = FALSE])


fwrite(sync, file = xargs$output, sep = "\t")
