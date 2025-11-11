#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments 
parser <- ArgumentParser(description= "")
parser$add_argument('--greSYNC', '-gre', help= 'sync file positions missing in anc vcf')
parser$add_argument('--allSYNC', '-all', help= 'sync file with North America, europe and africa samples')
#parser$add_argument('--popsize', '-pop', help= 'file with pop sizes (only use it for pop names)')
parser$add_argument('--output', '-o', help= 'sync file with all pops and corrected snps (no NA)')
xargs<- parser$parse_args()


#gre_sync <- fread("/dados/time_clines/data/database_seqs/call/all_variants_including_pool_samples_sync.sync")
gre_sync <- fread(xargs$greSYNC)


#fix names 
names <- 
colnames(gre_sync)[4:length(gre_sync)] |>
  sapply(function(x) 
    ifelse(x %like% "SU" | x %like% "FR",
           paste0(str_split_1(x, "\\."), "_EU"),
           paste0(str_split_1(x, "\\."), "_AFR")))

setnames(gre_sync,
         colnames(gre_sync)[4:length(gre_sync)],
         names)
setnames(gre_sync, "#chr", "chr")

melted_sync <- 
melt.data.table(gre_sync, id.vars = c("chr", "pos", "ref"),
                measure.vars = 4:length(gre_sync),
                variable.name = "sample",
                value.name = "genotype")
#A:T:C:G:N:del
melted_sync[, c("A","T","C","G","N","del") :=
              tstrsplit(genotype, ":")]

melted_sync[, pop := ifelse(sample %like% "AFR", "AFR", "EU")]

melted_sync[, c("A","T","C","G","N","del") := lapply(.SD, as.integer),
            .SDcols = c("A","T","C","G","N","del")]

#filter positions that don't look like homozygotes 
#there shouldn't be heterozygotes positions because the samples are haploid
melted_sync[, all_bases := rowSums(.SD, na.rm = TRUE),
            .SDcols = c("A","T", "C", "G", "N", "del")]

melted_sync[, maxdepth := max(A,`T`,C,G,N,del),
            by = 1:nrow(melted_sync)]

melted_sync[, is_25_bigger := ifelse(all_bases >= maxdepth*1.25,
                                     TRUE, FALSE)]

melted_sync[is_25_bigger == TRUE, 
            c("A","T","C","G","N","del") := 0]

#corrects for reading errors
melted_sync[(str_count(genotype, "(1:|2:|3:|4:|5:|6:|7:|8:|9:|[1-9]{1,2}0:)") > 1) &
              is_25_bigger == FALSE,
            c("A","T","C","G","N","del") := 
              lapply(.SD, function(x)
                ifelse(x != maxdepth, 0, x)),
            .SDcols = c("A","T","C","G","N","del")]


#masks sites with depth too high
depth_filter <- 
melted_sync[, .(mean(all_bases)), by = sample][, max_depth_filter := V1*8]

melted_sync <- 
melted_sync[depth_filter, on = "sample"]

melted_sync[all_bases >= max_depth_filter, 
            c("A","T","C","G","N","del") := 0]

melted_sync[all_bases >= max_depth_filter]

#transforms nº of reads in genotype counts
melted_sync[, c("A","T","C","G","N","del") := 
              lapply(.SD, function(x)
                ifelse(x != 0, 1, 0)),
            .SDcols = c("A","T","C","G","N","del")]


pre_sync <- 
melted_sync[, .(A = sum(A, na.rm = TRUE),
               `T` = sum(`T`,na.rm = TRUE),
                C = sum(C, na.rm = TRUE),
                G = sum(G, na.rm = TRUE),
                N = sum(N, na.rm = TRUE),
                del = sum(del, na.rm = TRUE)), by = c("chr", "pos", "ref", "pop")]


pre_sync[, genotype := paste(A,`T`,C,G,N,del, sep = ":")]

#pre_sync[(str_count(genotype, "(1:|2:|3:|4:|5:|6:|7:|8:|9:)") > 1)]
#dt[dt[, Reduce(`|`, lapply(.SD, `==`, 4)),.SDcols = sel.col], ..sel.col]
# pre_sync[pre_sync[, Reduce( `&`, lapply(.SD, `==`, 0)),
#                   .SDcols = c("A","T","C","G","N","del")]]

sync <- 
dcast.data.table(pre_sync, 
                 formula = chr + pos + ref ~ pop,
                 value.var = "genotype")

#all_samples <- fread("/dados/time_clines/analysis/ancestry/all_samples_sample_fixed.sync")
all_samples <- fread(xargs$allSYNC)
setnames(all_samples, colnames(all_samples)[1:3], colnames(sync)[1:3])
setkey(all_samples, chr, pos, ref)

merged_sync <- 
merge.data.table(all_samples, sync, all = TRUE)

merged_sync[, sync_EU := ifelse(is.na(sync_EU),
                                EU,
                                sync_EU)]

merged_sync[, sync_AFR := ifelse(is.na(sync_AFR),
                                AFR,
                                sync_AFR)]

merged_sync <- merged_sync[, !c("AFR", "EU")]

merged_sync <- merged_sync[!(is.na(sync_EU))]

# sample_names <- 
# fread("/dados/time_clines/analysis/ancestry/FST_ancestry/pool_sizes_autosome_with_anc_painel.csv")

#sample_names <- 
#  fread(xargs$popsize)

#setnames(merged_sync, colnames(merged_sync)[4:length(merged_sync)],
#         sample_names$V1)


fwrite(merged_sync,
       file = xargs$output,
       sep = "\t", col.names = FALSE)

