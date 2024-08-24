#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "")
parser$add_argument('--CoveragePath', '-cov', help= 'path to -samtools coverage- outputs')
parser$add_argument('--male_pools', '-mpools', help = 'name of the male pools "|" separated')
parser$add_argument('--output', '-o', help= 'tsv with fraq of reads to be removed')
xargs<- parser$parse_args()


#Stat_path <- xargs$StatsPath
Cov_path <- xargs$CoveragePath
#Cov_path <- "/dados/time_clines/data/seqs/processed/qltctrl/pre_downsample/samtools_coverage"

#male_pools <- "dlSC10_dlGA10_dlFL10_HFL97downto60mi"

male_pools <- xargs$male_pools
male_pools <- str_replace_all(male_pools, "_", "|")
male_pools <- gsub("([^|]+)", "\\1_L001", male_pools)
file_names <- paste0("(", male_pools,  ")_total_coverage.tsv")
coverage_files <- list.files(path = Cov_path,
                             pattern = file_names,
                             full.names = TRUE)

raw_coverages <- lapply(coverage_files, fread)
lapply(raw_coverages, setnames, "#rname", "chrom")

#gets the samples ids because files are out of order
sample_ids <- lapply(coverage_files, str_split_i, pattern = "/", i = 10) |>
  lapply(str_split_i, pattern = "_", i = 1)


#sample_names <- meta$population
coverages_pop_names <- vector("list", length(raw_coverages))
for (i in seq_along(sample_ids)) {
  coverages_pop_names[[i]] <- raw_coverages[[i]] %>%
    mutate(seq_label = sample_ids[[i]])
}

coverages_table <- rbindlist(coverages_pop_names)
downsample_fraq_x <- coverages_table[chrom == "X"][, downsample_fraq := 25/meandepth]


fwrite(downsample_fraq_x,
       xargs$output,
       sep = "\t")
