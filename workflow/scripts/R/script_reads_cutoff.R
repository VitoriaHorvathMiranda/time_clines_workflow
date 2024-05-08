#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "")
parser$add_argument('--MetaData', '-meta', help= 'metadata_table')
parser$add_argument('--StatsPath', '-stats', help= 'path to -samtools idxstats- outputs')
parser$add_argument('--CoveragePath', '-cov', help= 'path to -samtools coverage- outputs')
parser$add_argument('--output', '-o', help= 'tsv with fraq of reads to be removed')
xargs<- parser$parse_args()

#get metadata
meta <- fread(file = xargs$MetaData,
              sep = "\t")

Stat_path <- xargs$StatsPath
Cov_path <- xargs$CoveragePath

# get total reads --------------------------------------------------------------
files_stats <- list.files(path = Stat_path,
                          pattern = ".tsv", full.names = TRUE)

#gets the samples ids because files are out of order
sample_ids <- lapply(files_stats, str_split_i, pattern = "/", i = 10) |>
  lapply(str_split_i, pattern = "_", i = 1)

raw_stats <- lapply(files_stats, fread)
lapply(raw_stats,
       setnames, 
       paste("V", 1:3, sep = ""),
       c("chrom", "n_pos", "r_reads"))

#sample_names <- meta$population
stats_pop_names <- vector("list", length(raw_stats))
for (i in seq_along(sample_ids)) {
  stats_pop_names[[i]] <- raw_stats[[i]] %>%
    mutate(seq_label = sample_ids[[i]])
}

all_stats <- rbindlist(stats_pop_names)

all_stats <- all_stats %>%
  group_by(seq_label) %>%
  summarise(total_pos = sum(n_pos),
            total_reads = sum(r_reads))


# get mean depth ---------------------------------------------------------------

coverage_files <- list.files(path = Cov_path,
                             pattern = "_coverage.tsv", full.names = TRUE)

raw_coverages <- lapply(coverage_files, fread)
lapply(raw_coverages, setnames, "#rname", "chrom")

#sample_names <- meta$population
coverages_pop_names <- vector("list", length(raw_coverages))
for (i in seq_along(sample_ids)) {
  coverages_pop_names[[i]] <- raw_coverages[[i]] %>%
    mutate(seq_label = sample_ids[[i]])
}

coverages_table <- rbindlist(coverages_pop_names)

coverages_table <- coverages_table %>%
  group_by(seq_label) %>%
  summarise(total_pos = sum(endpos),
            total_covbases = sum(covbases),
            total_coverage = sum(covbases)/sum(endpos),
            total_meandepth = weighted.mean(meandepth, endpos)) #%>%
  # mutate(year = case_when(population %like% "97" ~ "1997",
  #                         population %like% "09" | population %like% "10" ~ "2009/2010",
  #                         population %like% "17" ~ "2017"))

# join -------------------------------------------------------------------------

#company_info <- meta %>% select(population, company)

all_stats <- all_stats %>%
  left_join(coverages_table, by = c("seq_label", "total_pos")) %>%
  #left_join(company_info, by = "population") %>%
  mutate(reads_cut = (25*total_reads)/total_meandepth) %>% 
  mutate(frac_base_mean = 25/total_meandepth) %>%
  mutate(frac_reads = reads_cut/total_reads) 

all_stats <- all_stats |>
  mutate(seq_label = if_else(seq_label == "09", "9", seq_label)) |>
  left_join(meta, by = "seq_label")

fwrite(all_stats,
       xargs$output,
       sep = "\t")

