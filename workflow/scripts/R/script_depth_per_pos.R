#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(argparse)

#parse arguments
parser <- ArgumentParser(description= "gets mean depth per 80kb windows for samples in the pipeline")
parser$add_argument('--PreDepthFile', '-predf', help= 'before downsample samtools depth file')
parser$add_argument('--PosDepthFile', '-posdf', help= 'after downsample samtools depth file')
parser$add_argument('--chrom', '-c', help= 'chromosome')
parser$add_argument('--output', '-o', help= 'tsv, with: chrom, sample, window, mean_depth before and after downsample')
xargs<- parser$parse_args()

# chrom <- "X"
# depth_pre_path <- "/dados/time_clines/data/seqs/processed/qltctrl/pre_downsample/samtools_depth_stats_include_unmerged_flag.tsv"
# depth_pos_path <- "/dados/time_clines/data/seqs/processed/qltctrl/pos_downsample/samtools_depth_stats_include_unmerged_flag.tsv"

chrom <- xargs$chrom
depth_pre_path <- xargs$PreDepthFile
depth_pos_path <- xargs$PosDepthFile

depth_pre_header <- fread(depth_pre_path,
                          nrows = 1)
depth_pre <- fread(cmd = paste("grep", chrom , depth_pre_path))

depth_pos_header <- fread(depth_pos_path,
                          nrow = 1)
depth_pos <- fread(cmd = paste("grep", chrom , depth_pos_path)) 

metadata <- fread("/dados/time_clines/data/meta/seq_metadata.tsv", fill = TRUE)

#fix variables names
sample_ids_pos <- colnames(depth_pos_header)[3:length(depth_pos_header)] |> 
  str_split_i(pattern = "/", i = 7) |> 
  str_split_i(pattern = "_", i = 1)

sample_ids_pre <- colnames(depth_pre_header)[3:length(depth_pre_header)] |> 
  str_split_i(pattern = "/", i = 7) |> 
  str_split_i(pattern = "_", i = 1)

setnames(depth_pre, colnames(depth_pre)[3:length(depth_pre)], sample_ids_pre)
setnames(depth_pos, colnames(depth_pos)[3:length(depth_pos)], sample_ids_pos)
setnames(depth_pre, c("V1", "V2"), c("CHROM", "POS"))
setnames(depth_pos, c("V1", "V2"), c("CHROM", "POS"))

#creat windowns
#80kb windowns
windows <- seq(1, 32079331, by = 80000)
depth_pre[, window := cut(POS, windows)]
depth_pos[, window := cut(POS, windows)]

#samples(pops) to row
depths_pre_tidy <- melt.data.table(data = depth_pre,
                                   measure = 3:(length(depth_pre)-1),
                                   value.name = "depth",
                                   variable.name = "seq_label")

depths_pos_tidy <- melt.data.table(data = depth_pos,
                                   measure = 3:(length(depth_pre)-1),
                                   value.name = "depth",
                                   variable.name = "seq_label")

#computes mean depth per window

mean_depth_window_pre <- 
  depths_pre_tidy[, .(mean_depth = mean(depth, na.rm = TRUE)), 
                  keyby = c("CHROM", "seq_label", "window")]

mean_depth_window_pos <- 
  depths_pos_tidy[, .(mean_depth = mean(depth, na.rm = TRUE)), 
                  keyby = c("CHROM", "seq_label", "window")]

rm(depths_pos_tidy, depths_pre_tidy, depth_pos, depth_pre)


setnames(mean_depth_window_pre, "mean_depth", "mean_depth_pre")
setnames(mean_depth_window_pos, "mean_depth", "mean_depth_pos")

all_depth_per_window <- 
  mean_depth_window_pre[mean_depth_window_pos,
                      on = c("CHROM", "seq_label", "window")]

fwrite(x = all_depth_per_window, file = xargs$output)

