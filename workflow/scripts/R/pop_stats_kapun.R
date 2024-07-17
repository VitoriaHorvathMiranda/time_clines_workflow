#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(ggthemes)

#parse arguments 
parser <- ArgumentParser(description= "Makes pi, w and tajima's D plots, all outputs are .jpeg")
parser$add_argument('--vcf', '-vcf', help= 'The vcf')
parser$add_argument('--meta', '-meta', help = "metadata table")
parser$add_argument('--stats', '-stats', help = "stats path and prefix")
parser$add_argument('--plots', '-plots', help= 'plots path and prexix')
xargs<- parser$parse_args()


small_vcf <- fread(xargs$vcf,
                   skip = "#CHROM",
                   nrows = 1,
                   drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                            "INFO","FORMAT"))


# small_vcf <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf",
#                   skip = "#CHROM",
#                   nrows = 1,
#                   drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
#                            "INFO","FORMAT"))

meta <- fread(xargs$meta,
              select = c("population", "latitude", "collection_year"))

# meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
#               select = c("population", "latitude", "collection_year"))

arg_stats <- xargs$stats
#arg_stats <- "/dados/time_clines/analysis/pop_stats/pop_stats_kapun_auto_200000_200000"
path <- str_split_i(arg_stats, "/pop_stats_", 1)
pattern_auto <- str_split_i(arg_stats, "/", 6)
pattern_X <-  str_replace(pattern_auto, "auto", "X")

new_col_names <- c("CHROM", "POS", colnames(small_vcf))
stats_files <- list.files(path,
                          pattern = pattern_auto,
                          full.names = T)

stats_files_X <- list.files(path,
                          pattern = pattern_X,
                          full.names = T)


stat_names <- c("D", "pi", "w")
  
stats <- lapply(stats_files, fread)
stats_X <- lapply(stats_files_X, fread)

stats <- lapply(stats, function(x) 
  setnames(x, colnames(x), new_col_names))

stats_X <- lapply(stats_X, function(x) 
  setnames(x, colnames(x), new_col_names))

all_stats <- vector("list", length(stat_names))
for (i in seq_along(stats)) {
  all_stats[[i]] <- rbindlist(list(stats[[i]], stats_X[[i]]))
}


melted_stats <- lapply(all_stats, function(x)
  melt.data.table(x, measure.vars = 3:length(x),
                  value.name = "stat",
                  variable.name = "pop"))

for (i in seq_along(stat_names)) {
  setnames(melted_stats[[i]], "stat", stat_names[[i]])
}

merged_data <- Reduce(merge.data.table, melted_stats)


merged_data <- merge.data.table(merged_data, meta, by.x = "pop", by.y = "population")

merged_data[, state := as.character(fcase(pop == "CMD97B", "MD",
                                          pop %like% "HFL" | pop %like% "dlFL", "sFL",
                                          rep(TRUE, .N), as.character(str_sub_all(pop, start = -4, -3))))]

merged_data[, state := fct_reorder(as.factor(state),latitude)]

merged_data[, plot_year := fcase(collection_year == "1997", "1997",
                                 collection_year == "2009" | collection_year == "2010", "2009/2010",
                                 collection_year == "2022" | collection_year == "2023" | collection_year == "2017", "2017/2022/2023")]

PI <- 
merged_data |> 
  ggplot() +
  geom_line(aes(x = POS, y = pi,
                color = state, group = pop),
            #show.legend = FALSE,
            linewidth = 0.5) +
  facet_grid(cols = vars(CHROM), rows = vars(plot_year), scales = "free_x") +
  scale_color_tableau() +
  theme_light() + 
  labs(y = "θπ", x = "Position",
       title = "", color = "US State") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", size = 12),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size = 12),
        title = element_text(size = 14)) 

W <- 
merged_data |> 
  ggplot() +
  geom_line(aes(x = POS, y = w,
                color = state, group = pop),
            #show.legend = FALSE,
            linewidth = 0.5) +
  facet_grid(cols = vars(CHROM), rows = vars(plot_year), scales = "free_x") +
  scale_color_tableau() +
  theme_light() + 
  labs(y = "θW", x = "Position",
       title = "", color = "US State") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", size = 12),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size = 12),
        title = element_text(size = 14)) 
D <- 
merged_data |> 
  ggplot() +
  geom_line(aes(x = POS, y = D,
                color = state, group = pop),
            #show.legend = FALSE,
            linewidth = 0.3) +
  facet_grid(cols = vars(CHROM), rows = vars(plot_year), scales = "free_x") +
  scale_color_tableau() +
  #coord_cartesian(ylim = c(-1, 0.7)) +
  theme_light() + 
  labs(y = "Tajima's D", x = "Position",
       title = "", color = "US State") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", size = 12),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size = 12),
        title = element_text(size = 14)) 

plot_path <- xargs$plots
#plot_path <- c("../results/pop_stats_")

pi <- paste0(plot_path, "pi.jpeg")
w <- paste0(plot_path, "w.jpeg")
d <- paste0(plot_path, "D.jpeg")

jpeg(filename = pi,
     width = 23,
     height = 17,
     units = "cm",
     res = 1200)

PI

dev.off()

jpeg(filename = w,
     width = 23,
     height = 17,
     units = "cm",
     res = 1200)

W

dev.off()

jpeg(filename = d,
     width = 23,
     height = 17,
     units = "cm",
     res = 1200)

D

dev.off()

