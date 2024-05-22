#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

#parse arguments 
parser <- ArgumentParser(description= "Makes pca plots, all outputs are .jpeg")
parser$add_argument('--thetaPi', '-pi', help = "thetaPI from grenegalf single")
#parser$add_argument('--vcf', '-vcf', help= 'The vcf')
parser$add_argument('--plot', '-plot', help= 'Plot of theta pi per 30kb windows')
xargs<- parser$parse_args()

stats <- fread(xargs$thetaPi)

# small_vcf <- fread(xargs$vcf,
#                    skip = "#CHROM",
#                    nrows = 1)

stats <- stats[, !(colnames(stats) %like% "snp_count"), with = FALSE]
stats <- stats[, !(colnames(stats) %like% "theta_pi_rel"), with = FALSE]
stats <- stats[, !(colnames(stats) %like% "coverage_fraction"), with = FALSE]

sample_name <- colnames(stats)[4:length(stats)] |> str_split_i(pattern = "\\.", 
                                                               i = 1)
setnames(stats, colnames(stats)[4:length(stats)], sample_name)
stats <- stats[, !"end", with = FALSE]

#creat windowns
#30kb windowns
windows <- seq(1, 32079331, by = 30000)
stats[, window := cut(start, windows)]

#samples(pops) to row
stats_tidy <- melt.data.table(data = stats,
                                   measure = 3:21,
                                   value.name = "thetaPI",
                                   variable.name = "pop")

#computes mean thetaPI per window
mean_thetaPI_window <- 
  stats_tidy[, .(mean_pi = mean(thetaPI)), 
                  keyby = c("chrom", "pop", "window")]

#extracts smallest position of each window
mean_thetaPI_window[,min_pos := tstrsplit(window, ",", keep = 1)][
  , min_pos := sub(pattern = "\\(",
                   replacement = "",
                   x = min_pos)]


#computes mid position of each window
mean_thetaPI_window[, min_pos := as.double(min_pos)][
  , mid_pos := (min_pos+15000)
]

mean_thetaPI_window[, year := fcase(pop %like% "10" | pop %like% "09", "2009/2010",
                                    pop %like% "97", "1997",
                                    pop %like% "17", "2017",
                                    pop %like% "22" | pop %like% "23", "2022/2023")]


mean_thetaPI_window[, pop := as.character(pop)]
mean_thetaPI_window[, pop := ifelse(pop == "HFL97_new", "ESC97_new", pop)]

THETA_PI_97 <- 
mean_thetaPI_window[year == "1997"] |>
  ggplot() +
  geom_line(aes(x = mid_pos, y = mean_pi,
                color = pop, group = pop),
            #show.legend = FALSE,
            linewidth = 0.3) +
  facet_wrap(vars(chrom), scales = "free_x") +
  scale_color_brewer(palette = "Set1") +
  theme_light() + 
  labs(y = "θπ", x = "Position",
       title = "1997", color = "Population") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", size = 12),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size = 12),
        title = element_text(size = 14)) 

THETA_PI_0910 <- 
  mean_thetaPI_window[year == "2009/2010"] |>
  ggplot() +
  geom_line(aes(x = mid_pos, y = mean_pi,
                color = pop, group = pop),
            #show.legend = FALSE,
            linewidth = 0.3) +
  facet_wrap(vars(chrom), scales = "free_x") +
  scale_color_brewer(palette = "Set1") +
  theme_light() + 
  labs(y = "θπ", x = "Position",
       title = "2009/2010", color = "Population") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", size = 12),
        legend.title = element_text(size = 12),
        title = element_text(size = 14)) 


jpeg(filename = xargs$plot,
     width = 17,
     height = 20,
     units = "cm",
     res = 1200)

THETA_PI_97/THETA_PI_0910

dev.off()
