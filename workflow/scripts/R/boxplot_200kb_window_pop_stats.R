library(data.table)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

#autossomos -------------------------------------------------
window_sizes <- fread("/dados/time_clines/analysis/pop_stats/truewindows-200000-200000.txt")
PI <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.pi")
theta <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.th")
#Dtajima <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.D")
window_sizes_X <- fread("/dados/time_clines/analysis/pop_stats/truewindows_X-200000-200000.txt")
PI_X <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.pi")
theta_X <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.th")

small_vcf <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf",
                   skip = "#CHROM",
                   nrows = 1,
                   drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                            "INFO","FORMAT"))
PI <- rbindlist(list(PI, PI_X))
theta <- rbindlist(list(theta, theta_X))
window_sizes <- rbindlist(list(window_sizes, window_sizes_X))

pop_info <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
                  select = c("population", "collection_year", "latitude"))

new_colnames <- c("CHROM", "window", colnames(small_vcf))

lapply(list(PI, theta, window_sizes),#, Dtajima, ),
       function(x) setnames(x, colnames(x), new_colnames))

melted_stats <- lapply(list(PI, theta, window_sizes),
                       function(x) melt.data.table(x, id.vars = c(1,2),
                                                   measure.vars = 3:length(x),
                                                   variable.name = "population",
                                                   value.name = "stat"))
melted_stats <- 
lapply(melted_stats, 
       function(x) x[pop_info, on = "population", nomatch = 0][
         ,plot_collection_year :=
           fcase(collection_year == "2009" | collection_year == "2010", "2009/2010", 
                 collection_year == "2022" | collection_year == "2023", "2022/2023",
                 rep(TRUE, .N), as.character(collection_year))])

lapply(melted_stats, function(x) x[, state := as.character(fcase(population == "CMD97B", "MD",
                                                                  population %like% "HFL" | population %like% "dlFL", "sFL",
                                                                  rep(TRUE, .N), as.character(str_sub_all(population, start = -4, -3))))])

lapply(melted_stats, function(x) x[,state := fct_reorder(as.factor(state), latitude)])
lapply(melted_stats, function(x) x[, population := fct_reorder(as.factor(population), latitude)])

PLOTS <- 
lapply(melted_stats, function(x) x |>
         ggplot(aes(y= stat, x= population, fill = plot_collection_year)) +
         geom_boxplot() +
         #geom_histogram(binwidth = 0.002) +
         facet_grid(#cols = vars(state),
                    cols = vars(plot_collection_year),
                    scales = "free", space='free') +
         labs(fill = "Colletion Year") +
         scale_fill_brewer(palette = "Dark2") +
         theme_light() +
         theme(strip.text = element_text(color = "black"),
               strip.background = element_rect(fill = "gray95"),
               legend.position = "none",
               axis.text.x = element_text(size = 8, angle = 35, hjust = 0.8)))

PLOTS[[1]]


#SANITY CHECK
melted_stats[[3]][CHROM == "X"] |>
  ggplot(aes(y= stat, x= population, fill = plot_collection_year)) +
  geom_boxplot() +
  #geom_histogram(binwidth = 0.002) +
  facet_grid(cols = vars(state), scales = "free", space='free') +
  labs(fill = "Colletion Year") +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "gray95"),
        legend.position = "none",
        axis.text.x = element_text(size = 8, angle = 35, hjust = 0.8))


stats <- c("θπ", "θw", "Tajima's D")

jpeg(filename = "boxplot_200kb_window_pop_stats.jpeg",
     width = 30,
     height = 20,
     units = "cm",
     res = 1200)

(PLOTS[[1]] + labs(y = "θπ", x = "") + 
    theme(legend.position = "right")) /
(PLOTS[[2]] + labs(y = "θw", x = "")) /
(PLOTS[[3]] + labs(y = "Tajima's D"))

dev.off()

#CHROMOSSOMO X -----------------------------------------------------------------
window_sizes <- fread("/dados/time_clines/analysis/pop_stats/truewindows_X-200000-200000.txt")
PI <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.pi")
theta <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.th")
Dtajima <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.D")

small_vcf <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf",
                   skip = "#CHROM",
                   nrows = 1,
                   drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                            "INFO","FORMAT"))

pop_info <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
                  select = c("population", "collection_year", "latitude"))

new_colnames <- c("CHROM", "window", colnames(small_vcf))

lapply(list(PI, theta, Dtajima, window_sizes),
       function(x) setnames(x, colnames(x), new_colnames))

summary(window_sizes$HFL97)

melted_stats <- lapply(list(PI, theta, Dtajima),
                       function(x) melt.data.table(x, id.vars = c(1,2),
                                                   measure.vars = 3:length(x),
                                                   variable.name = "population",
                                                   value.name = "stat"))
melted_stats <- 
  lapply(melted_stats, 
         function(x) x[pop_info, on = "population", nomatch = 0][
           ,plot_collection_year :=
             fcase(collection_year == "2009" | collection_year == "2010", "2009/2010", 
                   collection_year == "2022" | collection_year == "2023", "2022/2023",
                   rep(TRUE, .N), as.character(collection_year))])

lapply(melted_stats, function(x) x[, state := as.character(fcase(population == "CMD97B", "MD",
                                                                 population %like% "HFL" | population %like% "dlFL", "sFL",
                                                                 rep(TRUE, .N), as.character(str_sub_all(population, start = -4, -3))))])
lapply(melted_stats, function(x) x[,state := fct_reorder(as.factor(state), latitude)])
lapply(melted_stats, function(x) x[, population := fct_reorder(as.factor(population), collection_year)])

PLOTS <- 
  lapply(melted_stats, function(x) x |>
           ggplot(aes(y= stat, x= population, fill = plot_collection_year)) +
           geom_boxplot() +
           #geom_histogram(binwidth = 0.002) +
           facet_grid(cols = vars(state), scales = "free", space='free') +
           labs(fill = "Colletion Year") +
           scale_fill_brewer(palette = "Dark2") +
           theme_light() +
           theme(strip.text = element_text(color = "black"),
                 strip.background = element_rect(fill = "gray95"),
                 legend.position = "none",
                 axis.text.x = element_text(size = 8, angle = 35, hjust = 0.8)))

stats <- c("θπ", "θw")

jpeg(filename = "boxplot_200kb_window_pop_stats.jpeg",
     width = 35/1.7,
     height = 20/1.7,
     units = "cm",
     res = 1200)

(PLOTS[[1]] + labs(y = "θπ", x = "") + 
  theme(legend.position = "right")) /
  (PLOTS[[2]] + labs(y = "θw", x = ""))
 #+ theme(legend.position = "right")) /
#(PLOTS[[1]] + labs(y = "θπ", x = "")) / (PLOTS[[2]] + labs(y = "θw", x = "")) #/
  #(PLOTS[[3]] + labs(y = "Tajima's D"))

dev.off()

