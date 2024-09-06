library(data.table)
library(tidyverse)
library(argparse)

#parse arguments 
parser <- ArgumentParser(description= "computes pop stats (pi or theta) per chrom, based on 200kb windows")
parser$add_argument('--statA', '-statA', help = "stat from Poolgenvar (autosomes)")
parser$add_argument('--statX', '-statX', help= 'stat from Poolgenvar (X)')
parser$add_argument('--windowA', '-wA', help = "window sizes from Truewindows (autosomes)")
parser$add_argument('--windowX', '-wX', help= 'window sizes from Truewindows (X)')
parser$add_argument('--meta', '-meta', help= 'metadata table')
parser$add_argument('--vcf', '-vcf', help= 'vcf for header')
parser$add_argument('--output', '-out', help= 'global stat per chrom per pop')
parser$add_argument('--yearStat', '-y', help= 'global stat per chrom per year')
parser$add_argument('--statBar', '-bar', help= 'barplot of output, jpeg')
xargs<- parser$parse_args()



# pi_auto <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.pi")
# pi_X <- fread("/dados/time_clines/analysis/pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.pi")
# window_sizes <- fread("/dados/time_clines/analysis/pop_stats/truewindows-200000-200000.txt")
# window_sizes_X <- fread("/dados/time_clines/analysis/pop_stats/truewindows_X-200000-200000.txt")
# 
# pop_info <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
#                   select = c("population", "collection_year", "latitude"))
# 
# small_vcf <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf",
#                    skip = "#CHROM",
#                    nrows = 1,
#                    drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
#                             "INFO","FORMAT"))

##get data --------------------------------------
pi_auto <- fread(xargs$statA)
pi_X <- fread(xargs$statX)
window_sizes <- fread(xargs$windowA)
window_sizes_X <- fread(xargs$windowX)

pop_info <- fread(xargs$meta,
                  select = c("population", "collection_year", "latitude"))

small_vcf <- fread(xargs$vcf,
                   skip = "#CHROM",
                   nrows = 1,
                   drop = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                            "INFO","FORMAT"))

## parse -----------------------------------------
pi <- rbindlist(list(pi_auto, pi_X))
window_sizes <- rbindlist(list(window_sizes, window_sizes_X))

new_colnames <- c("CHROM", "window", colnames(small_vcf))

lapply(list(pi, window_sizes),#, Dtajima, ),
       function(x) setnames(x, colnames(x), new_colnames))

melted_stats <- lapply(list(pi, window_sizes),
                       function(x) melt.data.table(x, id.vars = c(1,2),
                                                   measure.vars = 3:length(x),
                                                   variable.name = "population",
                                                   value.name = "stat"))
setnames(melted_stats[[2]], "stat", "sizes")

merged_data <- 
merge.data.table(melted_stats[[1]], melted_stats[[2]],
                 by = c("CHROM", "population", "window"))

## Compute chromosomal stat --------------------------------
merged_data[, weighted_stat := stat*sizes]

merged_data[, sizes := ifelse(is.na(stat), NA, sizes)]

global_stat <- 
merged_data[, .(whole_size = sum(sizes, na.rm = TRUE),
                sum_stats = sum(weighted_stat, na.rm = TRUE)#,
                ),
            by = c("CHROM", "population")][, global_stat := sum_stats/whole_size]

global_stat <- 
global_stat[pop_info, on = "population", nomatch = 0][
  ,plot_collection_year :=
    fcase(collection_year == "2009" | collection_year == "2010", "2009/2010", 
          collection_year == "2022" | collection_year == "2023", "2022/2023",
          rep(TRUE, .N), as.character(collection_year))]


global_stat[, state := as.character(fcase(population == "CMD97B", "MD",
                               population %like% "HFL" | population %like% "dlFL", "sFL",
                               rep(TRUE, .N), as.character(str_sub_all(population, start = -4, -3))))]

global_stat[,state := fct_reorder(as.factor(state), latitude)]


global_stat[, population := fct_reorder(as.factor(population), collection_year)]

global_stats_supl <- 
dcast.data.table(global_stat, formula = CHROM ~ population,
                 value.var = "global_stat")

fwrite(global_stats_supl, xargs$output, sep = "\t")

BAX_GLOBA_stat <- 
global_stat |>
  ggplot() +
  geom_bar(aes(x = population, y = global_stat, fill = plot_collection_year),
           stat = "identity") +
  facet_grid(rows = vars(CHROM), scales = "free", space='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

jpeg(filename = xargs$statBar,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)

BAX_GLOBA_stat

dev.off()



global_stat_per_year <- 
global_stat[, .(mean(global_stat)),
            by = c("CHROM", "plot_collection_year")]
global_stat_per_year_noCMD97B <- 
global_stat[population != "CMD97B" & population %like% "97"][
  , .(mean(global_stat)),
  by = c("CHROM", "plot_collection_year")][
    , plot_collection_year := "1997_noCMD97B"
  ]

global_stat_per_year <- 
rbindlist(list(global_stat_per_year, global_stat_per_year_noCMD97B))

global_stat_per_year <-
dcast(global_stat_per_year,
      CHROM ~ plot_collection_year,
      value.var = "V1")

fwrite(global_stat_per_year, xargs$yearStat, sep = "\t")



