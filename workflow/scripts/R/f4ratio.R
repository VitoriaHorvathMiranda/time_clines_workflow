#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(poolfstat)
library(RColorBrewer)

#parse arguments 
parser <- ArgumentParser(description= "transformns ancestral pops vcf into sync")
parser$add_argument('--sync', '-sync', help = "sync file with all pops")
parser$add_argument('--sizes', '-sizes', help= 'a file with pop names and poolsizes in the same order as .sync')
parser$add_argument('--meta', '-meta', help= 'metadata table')
parser$add_argument('--output', '-o', help= 'path to output files')
parser$add_argument('--plot', '-plot', help= 'plot f4ratio per latitude with error bars')
xargs<- parser$parse_args()


# meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
# pop_info <- meta[, .(population, collection_year, collection_month, latitude)]
# 
# pool_sizes <- fread("~/time_clines_workflow/resources/pool_sizes_autosome.csv")
# 
# pool_data <- popsync2pooldata(sync.file = "/dados/time_clines/analysis/ancestry/f4ratio/all.sync",
#                               poolsizes = c(23, 23, 64, 32, pool_sizes$V2),
#                               poolnames = c("EU", "westAFR", "ZI", "EG", pool_sizes$V1),
#                               min.cov.per.pool = 10,
#                               min.maf = 0.05,
#                               nthreads = 5)
meta <- fread(xargs$meta)
pop_info <- meta[, .(population, collection_year, collection_month, latitude)]

pool_sizes <- fread(xargs$sizes)

pool_data <- popsync2pooldata(sync.file = xargs$sync,
                              poolsizes = c(23, 23, 64, 32, pool_sizes$V2),
                              poolnames = c("EU", "westAFR", "ZI", "EG", pool_sizes$V1),
                              min.cov.per.pool = 10,
                              min.maf = 0.05,
                              nthreads = 5)
res_fstats <- compute.fstats(pool_data, 
                             nsnp.per.bjack.block = 1e4,
                             return.F4.blockjackknife.samples = TRUE)

f4ratio <- 
  lapply(pool_sizes$V1, function(x) {
    compute.f4ratio(res_fstats,
                    num.quadruplet = paste0("EG,ZI;westAFR,",x),
                    den.quadruplet = "EG,ZI;westAFR,EU")
  })

f4table <- 
lapply(f4ratio, function(x){
  as.data.frame(x) |> 
    rownames_to_column() |>
    pivot_wider(names_from = rowname,
                values_from = x) |>
    data.table() 
}) |> rbindlist()


f4table[, nrow := .I]
pool_sizes[, nrow := .I]

f4table <- 
  f4table[pool_sizes, on = "nrow"][pop_info, on = .(V1 = population),
                                   nomatch = 0]

fwrite(f4table, file = xargs$output, sep = "\t")

f4table[, year := fcase(collection_year %in% c(2010,2009), "2009/2010",
                        collection_year == 1997, "1997",
                        collection_year %in% c(2022, 2023), "2022/2023",
                        collection_year == 2017, "2017")]
PLOT <- 
f4table|>
  ggplot(aes(x = latitude, y = Estimate, color = year)) +
  geom_point(aes(), size = 2) +
  geom_errorbar(aes(ymin = (Estimate - `bjack s.e.`),
                    ymax = (Estimate + `bjack s.e.`)),
                width = 0.25) +
  labs(x= "Latitude", y = "Alpha", color = "Collection Year") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() 

pdf(file = xargs$plot,
     width = 10,
     height = 6)

PLOT

dev.off()






