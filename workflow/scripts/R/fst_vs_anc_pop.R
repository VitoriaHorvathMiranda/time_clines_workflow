#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(RColorBrewer)

#parse arguments
parser <- ArgumentParser(description= "boxplots FST per year and lm fst ~year")
parser$add_argument('--FST', '-fst', help= 'FST pair list')
parser$add_argument('--outputLM', '-lm', help= 'output with lm summarys')
parser$add_argument('--plot', '-plot', help = 'scatter plot of FST per latitude per collection year per ref pop, jpeg')
parser$add_argument('--meta', '-meta')
xargs<- parser$parse_args()

# meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
#               select = c("population", "latitude", "longitude"))
# fst_2 <- fread("/dados/time_clines/analysis/ancestry/FST_ancestry/Genome_FST_autosome_fst-list.csv")

meta <- fread(xargs$meta,
              select = c("population", "latitude", "longitude"))
fst_2 <- fread(xargs$FST)

setnames(fst_2, colnames(fst_2), c("first", "second", "fst"))

fst <- 
fst_2[(first %in% c("EU", "AFR")) | (second %in% c("EU", "AFR")),
      .(first, second, fst)]

fst[, year := fcase(first %like% "10" | first %like% "09", "2009/2010",
                    first %like% "97", "1997",
                    first %like% "17", "2017",
                    first %like% "22" | first %like% "23", "2022/2023")]

fst <- 
  merge.data.table(fst, meta, by.x = "first", by.y = "population")

PLOT_LAT <- 
fst |>
  ggplot(aes(x = latitude, y = fst, shape = second, color = year)) +
  geom_point(size = 5) +
  geom_smooth(data = fst[year %in% c("1997", "2009/2010")], method = lm) +
  labs(shape = "", y = expression("F"[ST]), x = "Latitude") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

lm_eu <- 
lm(fst ~ latitude + year,
   data = fst[second == "EU" & year %in% c("2009/2010", "1997")])

lm_eu_int <- 
  lm(fst ~ latitude * year,
     data = fst[second == "EU" & year %in% c("2009/2010", "1997")])

anova(lm_eu, lm_eu_int)

lm_afr <- 
  lm(fst ~ latitude + year,
     data = fst[second == "AFR" & year %in% c("2009/2010", "1997")])

lm_afr_int <- 
  lm(fst ~ latitude * year,
     data = fst[second == "AFR" & year %in% c("2009/2010", "1997")])

sink(xargs$outputLM)
print(summary(lm_eu))
print(summary(lm_eu_int))
anova(lm_eu, lm_eu_int)


print(summary(lm_afr))
print(summary(lm_afr_int))
anova(lm_afr, lm_afr_int)
sink()

pdf(file = xargs$plot, #"fst_vs_eu_afr_per_lat_autossome.pdf", #
    width = 16/1.5, 
    height = 12/1.5)

PLOT_LAT

dev.off()

# jpeg(filename = "fst_vs_eu_afr_per_lat_autossome.jpeg", 
#      width = 27,
#      height = 15,
#      units = "cm",
#      res = 1200)
# 
# PLOT_LAT
# 
# dev.off()








