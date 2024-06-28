#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(RColorBrewer)


#parse arguments 
parser <- ArgumentParser(description= "makes inversion frequency plots")
parser$add_argument('--meta', '-meta', help = "mateadata table")
parser$add_argument('--InvFreq', '-invf', help= 'Inversion frequencies (output from InvFreq.py)')
parser$add_argument('--output', '-out', 
                    help= 'plot: inversion frequency x latitude, jpeg')
xargs<- parser$parse_args()


meta <- fread(xargs$meta)
pop_info <- meta[, .(population, collection_year, collection_month, latitude)]


inv_freq <- fread(xargs$InvFreq)

inv_freq <- 
melt.data.table(inv_freq, id.vars = "Inv",
                measure.vars = 2:length(inv_freq),
                variable.name = "population",
                value.name = "freq")



inv_freq <- merge.data.table(inv_freq, pop_info, by = "population")


inv_freq[, test_year :=
           fcase(collection_year == 2009 | collection_year == 2010, "2009/2010",
                 collection_year == 2022 | collection_year == 2023, "2022/2023",
                 rep(TRUE, .N), as.character(collection_year))]

INV_FREQ_PLOT <- 
inv_freq |>
  ggplot(aes(x = latitude, y = freq, color = test_year)) +
  geom_point() +
  geom_text(aes(label = population, y = freq*1.12), size = 2) +
  geom_smooth(method = "lm", level=0.9) +
  facet_wrap(vars(Inv), scales = "free_y") + 
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(y = "Frequency", x = "Latitude", color = "") +
  theme(strip.text = element_text(face = "italic"))



jpeg(filename = xargs$output,
     width = 26,
     height = 19,
     units = "cm",
     res = 1200)

INV_FREQ_PLOT

dev.off()

