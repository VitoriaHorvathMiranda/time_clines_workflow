#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)


#parse arguments 
parser <- ArgumentParser(description= "makes inversion frequency plots")
parser$add_argument('--meta', '-meta', help = "mateadata table")
parser$add_argument('--InvFreq', '-invf', help= 'Inversion frequencies (output from InvFreq.py)')
parser$add_argument('--output', '-out', 
                    help= 'plot: inversion frequency x latitude, jpeg')
parser$add_argument('--outputGLM', '-glm', 
                    help= 'summaries of all inversion freq ~ latitude + year glms')
xargs<- parser$parse_args()

#meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
meta <- fread(xargs$meta)
pop_info <- meta[, .(population, collection_year, collection_month,
                     latitude, n_females, flies_per_lin)]

pop_info[, chrom_n := ifelse(flies_per_lin == 2, n_females*1.5, n_females*2)]

#inv_freq <- fread("/dados/time_clines/analysis/inversions/inversion_freq.tsv")
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
  ggplot(aes(x = latitude, y = freq, color = test_year, weight=chrom_n)) +
  geom_point() +
  geom_text(aes(label = population, y = freq*1.12), size = 2) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
                                 #weights = chrom_n),
              level=0.9) +
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


##---------------- GLM summaries

inversions <- 
  unique(inv_freq$Inv)


models <- 
  lapply(inversions, function(x) glm(formula = 
                                       freq ~ latitude + test_year,
                                     data = inv_freq[Inv == x],
                                     family = binomial(),
                                     weights = chrom_n))

sink(xargs$outputGLM)
inversions
lapply(models, summary)
sink()


### ONLY 3R - for seged ---------------------------------------------------------
# INV_3RPayne <- 
#   inv_freq[Inv == "In(3R)Payne"] |>
#   ggplot(aes(x = latitude, y = freq, color = test_year)) +
#   geom_point(size = 5) +
#   geom_text_repel(aes(label = population), size = 3) +
#   geom_smooth(method = "lm", level=0.9) +
#   facet_wrap(vars(Inv), scales = "free_y") + 
#   scale_color_brewer(palette = "Dark2") +
#   theme_minimal() +
#   labs(y = "Frequency", x = "Latitude", color = "") +
#   theme(strip.text = element_text(face = "italic", size = 16),
#         axis.title = element_text(size = 16))
# 
# 
# 
# jpeg(filename = "in3RPayne_freq.jpeg",
#      width = 25,
#      height = 10,
#      units = "cm",
#      res = 1200)
# 
# INV_3RPayne

#dev.off()
