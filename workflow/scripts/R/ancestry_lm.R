#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggrepel)

#parse arguments 
parser <- ArgumentParser(description= "computes global ancestry based on Bargland et al 2016 and plots it")
parser$add_argument('--freqNE', '-ne', help = "table with SNP freqs and NEs")
parser$add_argument('--rawAnc', '-anc', help= 'path to raw ancestry painel files')
parser$add_argument('--meta', '-meta', help= 'metadata table')
parser$add_argument('--allModels', '-allmod', help= 'table with all models coefficients')
parser$add_argument('--FIGallModels', '-figAll', help= 'violoin plots with all models, jpeg')
parser$add_argument('--summary', '-s', help= 'table with summary (mean, mim, max) of all models per pop')
parser$add_argument('--FIGlatModels', '-figLat', help= 'plot 1997 pops and 2009/2010 pops per lat, jpeg')
parser$add_argument('--outputLM', '-lm', help= 'lm summaries of mean ancestry ~ latitude + year')
xargs<- parser$parse_args()

# painel_files <- list.files(path = "/dados/time_clines/analysis/ancestry",
#                            pattern = "raw_ancestry_painel_chrom_*",
#                            full.names = TRUE)

painel_files <- list.files(path = xargs$rawAnc,
                           pattern = "raw_ancestry_painel_chrom_*",
                           full.names = TRUE)

painels <- lapply(painel_files, fread)

painel <- rbindlist(painels)

#freqs <- fread("/dados/time_clines/data/seqs/calls/NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
#meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
freqs <- fread(xargs$freqNE)
meta <- fread(xargs$meta)

pop_info <- meta[, .(population, collection_year, latitude)]

#get frequencies from painel
painel[, freq_EU :=
         alt_allele_count_EU/(alt_allele_count_EU + ref_allele_count_EU)]

painel[, freq_AFR :=
         alt_allele_count_AFR/(alt_allele_count_AFR + ref_allele_count_AFR)]

#total genomes where there is info about that position
painel[, n_afr := (ref_allele_count_AFR + alt_allele_count_AFR)]
painel[, n_eu := (ref_allele_count_EU + alt_allele_count_EU)]

#get difference between freqs
painel[, freq_diff := abs(freq_EU - freq_AFR)]

s_painel <- painel[, !(5:8)]

#filter based on a minimum number of genomes in each painel
#filter based on freq difference between ancestral pops (get only informative positions)
s_painel <- s_painel[n_afr >=15 & n_eu >=15]
s_painel <- s_painel[freq_diff >= 0.2]

#só fica com os SNPs que estão presentes tanto no painel quanto nas minhas pops
available_snps <- freqs[s_painel, on = .(CHROM, POS, ALT, REF),
                        nomatch=0]

#deleta as colunas criadas para filtrar o painel
available_snps <- available_snps[, !(13:15)]

#fix new sample name
available_snps[, population := ifelse(population == "HFL97_new", "ESC97", population)]

#separar as pops em listas
pops <- unique(available_snps$population)
pops_list <- vector("list", length(pops))
for (i in seq_along(pops)) {
  pops_list[[i]] <- available_snps[population == pops[[i]]]
}

#criar indíce por lista
for (i in seq_along(pops_list)) {
  pops_list[[i]][, indice := .I]
}

#para cada pop sorteia 5000 snps 100x 
#para cada um desses sorteios roda um lm
models <- rep(list(list()), length(pops))
for (i in seq_along(pops_list)) {
  for (j in seq_along(1:100)) {
    test_data <- pops_list[[i]][sample(indice, 5000)]
    lm_test <- lm(freq~0+freq_EU+freq_AFR, data = test_data, weights = NE)
    models[[i]][[j]] <- data.table(ancestry = c("EU", "AFR"),
                                   value = lm_test$coefficients,
                                   pvalue = summary(lm_test)$coefficients[,4],
                                   tvalue = summary(lm_test)$coefficients[,3],
                                   stError = summary(lm_test)$coefficients[,2])
  }
  
}


joined_models <- vector("list", length(models))
for (i in seq_along(models)) {
  joined_models[[i]] <- rbindlist(models[[i]])
  joined_models[[i]][, pop := pops[[i]]]
}


all_pops_models <- rbindlist(joined_models)
fwrite(all_pops_models, file = xargs$allModels, sep = "\t")

all_pops_models <- all_pops_models[pop_info,
                                   on = .(pop = population),
                                   nomatch=0][, plot_collection_year :=
                                                fcase(collection_year == "2009" | collection_year == "2010", "2009/2010", 
                                                      collection_year == "2022" | collection_year == "2023", "2022/2023",
                                                      rep(TRUE, .N), as.character(collection_year))]

all_pops_models[, state := as.character(fcase(pop == "CMD97B", "MD",
                                              pop %like% "HFL" | pop %like% "JFL" | pop %like% "dlFL", "sFL",
                                              rep(TRUE, .N), as.character(str_sub_all(pop, start = -4, -3))))]

all_pops_models[,state := fct_reorder(as.factor(state), latitude)]

all_pops_models[, pop := fct_reorder(as.factor(pop), collection_year)]


VIOLIN_AFR <- 
all_pops_models[ancestry == "AFR"] |>
  ggplot(aes(y= value, x= pop, fill = plot_collection_year)) +
  geom_violin() +
  #geom_histogram(binwidth = 0.002) +
  facet_grid(cols = vars(state), scales = "free", space='free') +
  labs(y = "African Ancestry", x = "", fill = "Colletion Year") +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "gray95"),
        legend.position = "none",
        axis.text.x = element_text(size = 8, angle = 35, hjust = 0.8))

VIOLIN_EU <- 
all_pops_models[ancestry == "EU"] |>
  ggplot(aes(y= value, x= pop, fill = plot_collection_year)) +
  geom_violin() +
  #geom_histogram(binwidth = 0.002) +
  facet_grid(cols = vars(state), scales = "free", space='free') +
  labs(y = "European Ancestry", x = "", fill = "Colletion Year") +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "gray95"),
        axis.text.x = element_text(size = 8, angle = 35, hjust = 0.8))


jpeg(filename = xargs$FIGallModels,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
VIOLIN_EU / VIOLIN_AFR
dev.off()

summary_models <- all_pops_models[, .(mean= mean(value), sd = sd(value),
                                      min = min(value), max= max(value),
                                      median= median(value)),
                                  by = c("ancestry", "pop")]
 
fwrite(summary_models, xargs$summary, sep = "\t")

summary_models <- summary_models[pop_info,
                                 on = .(pop = population), nomatch=0]

summary_models[, collection_year :=
                 fcase(collection_year == "2009" | collection_year == "2010", "2009/2010", 
                       collection_year == "2022" | collection_year == "2023", "2022/2023",
                       rep(TRUE, .N), as.character(collection_year))]

summary_models[, SE := sd/sqrt(100)] #n = 100 testes

#summary_models_97_0910 <- summary_models[collection_year %in% c("2009/2010", "1997")]


mean_african_ancestry <- 
  summary_models[ancestry == "AFR"] |>
  ggplot(aes(x = latitude, y = mean, color = collection_year)) +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = max, ymin = min), 
                linewidth = 0.5,
                width = 0.2,
                color = "black") +
  geom_point(aes(), size = 3) +
  geom_text_repel(aes(label = pop), size = 3) + #dificultava a visualização do texto
  #facet_wrap(vars(collection_year)) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  #coord_cartesian(ylim = c(0,1)) +
  labs(y = "Mean African Ancestry", x = "Latitude", color = "") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size = 12),
        #legend.position = "none",
        axis.title = element_text(size = 16))

mean_european_ancestry <- 
  summary_models[ancestry == "EU"] |>
  ggplot(aes(x = latitude, y = mean, color = collection_year)) +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = max, ymin = min), 
                linewidth = 0.5,
                width = 0.2,
                color = "black") +
  geom_point(aes(), size = 3) +
  geom_text_repel(aes(label = pop), size = 3) + #dificultava a visualização do texto
  #facet_wrap(vars(collection_year)) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  #coord_cartesian(ylim = c(0,1)) +
  labs(y = "Mean European Ancestry", x = "Latitude", color = "") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size = 12),
        #legend.position = "none",
        axis.title = element_text(size = 16))

# global_ancestry_poolseq_97_0910_bergland_method <- 
# summary_models_97_0910 |>
#   ggplot(aes(x = latitude, y = mean, color = ancestry)) +
#   geom_smooth(method = "lm") +
#   geom_errorbar(aes(ymax = max, ymin = min), 
#                 linewidth = 0.5,
#                 width = 0.2,
#                 color = "black") +
#   geom_point(aes(), size = 1) +
#   geom_text(data = summary_models_97_0910[pop %like% "MCT"], #fiz isso porque a posição do MCT 
#             aes(x = latitude, y = mean-0.1, label = pop), size = 2) + #dificultava a visualização do texto
#   geom_text(data = summary_models_97_0910[!(pop %like% "MCT")],
#             aes(x = latitude, y = mean-0.05, label = pop), size = 2) +
#   facet_wrap(vars(collection_year)) +
#   scale_color_brewer(palette = "Dark2") +
#   theme_light() +
#   coord_cartesian(ylim = c(0,1)) +
#   labs(y = "mean ancestry") +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(color = "black", size = 12))

jpeg(xargs$FIGlatModels,
     width = 25,
     height = 13,
     units = "cm",
     res = 1200)
#global_ancestry_poolseq_97_0910_bergland_method
mean_african_ancestry / mean_european_ancestry

dev.off()


lm_afr <- 
  lm(data = summary_models[ancestry == "AFR"],
     formula = mean ~ latitude + collection_year)

lm_eu <- 
  lm(data = summary_models[ancestry == "EU"],
     formula = mean ~ latitude + collection_year)

sink(xargs$outputLM)
summary(lm_afr)
summary(lm_eu)
sink()

