#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(nasapower)
library(patchwork)

#parse arguments 
parser <- ArgumentParser(description= "Makes pca plots, all outputs are .jpeg")
parser$add_argument('--meta', '-m', help = "matadata table")
parser$add_argument('--freqs', '-f', help= 'output from freq_extraction_pop_ind.R')
parser$add_argument('--PC1x2all', '-p12All', help= 'PCA (1 x 2) of all pops in freq table')
parser$add_argument('--PC3x4all', '-p34All', help= 'PCA (3 X 4) of all pops in freq table')
parser$add_argument('--PCAtimes', '-pTimes', help= 'PCA of 1997 and 2009/2010 separated')
parser$add_argument('--PCAchrom', '-PCAchrom', help= 'PCA per chrom')
parser$add_argument('--lmALL', '-lm', help= 'Table with lm values (p, R, R-adjs) for first 4 PCs')
parser$add_argument('--Rsq', '-Rsq', help= 'Table with R-square for first 4 PCs')
xargs<- parser$parse_args()

#get envi variables ---------------------------------------------------------------
# seq_metadata <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
seq_metadata <- fread(xargs$meta)

meta <- seq_metadata |>
  select(population, collection_year, latitude, longitude, collection_month)
meta <- 
  meta[!(population %in% c("ESC97", "SNC10", "HFL97_new"))]

nasa_data <- vector("list", length = nrow(meta))
for (i in seq_along(1:nrow(meta))) {
  test_lonlat <- c(meta$longitude[i],meta$latitude[i])
  test_year <- c(meta$collection_year[i]-1,meta$collection_year[i]) 
  test_date <- c(paste((meta$collection_year[i]-1), meta$collection_month[i], "30",
                       sep = "-"),
                 paste((meta$collection_year[i]), meta$collection_month[i], "30",
                       sep = "-")) |>
    as_date()
  nasa_data[[i]] <- get_power(
    community = "ag",
    lonlat = test_lonlat,
    temporal_api = "daily",
    dates = test_date,
    pars = c("PRECTOTCORR","T2M")
  ) |> 
    mutate(population = meta$population[i])
}

nasa_data <- lapply(nasa_data, as.data.table)
nasa_data_tidy <- rbindlist(nasa_data)

envi_info <- 
  meta[nasa_data_tidy, on = c("population",
                                "latitude" = "LAT",
                                "longitude" = "LON")]

envi_info[, diff_to_c_month := ifelse(collection_year == YEAR,
                                      collection_month-MM,
                                      (collection_month-MM)+12)]


envi_6months <- 
  envi_info[diff_to_c_month <= 6 &
              diff_to_c_month >= 0][,
                                    .(tmin = min(T2M),
                                      tmax = max(T2M),
                                      tmean = mean(T2M),
                                      preci_mean = mean(PRECTOTCORR)),
                                    by = c("population","MM",
                                           "collection_month", "diff_to_c_month")]

envi_test_6months <- envi_6months[, .(tmin_6mday = min(tmin), tmax_6mday = max(tmax),
                                     tmin_6month = min(tmean), tmax_6month = max(tmean),
                                     tmean_6months = mean(tmean),
                                     preci_max_6month = max(preci_mean),
                                     preci_mean_6months = mean(preci_mean)),
                                 by = c("population")]


envi_3months <- 
  envi_info[diff_to_c_month <= 3 &
              diff_to_c_month >= 0][,
                                    .(tmin = min(T2M),
                                      tmax = max(T2M),
                                      tmean = mean(T2M),
                                      preci_mean = mean(PRECTOTCORR)),
                                    by = c("population","MM",
                                           "collection_month", "diff_to_c_month")]

envi_test_3months <- envi_3months[, .(tmin_3mday = min(tmin), tmax_3mday = max(tmax),
                                     tmin_3month = min(tmean), tmax_3month = max(tmean),
                                     tmean_3months = mean(tmean),
                                     preci_max_3month = max(preci_mean),
                                     preci_mean_3months = mean(preci_mean)),
                                 by = c("population")]



# gets freqs and all samples PCA --------------------------------------------------

freqs <- fread(xargs$freqs)

#freqs <- fread("/dados/time_clines/data/seqs/calls/freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")

chroms <- c("2L", "2R", "3L", "3R", "X")
freqs <- freqs[depth >= 20]


all_samples <- freqs[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
all_samples2 <- all_samples[, !c("CHROM", "POS")]
all_samples2[, freq := as.double(freq)]


# ALL pops PCA ---------------------------------------------------------------

pre_pca_table <- dcast.data.table(all_samples2, position2 ~ population,
                                  value.var = "freq")

fun_pca <- 
function(table){
  no_NA_table <- na.omit(table)
  pca_table <- no_NA_table |> 
    column_to_rownames(var="position2")
  
  pca_table <- t(pca_table)
  prcomp(pca_table)
}

pca_all <- fun_pca(pre_pca_table)

pca_var <- pca_all$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

pcs <- rownames_to_column(as.data.frame(pca_all$x), var = "population")
pcs <- as.data.table(pcs)


# lm dos pcs contra variáveis --------------------------------------------------

pc_info <- as.data.table(left_join(pcs, meta))
pc_info <- as.data.table(left_join(pc_info, envi_test_3months))
pc_info <- as.data.table(left_join(pc_info, envi_test_6months))

PC_names <- paste0("PC", c(1:length(pca_var)))

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

variables <- colnames(pc_info)[(length(pca_var)+2):length(pc_info)]

model_p_values <- vector("list", length(variables))
r_squared <- vector("list", length(variables))
r_squared_adj <- vector("list", length(variables))

for (i in seq_along(variables)) {
  
  test_models <-
    sapply(PC_names, function(x) paste(x, "~ ", variables[[i]])) |>
    lapply(function(x) lm(data = pc_info, formula = x))
  
  model_p_values[[i]] <- lapply(test_models, lmp)
  
  r_squared[[i]] <- lapply(test_models, summary) |>
    lapply(function(x) x$r.squared)
  
  r_squared_adj[[i]] <- lapply(test_models, summary) |>
    lapply(function(x) x$adj.r.squared)
  
}

models_info <- lapply(list(model_p_values, r_squared, r_squared_adj),
                      rbindlist) |>
  lapply(function(x) x[, test_variable := variables]) |>
  lapply(function(x) setcolorder(x, "test_variable")) 

test <- 
lapply(models_info, function(x) 
  melt.data.table(data = x,
                  id.vars = "test_variable",
                  measure.vars = 2:5,
                  variable.name = "resposta",
                  value.name = "value"))

setnames(test[[1]], "value", "pvalue")
setnames(test[[2]], "value", "r_squared")
setnames(test[[3]], "value", "r_squared_adj")

all_tests <- reduce(test, left_join)

all_tests[, r_sq := fcase(pvalue <= 0.01, paste0(r_squared, "***"),
                          pvalue <= 0.05, paste0(r_squared, "**"),
                          pvalue <= 0.1, paste0(r_squared, "*"),
                          pvalue <= 0.12, paste0(r_squared, "."),
                          rep(TRUE, .N), as.character(r_squared))]

fwrite(all_tests, file = xargs$lmALL, sep = "\t")

p_table <- 
dcast.data.table(data = all_tests,
                 formula = test_variable ~ resposta, 
                 value.var = "r_sq")

fwrite(p_table, xargs$Rsq, sep = "\t")

#plot PC1x2 e PC3x4 todas as variáveis -----------------------------------------

# graph_table <- tibble(population = rownames(pca_all$x), 
#                       pc1 = pca_all$x[,1], #gets the first four PCs
#                       pc2 = pca_all$x[,2],
#                       pc3 = pca_all$x[,3],
#                       pc4 = pca_all$x[,4]) %>% 
#   inner_join(pop_info, by = "population") %>%  #joins tables
#   mutate(collection_year = case_when( #treats collection_year as character
#     collection_year == 1997 ~ "1997", 
#     collection_year == 2009 | collection_year == 2010 ~ "2009/2010", #and join 2009 and 2010
#     collection_year == 2022 | collection_year == 2023 ~ "2022/2023",
#     collection_year == 2017 ~ "2017"
#   )) |>
#   mutate(population = if_else(population == "HFL97_new", "ESC97", population))

PC1x2_ALL <- 
pc_info |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.7) +
  scale_color_viridis_c(option = "B") +
  geom_text(aes(label = population), size = 1.5, nudge_y = 3) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_minimal()

PC3x4_ALL <- 
  pc_info |>
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = collection_year), size = 4, alpha = 0.5) +
  scale_color_viridis_b(option = "C") +
  geom_text(aes(label = population), size = 1.5, nudge_y = 3) +
  xlab(paste("PC3 - ", pca_var_per[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_var_per[4], "%", sep = "")) +
  labs(color = "Collection Year") +
  theme_minimal()


# 1997 pops PCA ----------------------------------------------------------------

pops_97 <- all_samples2[population %like% "97"]

pre_pca_table97 <- dcast.data.table(pops_97, position2 ~ population,
                                  value.var = "freq")

pca97 <- fun_pca(pre_pca_table97)

pca_var97 <- pca97$sdev^2
pca_var_per97 <- round(pca_var97/sum(pca_var97)*100, 1)


graph_table97 <- tibble(population = rownames(pca97$x), 
                      pc1 = pca97$x[,1], #gets the first four PCs
                      pc2 = pca97$x[,2],
                      pc3 = pca97$x[,3],
                      pc4 = pca97$x[,4]) %>% 
  inner_join(meta, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 | collection_year == 2010 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2022 | collection_year == 2023 ~ "2022/2023",
    collection_year == 2017 ~ "2017"
  )) |>
  mutate(population = if_else(population == "HFL97_new", "ESC97", population))

PCA_97 <- 
graph_table97 |>
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per97[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per97[2], "%", sep = "")) +
  labs(title = "1997") +
  theme_light() +
  theme(plot.title = element_text(size = 16))

# 2009/2010 pops PCA ----------------------------------------------------------------

pops_0910 <- all_samples2[population %like% "09" | population %like% "10"]


pre_pca_table0910 <- dcast.data.table(pops_0910, position2 ~ population,
                                    value.var = "freq")

pca0910 <- fun_pca(pre_pca_table0910)


pca_var0910 <- pca0910$sdev^2
pca_var_per0910 <- round(pca_var0910/sum(pca_var0910)*100, 1)


graph_table0910 <- tibble(population = rownames(pca0910$x), 
                        pc1 = pca0910$x[,1], #gets the first four PCs
                        pc2 = pca0910$x[,2],
                        pc3 = pca0910$x[,3],
                        pc4 = pca0910$x[,4]) %>% 
  inner_join(meta, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 | collection_year == 2010 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2022 | collection_year == 2023 ~ "2022/2023",
    collection_year == 2017 ~ "2017"
  )) |>
  mutate(population = if_else(population == "HFL97_new", "ESC97", population))

PCA_0910 <- 
graph_table0910 |>
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per0910[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per0910[2], "%", sep = "")) +
  labs(title = "2009/2010") +
  theme_light() +
  theme(plot.title = element_text(size = 16),
        legend.position = "none")

### ALL per chrom --------------------------------------------------------------


#separate each chrom
output <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  output[[i]] <- all_samples2[position2 %like% chroms[[i]],]
}

output <-
  lapply(output, dcast.data.table,position2 ~ population, value.var = "freq")

pcas_per_chrom <- lapply(output, fun_pca)

#gets pc1 and 2 for each chrom
pcs1e2_per_chrom <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  temp <- pcas_per_chrom[[i]]$x[,1:2]
  pcs1e2_per_chrom[[i]] <- as.data.frame(temp)
}

#change cols names
for (i in seq_along(pcs1e2_per_chrom)) {
  pcs1e2_per_chrom[[i]] <- setnames(pcs1e2_per_chrom[[i]],
                                    colnames(pcs1e2_per_chrom[[i]]),
                                    paste(colnames(pcs1e2_per_chrom[[i]]),
                                          chroms[[i]], sep = "_"))
}


all_chrom_pcs1e2 <- reduce(pcs1e2_per_chrom, cbind) |>
  rownames_to_column("population")

graph_table_per_chrom <- all_chrom_pcs1e2 |>
  inner_join(meta, by = "population") |>
  pivot_longer(cols = 2:11, names_to = "PCs",
               values_to = "var") %>%
  separate(col = PCs, into = c("PC", "chrom"), sep = "_") |>
  pivot_wider(names_from = PC, values_from = var)


PCA_PER_CHROM <- 
graph_table_per_chrom |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = latitude), size = 2, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  facet_wrap(~ chrom, nrow = 3) +
  labs(title = "") +
  theme_light() +
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "white"),
        title = element_text(size = 14))

### 2009/2010 per chrom --------------------------------------------------------------
#separate each chrom
# output0910 <- vector("list", length = length(chroms))
# for (i in seq_along(chroms)) {
#   output0910[[i]] <- pops_0910[position2 %like% chroms[[i]],]
# }
# 
# output0910 <-
#   lapply(output0910, dcast.data.table,position2 ~ population, value.var = "freq")
# 
# pcas_per_chrom_0910 <- lapply(output0910, fun_pca)
# 
# #gets pc1 and 2 for each chrom
# pcs1e2_per_chrom0910 <- vector("list", length = length(chroms))
# for (i in seq_along(chroms)) {
#   temp <- pcas_per_chrom_0910[[i]]$x[,1:2]
#   pcs1e2_per_chrom0910[[i]] <- as.data.frame(temp)
# }
# 
# #change cols names
# for (i in seq_along(pcs1e2_per_chrom0910)) {
#   pcs1e2_per_chrom0910[[i]] <- setnames(pcs1e2_per_chrom0910[[i]],
#                                       colnames(pcs1e2_per_chrom0910[[i]]),
#                                       paste(colnames(pcs1e2_per_chrom0910[[i]]),
#                                             chroms[[i]], sep = "_"))
# }
# 
# 
# all_chrom_pcs1e2_0910 <- reduce(pcs1e2_per_chrom0910, cbind) |>
#   rownames_to_column("population")
# 
# graph_table_per_chrom0910 <- all_chrom_pcs1e2_0910 |>
#   inner_join(pop_info, by = "population") |>
#   pivot_longer(cols = 2:11, names_to = "PCs",
#                values_to = "var") %>%
#   separate(col = PCs, into = c("PC", "chrom"), sep = "_") |>
#   pivot_wider(names_from = PC, values_from = var)
# 
# PCA_PER_CHROM_0910 <- 
# graph_table_per_chrom0910 |>
#   ggplot(aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = latitude), size = 2, alpha = 0.5) +
#   scale_color_viridis_c() +
#   geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
#   facet_wrap(~ chrom, nrow = 3) +
#   labs(title = "2009/2010") +
#   theme_light() +
#   theme(strip.text = element_text(size = 12, color = "black"),
#         strip.background = element_rect(fill = "white"),
#         title = element_text(size = 14))

####### saves images -----------------------------------------------

jpeg(filename = xargs$PC1x2all,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PC1x2_ALL
dev.off()


jpeg(filename = xargs$PC3x4all,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PC3x4_ALL
dev.off()

jpeg(filename = xargs$PCAtimes,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_97 / PCA_0910

dev.off()


jpeg(filename = xargs$PCAchrom,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_PER_CHROM

dev.off()






