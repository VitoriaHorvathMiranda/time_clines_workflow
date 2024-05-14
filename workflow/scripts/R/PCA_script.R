#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(patchwork)

#parse arguments 
parser <- ArgumentParser(description= "Makes pca plots, all outputs are .jpeg")
parser$add_argument('--meta', '-m', help = "matadata table")
parser$add_argument('--freqs', '-f', help= 'output from freq_extraction_pop_ind.R')
parser$add_argument('--PCAall', '-pAll', help= 'PCA of all pops in freq table')
parser$add_argument('--PCAtimes', '-pTimes', help= 'PCA of 1997 and 2009/2010 separated')
parser$add_argument('--PCAchrom97', '-PC97', help= 'PCA per chrom 1997')
parser$add_argument('--PCAchrom0910', '-PC0910', help= 'PCA per chrom 2009/2010')
xargs<- parser$parse_args()

seq_metadata <- fread(xargs$meta)
pop_info <- seq_metadata |>
  select(population, collection_year, latitude)

freqs <- fread(xargs$freqs)

chroms <- c("2L", "2R", "3L", "3R", "X")

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


graph_table <- tibble(population = rownames(pca_all$x), 
                      pc1 = pca_all$x[,1], #gets the first four PCs
                      pc2 = pca_all$x[,2],
                      pc3 = pca_all$x[,3],
                      pc4 = pca_all$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 | collection_year == 2010 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2022 | collection_year == 2023 ~ "2022/2023",
    collection_year == 2017 ~ "2017"
  )) |>
  mutate(population = if_else(population == "HFL97_new", "ESC97", population))

PCA_ALL <- 
graph_table |>
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 3) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_light()



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
  inner_join(pop_info, by = "population") %>%  #joins tables
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
  inner_join(pop_info, by = "population") %>%  #joins tables
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

### 1997 per chrom --------------------------------------------------------------


#separate each chrom
output97 <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  output97[[i]] <- pops_97[position2 %like% chroms[[i]],]
}

output97 <-
  lapply(output97, dcast.data.table,position2 ~ population, value.var = "freq")

pcas_per_chrom_97 <- lapply(output97, fun_pca)

#gets pc1 and 2 for each chrom
pcs1e2_per_chrom97 <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  temp <- pcas_per_chrom_97[[i]]$x[,1:2]
  pcs1e2_per_chrom97[[i]] <- as.data.frame(temp)
}

#change cols names
for (i in seq_along(pcs1e2_per_chrom97)) {
  pcs1e2_per_chrom97[[i]] <- setnames(pcs1e2_per_chrom97[[i]],
                                    colnames(pcs1e2_per_chrom97[[i]]),
                                    paste(colnames(pcs1e2_per_chrom97[[i]]),
                                          chroms[[i]], sep = "_"))
}


all_chrom_pcs1e2_97 <- reduce(pcs1e2_per_chrom97, cbind) |>
  rownames_to_column("population")

graph_table_per_chrom97 <- all_chrom_pcs1e2_97 |>
  inner_join(pop_info, by = "population") |>
  pivot_longer(cols = 2:11, names_to = "PCs",
               values_to = "var") %>%
  separate(col = PCs, into = c("PC", "chrom"), sep = "_") |>
  pivot_wider(names_from = PC, values_from = var)


PCA_PER_CHROM_97 <- 
graph_table_per_chrom97 |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = latitude), size = 2, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  facet_wrap(~ chrom, nrow = 3) +
  labs(title = "1997") +
  theme_light() +
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "white"),
        title = element_text(size = 14))

### 2009/2010 per chrom --------------------------------------------------------------


#separate each chrom
output0910 <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  output0910[[i]] <- pops_0910[position2 %like% chroms[[i]],]
}

output0910 <-
  lapply(output0910, dcast.data.table,position2 ~ population, value.var = "freq")

pcas_per_chrom_0910 <- lapply(output0910, fun_pca)

#gets pc1 and 2 for each chrom
pcs1e2_per_chrom0910 <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  temp <- pcas_per_chrom_0910[[i]]$x[,1:2]
  pcs1e2_per_chrom0910[[i]] <- as.data.frame(temp)
}

#change cols names
for (i in seq_along(pcs1e2_per_chrom0910)) {
  pcs1e2_per_chrom0910[[i]] <- setnames(pcs1e2_per_chrom0910[[i]],
                                      colnames(pcs1e2_per_chrom0910[[i]]),
                                      paste(colnames(pcs1e2_per_chrom0910[[i]]),
                                            chroms[[i]], sep = "_"))
}


all_chrom_pcs1e2_0910 <- reduce(pcs1e2_per_chrom0910, cbind) |>
  rownames_to_column("population")

graph_table_per_chrom0910 <- all_chrom_pcs1e2_0910 |>
  inner_join(pop_info, by = "population") |>
  pivot_longer(cols = 2:11, names_to = "PCs",
               values_to = "var") %>%
  separate(col = PCs, into = c("PC", "chrom"), sep = "_") |>
  pivot_wider(names_from = PC, values_from = var)

PCA_PER_CHROM_0910 <- 
graph_table_per_chrom0910 |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = latitude), size = 2, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  facet_wrap(~ chrom, nrow = 3) +
  labs(title = "2009/2010") +
  theme_light() +
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "white"),
        title = element_text(size = 14))

####### saves images -----------------------------------------------

jpeg(filename = xargs$PCAall,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_ALL
dev.off()


jpeg(filename = xargs$PCAtimes,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_97 / PCA_0910

dev.off()


jpeg(filename = xargs$PCAchrom97,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_PER_CHROM_97

dev.off()


jpeg(filename = xargs$PCAchrom0910,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
PCA_PER_CHROM_0910

dev.off()




