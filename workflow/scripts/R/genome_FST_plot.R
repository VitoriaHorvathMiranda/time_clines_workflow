#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(viridis)
library(patchwork)

#parse arguments
parser <- ArgumentParser(description= "Plots global FST")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--FSTX', '-x', help = 'FST matrix with X chrom')
parser$add_argument('--FSTauto', '-auto', help= 'FST matrix with all autosomes')
parser$add_argument('--outputAll', '-oall', help= 'A jpeg with FST between all pops')
parser$add_argument('--outputTime', '-otime', help= 'A jpeg with FST between pops from 1997 and 2009/2010 separated')
xargs<- parser$parse_args()


autosomes <- fread(xargs$FSTauto)
Xchrom <- fread(xargs$FSTX)
meta <- fread(xargs$metadata)
pop_info <- meta |>
  select(population, latitude, collection_year)

tidy_autosomes <- autosomes |>
  pivot_longer(cols = c(2:length(autosomes)), names_to = "sample2", values_to = "FST_auto")

tidy_X <- Xchrom |>
  pivot_longer(cols = c(2:length(Xchrom)), names_to = "sample2", values_to = "FST_X")

joint_table <- tidy_autosomes |>
  left_join(tidy_X) |> 
  left_join(pop_info, by = join_by(sample == population)) |>
  left_join(pop_info, by = join_by(sample2 == population)) |>
  mutate(sample = fct_reorder(factor(sample), latitude.x),
         sample2 = fct_reorder(factor(sample2), latitude.y))

table_pairs <- 
  joint_table |>
  arrange(sample) |>
  rowwise() |>
  mutate(pair = sort(c(sample, sample2)) %>% paste(collapse = ",")) %>% 
  group_by(pair) %>%
  distinct(pair, .keep_all = T)

pairs_auto <- table_pairs$pair

plot_table <- 
  joint_table |>
  mutate(pair = paste0(sample,",",sample2),
         upper = if_else(pair %in% pairs_auto, 1, 0))  |>
  mutate(label_auto = round(FST_auto, digits = 2),
         label_X = round(FST_X, digits = 2))

ALL_GLOBAL_FST <- 
plot_table |>
  ggplot(aes(x= sample, y = sample2)) +
  geom_tile(aes(fill = ifelse(upper==1, FST_auto, FST_X))) +
  geom_text(aes(label = ifelse(upper==1, label_auto, label_X)), size = 2) +
  scale_fill_viridis(option = "plasma", name  = "FST") +
  labs(y = "", x = "",
       caption = "Upper triangle: autosomes \n Botton triangule: X chromosome") +
  theme(axis.text = element_text(size = 8))
  
Global_97 <- 
plot_table |>
  filter(collection_year.x == "1997",
         collection_year.y == "1997") |>
  ggplot(aes(x= sample, y = sample2)) +
  geom_tile(aes(fill = ifelse(upper==1, FST_auto, FST_X))) +
  geom_text(aes(label = ifelse(upper==1, label_auto, label_X)), size = 2) +
  scale_fill_viridis(option = "plasma", name  = "FST") +
  labs(y = "", x = "",
       #caption = "Upper triangle: autosomes \n Botton triangule: X chromosome",
       title = "1997") +
  theme_minimal()

Global_0910 <- 
plot_table |>
  filter(collection_year.x %in% c("2009", "2010"),
         collection_year.y %in% c("2009", "2010")) |>
  ggplot(aes(x= sample, y = sample2)) +
  geom_tile(aes(fill = ifelse(upper==1, FST_auto, FST_X))) +
  geom_text(aes(label = ifelse(upper==1, label_auto, label_X)), size = 2) +
  scale_fill_viridis(option = "plasma", name  = "FST") +
  labs(y = "", x = "",
       caption = "Upper triangle: autosomes \n Botton triangule: X chromosome",
       title = "2009/2010") +
  theme_minimal()


jpeg(filename = xargs$outputAll,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
ALL_GLOBAL_FST

dev.off()

jpeg(filename = xargs$outputTime,
     width = 23,
     height = 19,
     units = "cm",
     res = 1200)
Global_97 / Global_0910

dev.off()

