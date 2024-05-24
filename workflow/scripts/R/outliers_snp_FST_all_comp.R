#!/usr/bin/env Rscript
# This script takes grenedalf SNP fst output, 

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "generates plots and file with maxSNPfst cutoffs from FST windows")
parser$add_argument('--sfst', '-sfst', help= 'output from grenedalf fst')
parser$add_argument('--tfst', '-tfst', help= 'output: tidy table with fst and cutoffs')
parser$add_argument('--FSTfig', '-fig', help= 'output: jpeg manhattan plot of fst')
xargs<- parser$parse_args()

SNP_fst <- fread(xargs$sfst)
file_name <- xargs$sfst

#fix colnames 
fst_cols_index <- c(5:length(SNP_fst))
fst_cols <- colnames(SNP_fst)[fst_cols_index]
new_col_names <- fst_cols |> str_replace(":", ".")
setnames(SNP_fst, fst_cols, new_col_names)

#prep

select_cols <- function(data, cols){
  data[, cols, with = FALSE]
}

windows <- seq(0, (nrow(SNP_fst)+250), by = 250) #daria para window ser um argumento
fst_cutoffs <- c(0.01, 0.005, 0.001)
SNP_fst_auto <- SNP_fst[chrom != "X"]
SNP_fst_X <- SNP_fst[chrom == "X"]

get_maxSNPfst_cutoffs <- 
function(data, cutoffs, windows) {
  #prep
  test_names <- colnames(data)[5:length(data)]
  cols_per_test <- map(test_names, c, colnames(data)[1:4])
  
  #get maxSNP fst
  data_per_test <- map(cols_per_test, select_cols, data = data)
  data_per_test <- lapply(data_per_test, function(x) x[, I := .I])
  data_per_test <- lapply(data_per_test, function(x) x[, window := cut(I, windows)])
  data_per_test <- lapply(data_per_test, function(x) x[, .SD[which.max(.SD[[1]])],
                                                       by = window])
  
  for (i in seq_along(data_per_test)) {
    setorderv(data_per_test[[i]], test_names[[i]], -1)
  }

  pos_cut <- vector("integer", length(cutoffs))
  for (j in seq_along(cutoffs)) {
    pos_cut[[j]] <- ceiling(nrow(data_per_test[[1]]) * cutoffs[[j]])
    
  }
  
  value_cut <- vector("integer", length(data_per_test))
  all_values_cut <- vector("list", length(pos_cut))
  for (j in seq_along(pos_cut)) {
    value_cut <- sapply(data_per_test, function(x) x[pos_cut[[j]]][,2])
    all_values_cut[[j]] <- value_cut
  }
  
  
  for (i in seq_along(data_per_test)) {
    data_per_test[[i]][, cutoff := fcase(.SD[[2]] >= all_values_cut[[3]][[i]], fst_cutoffs[[3]],
                                         .SD[[2]] >= all_values_cut[[2]][[i]], fst_cutoffs[[2]],
                                         .SD[[2]] >= all_values_cut[[1]][[i]], fst_cutoffs[[1]],
                                         default = 100)]
    
  }
  
  data_per_test <- 
    lapply(data_per_test, function(x) melt.data.table(x, measure.vars = 2,
                                                      variable.name = "test",
                                                      value.name = "FST"))
  rbindlist(data_per_test)
  
  
}

maxSNPfst_X <- get_maxSNPfst_cutoffs(data = SNP_fst_X, cutoffs = fst_cutoffs, windows = windows)
maxSNPfst_auto <- get_maxSNPfst_cutoffs(data = SNP_fst_auto, cutoffs = fst_cutoffs, windows = windows)

maxSNPfst <- merge.data.table(maxSNPfst_auto, maxSNPfst_X, all = TRUE)


fwrite(maxSNPfst, xargs$tfst, sep = "\t")


maxSNPfst[, plot_cutoff := fcase(cutoff == fst_cutoffs[[1]], "0.1%",
                                  cutoff == fst_cutoffs[[2]], "0.5%",
                                  cutoff == fst_cutoffs[[3]], "1%",
                                  cutoff == 100, "-")]

plot_table <- maxSNPfst[, start := (start/1e+06)]

local <- 
  str_split_i(string = file_name, pattern = "/", i = 7) |>
  str_split_i(pattern = "_", i = 1)


SNP_FST <-
plot_table |>
  ggplot(aes()) +
  geom_point(aes(x = start, y = FST, color = plot_cutoff),
             size = 0.01) +
  facet_grid(rows = vars(test),
             cols = vars(chrom),
             scales = "free_x") +
  labs(title = paste0(local, " maxSNP FST"),
       caption =  ("250 snp windows"),
       x = "Genomic position in Mb") +
  scale_color_manual(values = c("-" = "black",
                                "0.1%" = "red",
                                "0.5%" = "green",
                                "1%" = "blue")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.1, hjust=0),
        strip.background =element_rect(fill="white"),
        strip.text =element_text(color="black"))

jpeg(filename = xargs$FSTfig,
     width = 25,
     height = 25,
     units = "cm",
     res = 1200)

SNP_FST

dev.off()


