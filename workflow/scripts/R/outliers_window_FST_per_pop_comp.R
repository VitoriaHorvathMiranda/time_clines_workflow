#!/usr/bin/env Rscript
# This script takes grenedalf window fst output, 

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "generates plots and file with windown cutoffs from FST windows")
parser$add_argument('--wfst', '-wfst', help= 'output from grenedalf fst')
parser$add_argument('--chrom', '-chrom', help= 'chrom to analysis')
parser$add_argument('--tfst', '-tfst', help= 'output: tidy table with fst and cutoffs')
parser$add_argument('--wsize', '-wsize', help= 'output: table with summary window stats')
#parser$add_argument('--FSTfig', '-fig', help= 'output: jpeg manhattan plot of fst')
xargs<- parser$parse_args()

#w_fst <- fread("/dados/time_clines/analysis/fst/window/MA_100_hudson_fst.csv")
w_fst <- fread(xargs$wfst)
file_name <- xargs$wfst

keep_cols <- 
c(colnames(w_fst)[1:3],str_subset(colnames(w_fst),"fst"))
w_fst <- w_fst[, keep_cols, with = FALSE]

fst_cols_index <- c(4:length(w_fst))
fst_cols <- colnames(w_fst)[fst_cols_index]
new_col_names <- fst_cols |> str_split_i("\\.", i = 1) |> str_replace(":", ".")
setnames(w_fst, fst_cols, new_col_names)

#filter chrom
CHROM <- xargs$chrom
#CHROM <- '2L'
w_fst <- w_fst[chrom == CHROM]

#se for criar uma função, os argumentos seriam: fst_cutoffs; w_fst; new_names?
fst_cutoffs <- c(0.01, 0.005, 0.001)

get_fst_cutoffs <- 
function(fst_table, cutoffs){
  test_names <- colnames(fst_table)[4:length(fst_table)]
  print(test_names)
  w_fst_filters <- vector("integer", length(cutoffs))
  all_w_fst_cutoffs <- vector("list", length(test_names))
  for (i in seq_along(test_names)) {
    temp <- setorderv(fst_table, test_names[[i]], -1) #i é o cada par de pops
    
    for (j in seq_along(cutoffs)) {
      #pega a posição do top x% fst por cromossomo
      w_fst_filters[[j]] <- ceiling(nrow(temp)*cutoffs[[j]]) # j é cada cutoff
      #pega o ultimo fst dentro do cutoff
      w_fst_filters[[j]] <- temp[w_fst_filters[[j]],
                                 test_names[[i]],
                                 with = FALSE]
      all_w_fst_cutoffs[[i]] <- w_fst_filters
      print(all_w_fst_cutoffs[[i]][[j]][[1]])
      
    }
    
    fst_table[, paste0(test_names[[i]], "_cutoff") := 
            fcase(.SD[[i + 3]] >= all_w_fst_cutoffs[[i]][[3]][[1]], cutoffs[[3]], #0.01%
                  .SD[[i + 3]] >= all_w_fst_cutoffs[[i]][[2]][[1]], cutoffs[[2]], #0.5%
                  .SD[[i + 3]] >= all_w_fst_cutoffs[[i]][[1]][[1]], cutoffs[[1]], #1%
                  default = 100)]
    
  }
  
  pivoted_fst <- 
    melt.data.table(data = fst_table, measure.vars = test_names,
                    variable.name = "test", value.name = "FST") |> 
    melt.data.table(measure.vars = paste0(test_names, "_cutoff"),
                    variable.name = "test2", value.name = "cutoff")
   

  pivoted_fst[, test2 := str_split_i(test2, pattern = "_", i = 1)]
  pivoted_fst[test == test2]
  
}

#X_fst_complete <- get_fst_cutoffs(fst_table = w_fst_X, cutoffs = fst_cutoffs)
fst_complete <- get_fst_cutoffs(fst_table = w_fst, cutoffs = fst_cutoffs)

#plot_table <- merge.data.table(X_fst_complete, Auto_fst_complete, all = TRUE)
#fwrite(plot_table, file = xargs$tfst, sep = "\t")

fst_complete[, mid_pos_mb := ((start + end)/2)/1e+06]
#plot_table[, cutoff := as.character(cutoff)]
#fst_complete[, plot_cutoff := fcase(cutoff == 0.001, "0.1%",
#                                    cutoff == 0.005, "0.5%",
#                                    cutoff == 0.01, "1%",
#                                    cutoff == 100, "-")]

#saves all windows 
fst_pairs <- fst_complete$test |> unique()
#all_fst_windows_per_pair <- vector('list', length(fst_pairs))
for (i in seq_along(fst_pairs)) {
  all_fst <- 
    fst_complete[test == fst_pairs[[i]], .(chrom, start, end, FST, test, cutoff)]
  setnames(all_fst, colnames(all_fst), 
           c('Chrom', "Start", "End", "WindowFST", "pops", "cutoff"))
  setorder(all_fst)
  fwrite(all_fst, file = paste(xargs$tfst, fst_pairs[[i]], CHROM, '.tsv',
                               sep = '_'),
         sep = '\t')
  }



## stop
# local <- 
# str_split_i(string = file_name, pattern = "/", i = 7) |>
#   str_split_i(pattern = "_", i = 1)
# 
# size <- 
#   str_split_i(string = file_name, pattern = "/", i = 7) |>
#   str_split_i(pattern = "_", i = 2)
# 
#   
#   
# WINDOW_FST <- 
# fst_complete |>
#   ggplot(aes()) +
#   geom_point(aes(x = mid_pos_mb, y = FST, color = plot_cutoff),
#              size = 0.01) +
#   facet_grid(rows = vars(test),
#              cols = vars(chrom),
#              scales = "free_x") +
#   labs(#title = paste0(local, " Windown* FST"),
#        #caption =  paste0("*", size, " snp windows"),
#        x = "Genomic position in Mb") +
#   scale_color_manual(values = c("-" = "black",
#                                 "0.1%" = "red",
#                                 "0.5%" = "green",
#                                 "1%" = "blue")) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = -45, vjust = 0.1, hjust=0),
#         strip.background =element_rect(fill="white"),
#         strip.text =element_text(color="black"))
# 
# 
# jpeg(filename = xargs$FSTfig,
#      width = 25,
#      height = 25,
#      units = "cm",
#      res = 1200)
# WINDOW_FST
# 
# dev.off()
# 
# 
window_size <-
w_fst[, .(window_size_bp = end - start, chrom)][, .(mean = mean(window_size_bp),
                                                    max = max(window_size_bp),
                                                    min = min(window_size_bp),
                                                    sd = sd(window_size_bp)), by = chrom]


fwrite(x = window_size, file = xargs$wsize, sep = "\t")
