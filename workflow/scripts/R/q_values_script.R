#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(qvalue)
library(patchwork)

#parse arguments
parser <- ArgumentParser(description= "compute qvalues for all p-values and plot quality plots")
parser$add_argument('--SNPs', '-snps', help = 'path to outputs from glm_script.R')
parser$add_argument('--output97', '-o97', help = 'table with all SNPs and respective q-values 97, .tsv')
parser$add_argument('--output0910', '-o0910', help = 'table with all SNPs and respective q-values 0910, .tsv')
parser$add_argument('--histPlot', '-hp', help = 'p-value histogram output, .jpeg')
parser$add_argument('--qPlot97', '-qp97', help = 'p-plot output 97, .jpeg')
parser$add_argument('--qPlot0910', '-qp0910', help = 'p-plot output 0910, .jpeg')
xargs<- parser$parse_args()

files_97 <- list.files(xargs$SNPs,
                       full.names = TRUE, pattern = "*97_.{1,2}.tsv")

files_0910 <- list.files(xargs$SNPs,
                         full.names = TRUE, pattern = "*0910_.{1,2}.tsv")

results97 <- lapply(files_97, fread)
results0910 <- lapply(files_0910, fread)


#adiciona uma coluna com os cromossomos
chrom <- c("2L","2R","3L","3R","X")
add_chrom <- function(data, chrom){
  for (i in seq_along(data)) {
    data[[i]][, chrom := chrom[[i]]]
  }
  
}
add_chrom(results97, chrom)
add_chrom(results0910, chrom)

# junta os chroms
results97 <- rbindlist(results97)
results0910 <- rbindlist(results0910)

fdr_cutoffs <- function(data) {
  p_values <- data$p_value
  qobj <- qvalue(p = p_values)
  data[, qvalue := qobj$qvalues]
  qobj #o output é o qobj
  
}

qobj97 <- fdr_cutoffs(results97)
qobj0910 <- fdr_cutoffs(results0910)

fwrite(results97, file = xargs$output97,sep = "\t")
fwrite(results0910, file = xargs$output0910,sep = "\t")

jpeg(filename = xargs$histPlot,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

hist(qobj97) / hist(qobj0910)

dev.off()


jpeg(filename = xargs$qPlot97,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj97) 

dev.off()


jpeg(filename = xargs$qPlot0910,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj0910) 

dev.off()


