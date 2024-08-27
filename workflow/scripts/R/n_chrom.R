#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "gets latitude and effective Nº of chrom")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--DepthFreq', '-df', help= 'output from freq_extraction.R')
parser$add_argument('--output', '-o', help= 'tsv with depth, latitude and n of chrom for each snp')
xargs<- parser$parse_args()

#get data
#metadata <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
metadata <- fread(file = xargs$metadata)
#depths <- fread("/dados/time_clines/data/seqs/calls/freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
depths <- fread(file = xargs$DepthFreq)

depths <- depths[depth != ".", ]
depths[, depth := as.double(depth)]

#get relevant columns 
I_info <- metadata[, .(population, n_females, fly_sex, flies_per_lin, latitude)]

#compute n of chrom
# #multiply by 1.5 the number of females because we got 2 f1 daughters per line
I_info <- I_info[, n_chrom := fcase(flies_per_lin == 2, n_females*1.5,
                                    flies_per_lin == 1, n_females*2)]#[,!c("n_females")]

#join 
full <- merge.data.table(I_info, depths, by = "population")
full[, n_chrom := ifelse(CHROM == "X" & fly_sex == "male", n_chrom/2, n_chrom)]
full <- full[,!c("n_females", "fly_sex", "flies_per_lin")]

#calcula o NE 
full[, NE := ((1/depth) + (1/n_chrom))^-1]

#save
output <- xargs$output
write.table(full, file = output, sep = "\t", row.names = FALSE)



