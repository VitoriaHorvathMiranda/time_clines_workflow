#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Separate pops from different times")
parser$add_argument('--NEtable', '-NE', help = 'outputs from NE.R')
parser$add_argument('--outputPath', '-oPh', 
                    help= 'the path to the outputs + prefix')

xargs<- parser$parse_args()
pops <- xargs$populations

#reads all pops
ALL_POP <- fread(file = xargs$NEtable)

#separate pops based on year
POP_97 <- ALL_POP[population %like% "97"]

POP_09_10 <- ALL_POP[population %like% "09" | population %like% "10"]

POP_17 <- ALL_POP[population %like% "17"]

POP_22_23 <- ALL_POP[population %like% "22" | population %like% "23"]

otp97 <- paste0(xargs$outputPath, "97.tsv")
otp0910 <- paste0(xargs$outputPath, "0910.tsv")
otp17 <- paste0(xargs$outputPath, "17.tsv")
otp2223 <- paste0(xargs$outputPath, "2223.tsv")

fwrite(POP_97, otp97, sep = "\t")
fwrite(POP_09_10, otp0910, sep = "\t")
fwrite(POP_17, otp17, sep = "\t")
fwrite(POP_22_23, otp2223,  sep = "\t")
