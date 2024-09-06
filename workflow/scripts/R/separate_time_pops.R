#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Separate pops from different times")
parser$add_argument('--NEtable', '-NE', help = 'outputs from NE.R')
parser$add_argument('--outputPath', '-oPh', 
                    help= 'the path to the outputs + prefix')

xargs<- parser$parse_args()
#pops <- xargs$populations

#reads all pops
ALL_POP <- fread(file = xargs$NEtable)

#separate pops based on year
POP_97 <- ALL_POP[population %like% "97"]

POP_09_10 <- ALL_POP[(population %like% "09" | population %like% "10") & population != "dlFL10"]

POP_09_10dlFL <- ALL_POP[(population %like% "09" | population %like% "10") & population != "JFL10"]

POP_09_10dlFL_nodls <- ALL_POP[(population %like% "09" | population %like% "10") & 
                                 !(population %in% c("JFL10", "dlGA10", "dlSC10"))]

POP_97_noMFL <- ALL_POP[population %like% "97" & population != "MFL97"]

POP_09_10_all <- ALL_POP[(population %like% "09" | population %like% "10")]


otp97 <- paste0(xargs$outputPath, "97.tsv")
otp0910 <- paste0(xargs$outputPath, "0910.tsv")
otp0910FL2 <- paste0(xargs$outputPath, "0910FL2.tsv")
otp0910onlyFL2 <- paste0(xargs$outputPath, "0910FL2noDL.tsv")
otp97FL <- paste0(xargs$outputPath, "97noMFL.tsv")
otp0910all <- paste0(xargs$outputPath, "0910all.tsv")


fwrite(POP_97, otp97, sep = "\t")
fwrite(POP_09_10, otp0910, sep = "\t")
fwrite(POP_09_10dlFL, otp0910FL2, sep = "\t")
fwrite(POP_09_10dlFL_nodls, otp0910onlyFL2, sep = "\t")
fwrite(POP_97_noMFL, otp97FL, sep = "\t")
fwrite(POP_09_10_all, otp0910all, sep = "\t")

