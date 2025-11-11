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
ALL_POP <- fread("/dados/time_clines/data/seqs/calls/NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
ALL_POP <- fread(file = xargs$NEtable)

#separate pops based on year
# 1997 -----------------------------------------------
POP_97 <- ALL_POP[population %like% "97"]

POP_97_noMFL <- ALL_POP[population %like% "97" & population != "MFL97"]

POP_97_noCMD97B <- ALL_POP[population %like% "97" & population != "CMD97B"]


# 2009/2010  -----------------------------------------------
#POP_09_10dlFL <- ALL_POP[(population %like% "09" | population %like% "10") & population != "JFL10"]
#POP_09_10 <- ALL_POP[(population %like% "09" | population %like% "10") & population != "dlFL10"]

POP_09_10_all <- ALL_POP[(population %like% "09" | population %like% "10")]

POP_09_10HFL_nodls_noJ <- ALL_POP[(population %like% "09" | population %like% "10") & 
                                 !(population %in% c("JFL10", "dlGA10", "dlSC10"))]


POP_09_10HJFL_nodls <- ALL_POP[(population %like% "09" | population %like% "10") & 
                                 !(population %in% c("dlGA10", "dlSC10"))]


otp97 <- paste0(xargs$outputPath, "97.tsv")
otp97FL <- paste0(xargs$outputPath, "97noMFL.tsv")
otp97cmd <- paste0(xargs$outputPath, "97noCMD97B.tsv")
otp0910all <- paste0(xargs$outputPath, "0910all.tsv")
otp0910FL2nodlJ <- paste0(xargs$outputPath, "0910FL2noDLnoJ.tsv")
otp0910onlyFL2 <- paste0(xargs$outputPath, "0910FL2noDL.tsv")



fwrite(POP_97, otp97, sep = "\t")
fwrite(POP_97_noMFL, otp97FL, sep = "\t")
fwrite(POP_97_noCMD97B, otp97cmd, sep = "\t")
fwrite(POP_09_10_all, otp0910all, sep = "\t")
fwrite(POP_09_10HFL_nodls_noJ, otp0910FL2nodlJ, sep = "\t")
fwrite(POP_09_10HJFL_nodls, otp0910onlyFL2, sep = "\t")



