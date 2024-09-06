#!/usr/bin/env Rscript

library(argparse)
library(data.table)
#library(tidyverse)

#parse arguments 
parser <- ArgumentParser(description= "")
parser$add_argument('--sync', '-sync', help= 'sync file with a lot of missing data')
parser$add_argument('--bs', '-bs', help= 'bed sites file')
parser$add_argument('--output', '-o', help= 'sync file with missing data from my samples fixed (but not from the painel!!)')
xargs<- parser$parse_args()

#all_sync <- fread("/dados/time_clines/analysis/ancestry/all_samples.sync")
all_sync <- fread(xargs$sync)
#setnames(all_sync, colnames(all_sync)[1:3], c("CHROM", "POS", "REF"))

# BS <- fread("importand_bed_sites.tsv",
#             col.names = c("CHROM", "POS", "bs"),
#             colClasses = c("character", "integer", "character"))

BS <- fread(xargs$bs,
            col.names = c("CHROM", "POS", "bs"),
            colClasses = c("character", "integer", "character"))


BS[, paste0("V", 4:22, "bs") := tstrsplit(bs, "", fixed = TRUE) ]

par_fixed <- merge.data.table(all_sync, BS, by = c("CHROM", "POS"), all.x = TRUE)

# Get the column names
V_cols <- paste0("V", 4:22)        # V4 to V22
Vbs_cols <- paste0("V", 4:22, "bs") # V4bs to V22bs

# Loop through each pair of V and Vbs columns
for (i in seq_along(V_cols)) {
  V_col <- V_cols[i]
  Vbs_col <- Vbs_cols[i]
  
  # Apply the logic to each pair
  par_fixed[, (V_col) := fcase(
    get(Vbs_col) == "0" & is.na(get(V_col)), "fill", # posição sem variante
    get(Vbs_col) == "0" & !is.na(get(V_col)), get(V_col), # posição com variante
    get(Vbs_col) == "1", "-", # posição sem info
    is.na(get(Vbs_col)) & is.na(get(V_col)), "fill", # posição sem variante
    is.na(get(Vbs_col)) & !is.na(get(V_col)), get(V_col) # posição com variante
  )]
  
  #sync files are A:T:C:G:N:del format
  # 23 is the mean depth
  par_fixed[, (V_col) := fcase(get(V_col) == "fill" & REF == "A", "23:0:0:0:0:0",
                          get(V_col) == "fill" & REF == "T", "0:23:0:0:0:0",
                          get(V_col) == "fill" & REF == "C", "0:0:23:0:0:0",
                          get(V_col) == "fill" & REF == "G", "0:0:0:23:0:0",
                          get(V_col) == "-", "0:0:0:0:0:0",
                     rep(TRUE, .N), as.character(get(V_col)))]
  
}


par_fixed <- par_fixed[, !c(Vbs_cols, "bs"), with = FALSE]

#df[df[, Reduce(`|`, lapply(.SD, `%like%`, "test")), .SDcols = cols]]
par_fixed <- par_fixed[!(par_fixed[, Reduce(`&`, lapply(.SD, `==`, "0:0:0:0:0:0")), .SDcols = V_cols])]


fwrite(par_fixed,
       file = xargs$output,
       sep = "\t", col.names = TRUE,
       na = "NA")
