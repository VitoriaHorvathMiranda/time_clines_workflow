#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(biomaRt)

#parse arguments

parser <- ArgumentParser(description= "gets association GO file")
parser$add_argument('--dmelgtf', '-gtf', help = 'gtf file from flybase')
parser$add_argument('--output', '-o', help = 'table with go assossiations')
xargs<- parser$parse_args()

# parse gtf -----------------------------------------------------------------
dmel_R6.55_gtf <- fread(xargs$dmelgtf,
                        drop = c(2, 6:8),
                        col.names = c("chrom", "function", "start",
                                      "end", "gene_id"))


dmel_R6.55_gtf[, c("gene_ID","gene_symbol",
                   "transcript_id", "transcript_symbol", "a") := 
                 tstrsplit(gene_id, ";", fixed = TRUE)]

dmel_R6.55_gtf[, gene_ID := str_split_i(gene_ID, pattern = " ", i = 2)]

dmel_R6.55_gtf[, c("gene_symbol","transcript_id", "transcript_symbol") := 
                 lapply(.SD, function(x) str_split_i(x, pattern = " ", i = 3)),
               .SDcols = c("gene_symbol",
                           "transcript_id", "transcript_symbol")]

dmel_R6.55_gtf <- dmel_R6.55_gtf[, !c("gene_id", "a"), with = FALSE]
dmel_R6.55_gtf <- dmel_R6.55_gtf[chrom %in% c("2L", "2R", "3L", "3R", "X")]

dmel_R6.55_gtf[, gene_ID := str_remove_all(gene_ID, "\"")]

vector_gene_IDs <- unique(dmel_R6.55_gtf$gene_ID)


# get database ------------------------------------------------------------

ensembl <- useEnsembl(biomart = "genes")
searchDatasets(mart = ensembl, pattern = "melanogaster")
ensembl <- useDataset(dataset = "dmelanogaster_gene_ensembl", mart = ensembl)

searchAttributes(mart = ensembl, pattern = "GO term")

GOs <- getBM(attributes = c("ensembl_gene_id","go_id","name_1006"),
             filters = "ensembl_gene_id",
             values = vector_gene_IDs,
             mart = ensembl)

head(GOs)

#data.table function for nesting:
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

Gos_nested <- as.data.table(GOs)[, group_nest_dt(.SD, go_id, name_1006)]


#all <- overlap2_gtf[Gos_nested, on = c("gene_ID" = "ensembl_gene_id")]

Gos_nested[, gene_id := lapply(data, function(x) paste(x$ensembl_gene_id, collapse = " "))]
Gos_nested[, gene_id := as.character(gene_id)]
# Gos_nested[, GO_terms := lapply(data, function(x) paste(x$name_1006, collapse = ""))]
# Gos_nested[, GO_terms := as.character(GO_terms)]

Gos_nested <- Gos_nested[, !("data"), with = FALSE]

glimpse(Gos_nested)
Gos_nested <- Gos_nested[!(go_id == "")]

fwrite(Gos_nested, file = xargs$output,
       sep = "\t",
       col.names = FALSE)


