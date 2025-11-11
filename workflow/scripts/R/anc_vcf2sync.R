#!/usr/bin/env Rscript
#como o vcf dos dos "ancestrais" é individual é ruim usar um script pronto para criar o .sync
library(argparse)
library(data.table)
library(tidyverse)


#parse arguments 
parser <- ArgumentParser(description= "transformns ancestral pops vcf into sync")
parser$add_argument('--vcf', '-vcf', help= 'clean vcf with ancetral pops')
parser$add_argument('--chrom', '-chrom', help = "chromossome number to filter data")
parser$add_argument('--mask', '-m', help = "Regions with admixture in westAFR pop to mask")
parser$add_argument('--output', '-o', help= 'sync file with ancetral pops')
xargs<- parser$parse_args()

# chrom <- "2L"
# vcf <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
#                          "/dados/time_clines/data/database_seqs/call/database_all_zi_eg/zi_eg_eu_westafr_2L_biallelic_clean.vcf"),
#              drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT", "EG31nWGS"))


#admix_westAFR <- fread("~/time_clines_workflow/resources/Admixture_segments_westAFR_R6.txt")

admix_westAFR <- fread(xargs$mask)
admix_westAFR[, start := Start]
admix_westAFR[, end := End]


chrom <- as.character(xargs$chrom)
chrom
paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
                         xargs$vcf)
vcf <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
                         xargs$vcf),
             drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT", "EG31nWGS"))
#tirando "EG31nWGS" porque tem muito IBD com "EG26nWGS"

head(vcf)
setnames(vcf, "#CHROM", "CHROM")

ancestries <- c("EU", "westAFR", "ZI", "EG")
bases <- c("A", "T", "C", "G")
#vcf <- vcf[1:100000]
# adds a column that divides vcf into 10 parts so it is possible to loop over smaller parts
vcf[, chrom_frac := rep(1:10, each = ceiling(nrow(vcf)/10))[1:nrow(vcf)]]


sync <- vector("list", length = 10)
for (i in seq_along(1:10)) {
  vcf_tidy <- melt.data.table(vcf[chrom_frac == i],
                              id.vars = c(1:4, length(vcf)),
                              variable.name = "sample",
                              value.name = "genotype")
  vcf_tidy[, pos := POS]
  
  vcf_tidy <- 
  admix_westAFR[vcf_tidy, on = .(chrom = CHROM,
                                 Line = sample,
                                 Start <= pos,
                                 End >= pos)]
  
  vcf_tidy[!is.na(start), genotype := ".:.:.:.:."]
  
  vcf_tidy <- 
  vcf_tidy[, !c("Start", "start", "End", "end"), with = FALSE]
  
  setnames(vcf_tidy, c("Line", "chrom"), c("sample", "CHROM"))
  
  vcf_tidy[, c("GT","AD","DP","GQ","PL") := 
             tstrsplit(genotype, ":", fixed = TRUE)]
  
  #GT = genotype
  #DP = Depth
  vcf_tidy[, DP := as.integer(DP)]
  vcf_tidy <- vcf_tidy[DP > 2]
  
  vcf_tidy[, ancestry := fcase(sample %like% "FR|SU", "EU",
                               sample %like% "CO|NG|GU|GA", "westAFR",
                               sample %like% "ZI", "ZI",
                               sample %like% "EG", "EG")]
  
  #vcf_tidy <- vcf_tidy[!(GT %in% c("0/1", "0|1"))]
  
  vcf_tidy[, alt_gt := fcase(GT == 1, 1,
                             GT == 0, 0,
                             GT ==  "0/0", 0,
                             GT == "1|1", 1,
                             GT == "1/1", 1,
                             GT == "0|0", 0,
                             GT == "0/1", sample(c(0,1), 1),
                             GT == "0|1", sample(c(0,1), 1))]
  
  vcf_tidy[, ref_gt := fcase(GT == 1, 0,
                             GT == 0, 1,
                             GT ==  "0/0", 1,
                             GT == "1|1", 0,
                             GT == "1/1", 0,
                             GT == "0|0", 1,
                             GT == "0/1" & alt_gt == 1, 0,
                             GT == "0/1" & alt_gt == 0, 1,
                             GT == "0|1" & alt_gt == 1, 0,
                             GT == "0|1" & alt_gt == 0, 1,
                             default = NA)]
  
  painel_count <- 
    vcf_tidy[, .(alt_allele_count = sum(alt_gt),
                 ref_allele_count = sum(ref_gt)),
             by = c("CHROM", "POS", "REF", "ALT", "ancestry")]
  
  ref_count <- lapply(ancestries, function(x) {
    dcast.data.table(painel_count[ancestry == x], CHROM + POS + REF ~ REF,
                     value.var = "ref_allele_count",
                     fill = 0)}
  )
  
  alt_count <- lapply(ancestries, function(x) {
    dcast.data.table(painel_count[ancestry == x], CHROM + POS + ALT ~ ALT,
                     value.var = "alt_allele_count",
                     fill = 0)}
  )
  
  sync_list <- 
  lapply(1:4, function(x) {
    all_count <- merge.data.table(alt_count[[x]], ref_count[[x]],
                                  by = c("CHROM", "POS"))
    
    all_count[, (bases) := lapply(bases, function(base) {
      get(paste0(base, ".x")) + get(paste0(base, ".y"))
    })]
    
    all_count <- all_count[, !(str_subset(colnames(all_count), ".x|.y")), with = FALSE]
    
    all_count[, ancestries[[x]] := do.call(paste, c(.SD, "0", "0", sep = ":")), .SDcols = bases]
    
    all_count[, !bases, with = FALSE]
    
  })
  
  sync[[i]] <- reduce(sync_list, merge.data.table, by = c("CHROM", "POS", "ALT", "REF"))
  print(sync[[i]])
}

sync <- rbindlist(sync)

fwrite(sync, xargs$output, sep = "\t")

