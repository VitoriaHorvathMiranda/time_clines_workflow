#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)

#parse arguments
parser <- ArgumentParser(description= "")
parser$add_argument('--vcf', '-vcf',
                    help = 'annotated vcf')
parser$add_argument('--effects', '-ef',
                    help = 'table with sequence ontology and effects')
parser$add_argument('--qValue', '-qv',
                    help = 'table with q-values for each tested snp')
parser$add_argument('--introns', '-i',
                    help = 'bed file with short introns')
parser$add_argument('--rec', '-r',
                    help = 'bed file with genome-wide recombination rate')
parser$add_argument('--output', '-o',
                    help = "table (tsv) with snps positions, q-values
                    and it's annotation")
xargs<- parser$parse_args()

# effects <- fread("~/time_clines_workflow/resources/effects.tsv")
# VCF <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.ann.vcf",
#               select = c(1,2,5,8))
# short_introns <- fread("/home/vitoria/time_clines_workflow/resources/dmel-6.32_neutral_short_introns.bed")
# 
# q_values <- fread("/dados/time_clines/analysis/time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910all.tsv")
#rec_rate <- fread("~/time_clines_workflow/resources/Dmel_recombination_rate_R6.bed")


effects <- fread(xargs$effects)
VCF <- fread(xargs$vcf,
             select = c(1,2,5,8))

short_introns <- fread(xargs$introns)

q_values <- fread(xargs$qValue)

rec_rate <- fread(xargs$rec)

setnames(VCF, "#CHROM", "chrom")
setnames(q_values, "CHROM", "chrom")

rec_rate <- rec_rate[rec_rate >= 3]
rec_rate[, start := as.integer(start)]
# filter for biallelic snps
VCF <- VCF[!(ALT %like% ","),]

VCF[, ANN := gsub(".*ANN=[A-Z].(.+?)\\|.*", "\\1", INFO)]
# VCF[, ANN2 := tstrsplit(INFO, ",", keep = 2)]
# VCF[, ANN2 := tstrsplit(ANN2, "\\|", keep = 2)]
# VCF[, ANN3 := tstrsplit(INFO, ",", keep = 3)]
# VCF[, ANN3 := tstrsplit(ANN3, "\\|", keep = 2)]


VCF[, impact := str_split_i(INFO, pattern = "\\|", i = 3)]
VCF <- VCF[, !"INFO", with = FALSE]

VCF_effects <- merge.data.table(VCF, effects, by = "ANN")

short_introns[, start := neutral_start]
short_introns[, end := neutral_end]

all_effects <- 
  short_introns[VCF_effects, on = .(chrom,
                                    start <= POS,
                                    end >= POS)]#[!is.na(start)]

all_effects <- all_effects[, !"end"]
setnames(all_effects, "start", "POS")

all_effects <- rec_rate[all_effects, on = .(chrom,
                                            start <= POS,
                                            end >= POS)]#[!is.na(rec_rate) &
                                          #  !is.na(neutral_start) &
                                          #  ANN == 'intron_variant']


all_effects[!is.na(neutral_start) & !is.na(rec_rate), effect := ifelse(ANN == 'intron_variant',
                                                                       "SHORT_INTRON",
                                                                       effect)]

#all_effects[effect == "SHORT_INTRON"]

all_effects <- all_effects[, !c("neutral_start", "neutral_end", "start", "rec_rate"), with = FALSE]
setnames(all_effects, "end", "POS")

# test2 <-
#   short_introns[VCF, on = .(chrom,
#                             start <= POS,
#                             end >= POS)][!is.na(neutral_start)]
# 
# test2[, .(nrow(.SD)), by = ANN][ANN != "intron_variant", !"ANN", with = FALSE] |> sum()
# test2[, .(nrow(.SD)), by = c("ANN","ANN2", "ANN3")] |> View()
# test3 <- test2[ANN == "intron_variant"]
# rec_rate[, start := as.integer(start)]
# 
# rec_rate[rec_rate >= 3][test3, on = .(chrom,
#                                       start <= start,
#                                       end >= start)][!is.na(rec_rate)]


q_value_eff <- merge.data.table(q_values, all_effects, by = c("POS", "chrom"))


fwrite(q_value_eff, file = xargs$output, sep = "\t")



