#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(viridis)
library(ggrepel)

#parse arguments 
parser <- ArgumentParser(description= "transformns ancestral pops vcf into sync")
parser$add_argument('--vcf', '-vcf', help= 'clean vcf with ancetral pops')
parser$add_argument('--chrom', '-chrom', help = "chromossome number to filter data")
parser$add_argument('--sync', '-sync', help = "filtered sync file with ancetral pops")
parser$add_argument('--inv', '-inv', help = "table with inversions")
#parser$add_argument('--admix', '-admix', help = "table with admixed regions in westAFR")
parser$add_argument('--pca', '-pca', help= 'pca plot')
xargs<- parser$parse_args()


#read data ---------------------------------------------------------------------
# chrom <- "3R"
# vcf <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
#                          "/dados/time_clines/data/database_seqs/call/database_all_zi_eg/zi_eg_eu_westafr_3R_biallelic_clean.vcf"),
#              drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT"),
#              nrows = 100000)
# admix_westAFR <- fread("~/time_clines_workflow/resources/Admixture_segments_westAFR_R6.txt")
# inv <- fread("~/time_clines_workflow/resources/inversions_anc_samples.csv")

#admix_westAFR <- fread(xargs$admix)
inv <- fread(xargs$inv)
chrom <- xargs$chrom
vcf <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
                          xargs$vcf),
             drop = c("ID", "QUAL", "FILTER", "INFO", "FORMAT"))

#setnames(admix_westAFR, c("chrom"), c("CHROM"))
setnames(vcf, "#CHROM", "CHROM")

# sync <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
#                            "/dados/time_clines/data/database_seqs/call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean_minCount2_minDepth10_totalMAF0.05sync.sync"))

sync <- fread(cmd = paste0("awk '$1 == \"#CHROM\" || $1 == \"", chrom, "\"' ",
                           xargs$sync))

setnames(sync, colnames(sync), c("CHROM", "POS", "REF", "EU", "westAFR", "ZI", "EG"))


#filter vcf based on filered .sync ---------------------------------------------
#this sync was filtered with grenedalf
sync <- 
sync[ZI != ".:.:.:.:.:." & EG != ".:.:.:.:.:." & westAFR != ".:.:.:.:.:." & EU != ".:.:.:.:.:."]

vcf <- 
vcf[sync[, .(CHROM, POS)], on = .(CHROM, POS), nomatch = NULL]

# vcf <- 
# vcf[!admix_westAFR, on = .(CHROM, POS >= Start, POS <= End)]

#Get genotype ------------------------------------------------------------------
# adds a column that divides vcf into 20 parts so it is possible to loop over smaller parts
vcf[, chrom_frac := rep(1:20, each = ceiling(nrow(vcf)/20))[1:nrow(vcf)]]

pre_pca_table <- vector("list", length = 20)
for (i in seq_along(1:20)) {
  vcf_tidy <- melt.data.table(vcf[chrom_frac == i],
                              id.vars = c(1:4, length(vcf)),
                              variable.name = "sample",
                              value.name = "genotype")
  
  vcf_tidy[, c("GT","AD","DP","GQ","PL") := 
             tstrsplit(genotype, ":", fixed = TRUE)]
  
  #GT = genotype
  #DP = Depth
  vcf_tidy[, DP := as.integer(DP)]
  vcf_tidy <- vcf_tidy[DP > 2]
  
  vcf_tidy[, gt := fcase(GT == 1, 1,
                         GT == 0, 0,
                         GT ==  "0/0", 0,
                         GT == "1|1", 1,
                         GT == "1/1", 1,
                         GT == "0|0", 0,
                         GT == "0/1", sample(c(0,1), 1),
                         GT == "0|1", sample(c(0,1), 1))]
  
  vcf_tidy[, position2 := paste(CHROM, POS, sep = ":")]
  
  vcf_tidy <- 
    vcf_tidy[, .(position2, sample, gt)]
  
  pre_pca_table[[i]] <- dcast.data.table(vcf_tidy, position2 ~ sample,
                                         value.var = "gt")
  
}

pre_pca_table <- 
rbindlist(pre_pca_table)

#pca ---------------------------------------------------------------------------
fun_pca <- 
  function(table){
    no_NA_table <- na.omit(table)
    pca_table <- no_NA_table |> 
      column_to_rownames(var="position2")
    
    pca_table <- t(pca_table)
    prcomp(pca_table)
  }

pca_all <- fun_pca(pre_pca_table)

pca_var <- pca_all$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

pcs <- rownames_to_column(as.data.frame(pca_all$x), var = "sample")
pcs <- as.data.table(pcs)


#plot pca ----------------------------------------------------------------------
pcs_main <- pcs[, c("sample", "PC1", "PC2", "PC3", "PC4")][
  , population := fcase(sample %like% "FR|SU", "EU",
                        sample %like% "CO|NG|GU|GA", "westAFR",
                        sample %like% "ZI", "ZI",
                        sample %like% "EG", "EG")]

inv[, sample := ifelse(sample %like% "EG", 
                       str_replace(sample, "N", "nWGS"),
                       sample)]

test <- 
inv[pcs_main, on = .(sample)][
  , c("2L", "2R", "3L", "3R") := lapply(.SD, function(x){
    ifelse(is.na(x), "ST", x)
  }), .SDcols = c("2L", "2R", "3L", "3R")]

test[, `2L` := ifelse(sample %in% c("SU08", "SU37n", "SU29", "SU26n"),
                      "In(2L)t", `2L`)]

test[, X := "X"]

nshape <- test[, chrom, with = FALSE] |> unique() |> nrow()

PC1x2_ALL <- 
  test |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = population, shape = get(chrom)), size = 2, alpha = 0.7) +
  scale_shape_manual(values=1:nshape) +
  geom_text_repel(aes(label = sample), size = 1,
                  arrow = arrow(length = unit(0, "npc"), type = "open"),
                  segment.size = 0.2,       # Adjusts the line width
                  segment.linetype = "solid", # Ensures lines are solid
                  force = 5,               # Strengthens repulsion
                  max.overlaps = Inf,      # Prevents removal due to overlaps
                  box.padding = 0.5,       # Space around the label box
                  point.padding = 0.3,      # Space around the data point
                  #ylim = c(0.2, 0.5),
                  min.segment.length = 0, 
                  seed = 15) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  labs(title = chrom, shape = "Inversion") +
  theme_minimal()



pdf(file = xargs$pca,
    width = 7,
    height = 4)

PC1x2_ALL

dev.off()


