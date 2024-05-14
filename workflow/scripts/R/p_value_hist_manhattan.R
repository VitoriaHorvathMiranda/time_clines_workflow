#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(qvalue)
library(tidyverse)


#parse arguments 
parser <- ArgumentParser(description= "Makes hist of p-values and manhattan plots")
parser$add_argument('--P97Path', '-p97', help= 'Path + prefix of glm_script output from 1997')
parser$add_argument('--P0910Path', '-p0910', help= 'Path + prefix of glm_script output from 2009/2010')
parser$add_argument('--hist', '-ht', help= 'p-values histogram')
parser$add_argument('--m97', '-m97', help= 'Manhattan plot 1997')
parser$add_argument('--m0910', '-m0910', help= 'Manhattan plot 2009/2010')
xargs<- parser$parse_args()

input97 <- xargs$P97Path
input0910 <- xargs$P0910Path
pattern_97 <- str_split_i(input97, "/", i=6)
pattern_0910 <- str_split_i(input0910, "/", i=6)
path <- str_split_i(input97, "/p", 1)

path_97 <- list.files(path,
           pattern = pattern_97,
           full.names = TRUE)
path_0910 <- list.files(path,
                        pattern = pattern_0910,
                        full.names = TRUE)

p_values_97 <- lapply(path_97, fread)
p_values_0910 <- lapply(path_0910, fread)

chrom <- c("2L","2R","3L","3R","X")

for (i in seq_along(chrom)) {
  p_values_97[[i]] <- p_values_97[[i]][, CHROM := chrom[[i]]]
}

for (i in seq_along(chrom)) {
  p_values_0910[[i]] <- p_values_0910[[i]][, CHROM := chrom[[i]]]
}

p_values_97 <- rbindlist(p_values_97)
p_values_0910 <- rbindlist(p_values_0910)

p_values_97[, year := "1997"]
p_values_0910[, year := "2009/2010"]

all_ps <- merge.data.table(p_values_97, p_values_0910, all = TRUE)

HIST <- 
all_ps |>
  ggplot() +
  geom_histogram(aes(p_value), boundary = 0, binwidth = 0.02) +
  facet_wrap(vars(year)) +
  labs(x = "p-value",
       caption = "min-freq = 0.001, min-cov = 15, min-count = 5") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14))


######get FDR cutof 1997 ---------------
Pvalue_97 <- p_values_97$p_value #all_ps[year == "1997"]$p_value
qobj_97 <- qvalue(p = Pvalue_97)

#cutoff 10%
p_value97_cutoff_0.1 <- max(qobj_97$pvalues[qobj_97$qvalues <= 0.1])

#cutoff 0.05%
p_value97_cutoff_0.05 <- max(qobj_97$pvalues[qobj_97$qvalues <= 0.05])

#cutoff 1%
p_value97_cutoff_0.01 <- max(qobj_97$pvalues[qobj_97$qvalues <= 0.01])

######get FDR cutof 2009/2010 ---------------
Pvalue_0910 <- p_values_0910$p_value #all_ps[year == "2009/2010"]$p_value
qobj_0910 <- qvalue(p = Pvalue_0910)

#cutoff 10%
p_value0910_cutoff_0.1 <- max(qobj_0910$pvalues[qobj_0910$qvalues <= 0.1])

#cutoff 0.05%
p_value0910_cutoff_0.05 <- max(qobj_0910$pvalues[qobj_0910$qvalues <= 0.05])

#cutoff 1%
p_value0910_cutoff_0.01 <- max(qobj_0910$pvalues[qobj_0910$qvalues <= 0.01])

####### 1997
setkey(p_values_97, CHROM, POS)

axisdf_97 <-
  p_values_97[, ID := .I] |> 
  group_by(CHROM) |> 
  summarize(center=( max(ID) + min(ID) ) / 2 )

MANHATTAN_97 <- 
p_values_97 |>
  ggplot(aes(x = ID, y =-log10(p_value))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.7) +
  geom_hline(aes(yintercept = -log10(p_value97_cutoff_0.1)), color = "darkgreen") +
  geom_hline(aes(yintercept = -log10(p_value97_cutoff_0.01)), color = "blue") +
  geom_hline(aes(yintercept = -log10(p_value97_cutoff_0.05)), color = "red") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value97_cutoff_0.1) + 0.25,
                label = "FDR 10%"),
            size = 2.5, color = "green") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value97_cutoff_0.01) + 0.25,
                label = "FDR 1%"),
            size = 2.5, color = "blue") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value97_cutoff_0.05) + 0.25,
                label = "FDR 5%"),
            size = 2.5, color = "red") +
  scale_color_manual(values = rep(c("lightblue3", "lightcoral"), 5)) +
  scale_x_continuous(label = axisdf_97$CHROM, breaks= axisdf_97$center) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,16.3)) +
  theme_light() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 15)
  ) + 
  labs(title = "1997", y = "-log10(p-value)", x = "Chromosome")

####### 2009/2010
setkey(p_values_0910, CHROM, POS)

axisdf_0910 <-
  p_values_0910[, ID := .I] |> 
  group_by(CHROM) |> 
  summarize(center=( max(ID) + min(ID) ) / 2 )


MANHATTAN_0910 <- 
p_values_0910 |>
  ggplot(aes(x = ID, y =-log10(p_value))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.7) +
  geom_hline(aes(yintercept = -log10(p_value0910_cutoff_0.1)), color = "darkgreen") +
  geom_hline(aes(yintercept = -log10(p_value0910_cutoff_0.01)), color = "blue") +
  geom_hline(aes(yintercept = -log10(p_value0910_cutoff_0.05)), color = "red") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value0910_cutoff_0.1) + 0.25,
                label = "FDR 10%"),
            size = 2.5, color = "green") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value0910_cutoff_0.01) + 0.25,
                label = "FDR 1%"),
            size = 2.5, color = "blue") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value0910_cutoff_0.05) + 0.25,
                label = "FDR 5%"),
            size = 2.5, color = "red") +
  scale_color_manual(values = rep(c("lightblue3", "lightcoral"), 5)) +
  scale_x_continuous(label = axisdf_0910$CHROM, breaks= axisdf_0910$center) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,16.3)) +
  theme_light() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 15)
  ) + 
  labs(title = "2009/2010", y = "-log10(p-value)", x = "Chromosome")


#### save images -----------------------------------------------------

jpeg(filename = xargs$m97,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
MANHATTAN_97

dev.off()

jpeg(filename = xargs$m0910,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
MANHATTAN_0910

dev.off()


jpeg(filename = xargs$hist,
     width = 26,
     height = 14,
     units = "cm",
     res = 1200)
HIST

dev.off()
