#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(qvalue)
library(tidyverse)


#parse arguments 
parser <- ArgumentParser(description= "Makes hist of p-values and manhattan plots")
parser$add_argument('--path', '-path', help= 'Path + prefix of glm_script outputs')
parser$add_argument('--Opath', '-Opath', help= 'Path and prefix to q-value output')
parser$add_argument('--histPlots', '-qhist', help = 'p-value histogram output from qvalue package, .jpeg')
parser$add_argument('--hist', '-hist', help= 'p-values histogram, .jpeg')
parser$add_argument('--m97', '-m97', help= 'Manhattan plot 1997')
parser$add_argument('--m0910', '-m0910', help= 'Manhattan plot 2009/2010')
parser$add_argument('--m0910FL2', '-m0910FL2', help= 'Manhattan plot 2009/2010 with FL2')
parser$add_argument('--qPlot97', '-qp97', help = 'p-plot output 97, .jpeg')
parser$add_argument('--qPlot0910', '-qp0910', help = 'p-plot output 0910, .jpeg')
parser$add_argument('--qPlot0910FL2', '-qp0910FL2', help = 'p-plot output 0910, .jpeg')
xargs<- parser$parse_args()


input_general <- xargs$path
input97 <- paste0(input_general, "97")
input0910 <- paste0(input_general, "0910")
input0910FL2 <- paste0(input_general, "0910FL2")

pattern_97 <- str_split_i(input97, "/", i=6)
pattern_0910 <- str_split_i(input0910, "/", i=6)
pattern_0910FL2 <- str_split_i(input0910FL2, "/", i=6)

path <- str_split_i(input97, "/p", 1)

path_97 <- list.files(path,
           pattern = pattern_97,
           full.names = TRUE)
path_0910 <- list.files(path,
                        pattern = pattern_0910,
                        full.names = TRUE)

path_0910FL2 <- list.files(path,
                        pattern = pattern_0910FL2,
                        full.names = TRUE)

p_values_97 <- lapply(path_97, fread)
p_values_0910 <- lapply(path_0910, fread)
p_values_0910FL2 <- lapply(path_0910FL2, fread)

chrom <- c("2L","2R","3L","3R","X")

for (i in seq_along(chrom)) {
  p_values_97[[i]] <- p_values_97[[i]][, CHROM := chrom[[i]]]
}

for (i in seq_along(chrom)) {
  p_values_0910[[i]] <- p_values_0910[[i]][, CHROM := chrom[[i]]]
}

for (i in seq_along(chrom)) {
  p_values_0910FL2[[i]] <- p_values_0910FL2[[i]][, CHROM := chrom[[i]]]
}

p_values_97 <- rbindlist(p_values_97)
p_values_0910 <- rbindlist(p_values_0910)
p_values_0910FL2 <- rbindlist(p_values_0910FL2)

p_values_97[, year := "1997"]
p_values_0910[, year := "2009/2010"]
p_values_0910FL2[, year := "2009/2010-FL2"]

all_ps <- merge.data.table(p_values_97, p_values_0910, all = TRUE)
all_ps <- merge.data.table(all_ps, p_values_0910FL2, all = TRUE)

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

######get FDR cutof FL2 ---------------
Pvalue_FL2 <- p_values_0910FL2$p_value #all_ps[year == "1997"]$p_value
qobj_FL2 <- qvalue(p = Pvalue_FL2)

#cutoff 10%
p_valueFL2_cutoff_0.1 <- max(qobj_FL2$pvalues[qobj_FL2$qvalues <= 0.1])

#cutoff 0.05%
p_valueFL2_cutoff_0.05 <- max(qobj_FL2$pvalues[qobj_FL2$qvalues <= 0.05])

#cutoff 1%
p_valueFL2_cutoff_0.01 <- max(qobj_FL2$pvalues[qobj_FL2$qvalues <= 0.01])


###### adds q-value column
p_values_0910FL2[, qvalue := qobj_FL2$qvalues]
p_values_0910[, qvalue := qobj_0910$qvalues]
p_values_97[, qvalue := qobj_97$qvalues]

out0910 <- paste0(xargs$Opath, "0910.tsv")
out0910FL2 <- paste0(xargs$Opath, "0910FL2.tsv")
out97 <- paste0(xargs$Opath, "97.tsv")

fwrite(p_values_0910FL2, file = out0910FL2, sep = "\t")
fwrite(p_values_0910, file = out0910, sep = "\t")
fwrite(p_values_97, file = out97, sep = "\t")

#### q-plots -----------------------------------------------

jpeg(filename = xargs$histPlots,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

hist(qobj_97) / hist(qobj_0910) / hist(qobj_FL2)

dev.off()

jpeg(filename = xargs$qPlot97,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj_97) 

dev.off()


jpeg(filename = xargs$qPlot0910,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj_0910) 

dev.off()

jpeg(filename = xargs$qPlot0910FL2,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj_FL2) 

dev.off()

################ Manhattan plots ------------------------------------
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
            size = 2.5, color = "darkgreen") +
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
            size = 2.5, color = "darkgreen") +
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

####### 2009/2010 - FL2
setkey(p_values_0910FL2, CHROM, POS)

axisdf_FL2 <-
  p_values_0910FL2[, ID := .I] |> 
  group_by(CHROM) |> 
  summarize(center=( max(ID) + min(ID) ) / 2 )


MANHATTAN_FL2 <- 
  p_values_0910FL2 |>
  ggplot(aes(x = ID, y =-log10(p_value))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.7) +
  geom_hline(aes(yintercept = -log10(p_valueFL2_cutoff_0.1)), color = "darkgreen") +
  geom_hline(aes(yintercept = -log10(p_valueFL2_cutoff_0.01)), color = "blue") +
  geom_hline(aes(yintercept = -log10(p_valueFL2_cutoff_0.05)), color = "red") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_valueFL2_cutoff_0.1) + 0.25,
                label = "FDR 10%"),
            size = 2.5, color = "darkgreen") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_valueFL2_cutoff_0.01) + 0.25,
                label = "FDR 1%"),
            size = 2.5, color = "blue") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_valueFL2_cutoff_0.05) + 0.25,
                label = "FDR 5%"),
            size = 2.5, color = "red") +
  scale_color_manual(values = rep(c("lightblue3", "lightcoral"), 5)) +
  scale_x_continuous(label = axisdf_FL2$CHROM, breaks= axisdf_FL2$center) +
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
  labs(title = "2009/2010 - FL2", y = "-log10(p-value)", x = "Chromosome")


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


jpeg(filename = xargs$m0910FL2,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
MANHATTAN_FL2

dev.off()

jpeg(filename = xargs$hist,
     width = 26,
     height = 14,
     units = "cm",
     res = 1200)
HIST

dev.off()
