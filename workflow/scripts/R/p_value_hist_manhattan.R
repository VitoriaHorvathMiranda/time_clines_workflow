#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(qvalue)
library(tidyverse)
library(patchwork)


#parse arguments 
parser <- ArgumentParser(description= "Makes hist of p-values and manhattan plots")
parser$add_argument('--glm', '-glm', nargs = 5, help= 'glm_script outputs, one for each chrom')
parser$add_argument('--histPlots', '-qhist', help = 'p-value histogram output from qvalue package, .jpeg')
parser$add_argument('--hist', '-hist', help= 'p-values histogram, .jpeg')
parser$add_argument('--man', '-man', help= 'Manhattan plot')
parser$add_argument('--qPlot', '-qp', help = 'q-plot output, .jpeg')
parser$add_argument('--output', '-out', help = 'table with q-values')

xargs<- parser$parse_args()

# p_values_names <- c("/dados/time_clines/analysis/time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97_2L.tsv",
#               "/dados/time_clines/analysis/time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97_2R.tsv",
#               "/dados/time_clines/analysis/time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97_3L.tsv",
#               "/dados/time_clines/analysis/time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97_3R.tsv",
#               "/dados/time_clines/analysis/time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97_X.tsv")
xargs$glm
#p_values <- lapply(p_values_names, fread)
p_values <- lapply(xargs$glm, fread)

chrom <- c("2L","2R","3L","3R","X")

for (i in seq_along(chrom)) {
  p_values[[i]] <- p_values[[i]][, CHROM := chrom[[i]]]
}


p_values <- rbindlist(p_values)


HIST <- 
p_values |>
  ggplot() +
  geom_histogram(aes(p_value), boundary = 0, binwidth = 0.02) +
  coord_cartesian(ylim = c(0, 2.2e5)) +
  labs(x = "p-value",
       caption = "min-freq = 0.001, min-cov = 15, min-count = 5") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14))


######get FDR cutof 1997 ---------------
Pvalue <- p_values$p_value #all_ps[year == "1997"]$p_value
qobj <- qvalue(p = Pvalue)

#cutoff 10%
p_value_cutoff_0.1 <- max(qobj$pvalues[qobj$qvalues <= 0.1])

#cutoff 0.05%
p_value_cutoff_0.05 <- max(qobj$pvalues[qobj$qvalues <= 0.05])

#cutoff 1%
p_value_cutoff_0.01 <- max(qobj$pvalues[qobj$qvalues <= 0.01])

###### adds q-value column
p_values[, qvalue := qobj$qvalues]

fwrite(p_values, file = xargs$output, sep = "\t")

#### q-plots -----------------------------------------------

jpeg(filename = xargs$histPlots,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

hist(qobj)

dev.off()

jpeg(filename = xargs$qPlot,
     width = 23,
     height = 20,
     units = "cm",
     res = 1200)

plot(qobj) 

dev.off()

################ Manhattan plots ------------------------------------
####### 1997
setkey(p_values, CHROM, POS)

axisdf <-
  p_values[, ID := .I] |> 
  group_by(CHROM) |> 
  summarize(center=( max(ID) + min(ID) ) / 2 )

MANHATTAN <- 
p_values |>
  ggplot(aes(x = ID, y =-log10(p_value))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.7) +
  geom_hline(aes(yintercept = -log10(p_value_cutoff_0.1)), color = "darkgreen") +
  geom_hline(aes(yintercept = -log10(p_value_cutoff_0.01)), color = "blue") +
  geom_hline(aes(yintercept = -log10(p_value_cutoff_0.05)), color = "red") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.1) + 0.25,
                label = "FDR 10%"),
            size = 2.5, color = "darkgreen") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.01) + 0.25,
                label = "FDR 1%"),
            size = 2.5, color = "blue") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.05) + 0.25,
                label = "FDR 5%"),
            size = 2.5, color = "red") +
  scale_color_manual(values = rep(c("lightblue3", "lightcoral"), 5)) +
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center) +
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


#### save images -----------------------------------------------------

jpeg(filename = xargs$man,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
MANHATTAN

dev.off()


jpeg(filename = xargs$hist,
     width = 26,
     height = 14,
     units = "cm",
     res = 1200)
HIST

dev.off()
