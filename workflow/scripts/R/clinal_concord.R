
library(argparse)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(VennDiagram)
library(patchwork)
library(cowplot)


parser <- ArgumentParser(description= "compares clinal snps")
#inputs:
parser$add_argument('--q97', '-q97', help= 'q-values_97 file')
parser$add_argument('--q0910', '-q0910', help= 'q-values_0910 file')
#outputs:
parser$add_argument('--painel', '-painel', help= 'painel of 4 plots')
parser$add_argument('--concord', '-concord', help = 'table with concordant ratio')
parser$add_argument('--qsu', '-qsu', help = 'get summaries of qvalues in 97 for snps that are clinal in 0910 but not in 97')
parser$add_argument('--meanSdiff', '-meanSdiff', help = "plot with slopes mean diff")
parser$add_argument('--meanSdiffCoef', '-meanSdiffCoef', help = "table with weakening vs strengthening slopes lms")
parser$add_argument('--coeff', '-coeffs', help = "coefficients of slopes lm")
parser$add_argument('--incr', '-incr', help = "table with 10% highest increasing slope snps")
xargs<- parser$parse_args()

# q97 <- fread("/dados/time_clines/analysis/time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv")
# q0910 <- fread("/dados/time_clines/analysis/time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910all.tsv")

q97 <- fread(xargs$q97)
q0910 <- fread(xargs$q0910)
head(q97)
head(q0910)

setkey(q97, CHROM, POS)
setkey(q0910, CHROM, POS)

setcolorder(q97, c("CHROM", "POS"))
setcolorder(q0910, c("CHROM", "POS"))

newnames97 <- paste0(colnames(q97)[-1:-2], "_97")
newnames0910 <- paste0(colnames(q0910)[-1:-2], "_0910")

setnames(q97, colnames(q97)[-1:-2], newnames97)
setnames(q0910, colnames(q0910)[-1:-2], newnames0910)

all_q <- 
merge.data.table(q97, q0910)

# all_q[, slope_coefficients_97 := exp(slope_coefficients_97)/ 1 + exp(slope_coefficients_97)]
# all_q[, slope_coefficients_0910 := exp(slope_coefficients_0910)/ 1 + exp(slope_coefficients_0910)]


#proportion of concordant clinal SNPs ------------------------------------------

bins <- seq(0.05,0.2, by = 0.01)

clinal97 <- vector("double", length(bins))
clinal0910 <- vector("double", length(bins))
clinal0910_shifted <- vector("double", length(bins))
total <- vector("double", length(bins))

for (i in seq_along(bins)) {
  total[i] <- nrow(all_q[qvalue_97 <= bins[i] & qvalue_0910 <= bins[i]])
  clinal97[[i]] <- 
    nrow(all_q[qvalue_97 <= bins[i] & qvalue_0910 <= bins[i]])/nrow(all_q[qvalue_97 <= bins[i]])
  clinal0910[[i]] <- 
    nrow(all_q[qvalue_0910 <= bins[i] & qvalue_97 <= bins[i]])/nrow(all_q[qvalue_0910 <= bins[i]])
  # clinal0910_shifted[[i]] <- 
  #   nrow(all_q[qvalue_0910 <= bins[i] & qvalue_97 <= bins[i] + 0.05])/nrow(all_q[qvalue_0910 <= bins[i]])
}


# all_q[qvalue_0910  <= 0.1 & qvalue_97 >= 0.1,
#       .(median = median(qvalue_97),
#         mean = mean(qvalue_97))]


concor <- paste("concor", bins, sep = "_")

for (i in seq_along(bins)) {
  all_q[qvalue_97 <= bins[i] & qvalue_0910 <= bins[i],
        concor[i] := ifelse((slope_coefficients_97 > 0 & slope_coefficients_0910 > 0) |
                           (slope_coefficients_97 < 0 & slope_coefficients_0910 < 0),
                         1, 0)]
  
}

N_concord_slope <- 
  sapply(concor, function(x) 
    all_q[get(x) == 1, .N][1])

table <- 
data.table(bin = bins, clinal97 = clinal97, 
           clinal0910 = clinal0910,
           #clinal0910_shifted = clinal0910_shifted,
           total = total,
           N_concord_slope = N_concord_slope)

table[, prop_concord_slope := N_concord_slope/total]

fwrite(table,
       file = xargs$concord,
       sep = "\t")

Prop_plot <- 
melt.data.table(table,
                id.vars = c("bin", "total", "N_concord_slope", "prop_concord_slope"),
                measure.vars = c("clinal0910", "clinal97"),
                variable.name = "year", value.name = "prop_clinal") |>
  ggplot(aes(x = bin, y = prop_clinal, color = year)) +
  geom_line() + #stat = "identity", position = "dodge" 
  geom_point() +
  scale_color_brewer(palette = "Dark2", labels = c("2009/2010", "1997")) +
  #scale_color_manual(labels = c("2009/2010", "1997")) +
  labs(y = "Proportion of clinal SNPs on both times", x = "q-value",
       color = "In respect to:") + 
  theme_minimal()

#Venn diagram ------------------------------------------------------------------- 
grid.newpage() 

Venn_q01 <- draw.pairwise.venn(
    area1 = nrow(all_q[qvalue_97 <= 0.1]),    # Area of the first set
    area2 = nrow(all_q[qvalue_0910 <= 0.1]),  # Area of the second set
    cross.area = table[bin == 0.1]$total,     # Overlap area
    category = c("1997", "2009/2010"),        # Category labels
    fill = c("#1B9E77", "#D95F02"),           # Fill colors
    alpha = 0.5,                              # Transparency of the circles
    cex = 1,                                  # Font size for the numbers
    cat.cex = 2,                            # Font size for the category labels
    label.col = "white",                      # Color for the labels inside the circles
    xlim = c(-2, 2),                          # Control x-axis limits (adjust as needed)
    ylim = c(-2, 2),                          # Control y-axis limits (adjust as needed)
    cat.pos = c(0, 180),                      # Position of category labels (adjust angles)
    cat.dist = c(0.05, 0.06)                  # Distance of labels from the circles
  )

# Convert Venn diagram to a ggplot object using cowplot::ggdraw
Venn_gg <- ggdraw() + draw_grob(Venn_q01)

# Add a title
#grid.text("qvalue = 0.1", y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))

#get summaries of qvalues in 97 for snps that are clinal in 0910 but not in 97 -----
qmeans <- vector("list", length(bins))
for (i in seq_along(bins)) {
  qmeans[[i]] <- 
    summary(all_q[qvalue_0910  <= bins[i] & qvalue_97 >= bins[i]]$qvalue_97)
}

qmeans_dt <- rbindlist(lapply(qmeans, function(x) as.data.table(as.list(x))), 
                       use.names = TRUE, fill = TRUE)

qmeans_dt[, bin := bins]

fwrite(qmeans,
       file = xargs$qsu,
       sep = "\t")


# clinal on both ---------------------------------------------------------------
all_q[, diff_slope_abs :=
        abs(slope_coefficients_97 - slope_coefficients_0910)]

all_q[, diff_slope :=
        (abs(slope_coefficients_97) - abs(slope_coefficients_0910))]

#diferenças > 0 estão diminuindo; < 0 aumentando
all_q[, diff_sign := ifelse( diff_slope > 0, "-", "+")]


setorder(all_q, -diff_slope_abs)

#get top increasing slopes 

increasing <- 
  all_q[diff_sign == "+" & concor_0.1 == 1][
    , !c(paste0("concor_", bins)), with = FALSE]

#increasing[, pec_diff := ifelse(.I < .N/2, 0.05, 0.1)]

fwrite(increasing, xargs$incr, sep = "\t")

top_bins <- seq(0.05,1, by = 0.05)
diff_tables <- vector("list", length(concor))
means_test_coeffs <- vector("list", length(concor))
means_test_rsq <- vector("list", length(concor))
diff_n <- vector("list", length(concor))
diff_mean <- vector("list", length(concor))
change_ratio <- vector("double", length(top_bins)*length(concor))
for (i in seq_along(concor)) {
  for (j in seq_along(top_bins)) {
    diff_tables[[i]][[j]] <- all_q[get(concor[i]) == 1][
      1:ceiling(.N * top_bins[j])]
    
    diff_tables[[i]][[j]] <- diff_tables[[i]][[j]][
      , diff_sign := factor(diff_sign, levels = c('+', '-'))
    ]
    
    diff_tables[[i]][[j]][, qvalue := bins[i]]
    diff_tables[[i]][[j]][, top := top_bins[j]]
    
    means_test <- summary(lm(diff_slope_abs ~ diff_sign,
                               diff_tables[[i]][[j]]))
    
    means_test_coeffs[[i]][[j]] <- data.table(means_test$coefficients)
    means_test_rsq[[i]][[j]] <- data.table(rsq = means_test$r.squared)
    
    # Calculate the counts and mean based on 'diff_sign'
    diff_n[[i]][[j]] <- diff_tables[[i]][[j]][, .(total = .N), by = diff_sign]
    diff_mean[[i]][[j]] <- diff_tables[[i]][[j]][, .(mean = mean(diff_slope_abs)),
                                                 by = diff_sign]
    
    
    # Calculate the percentage of negative values
    total_neg <- diff_n[[i]][[j]][diff_sign == "-", total]
    total_posi <- diff_n[[i]][[j]][diff_sign == "+", total]
  
    # Use the correct index in 'perc_neg'
    change_ratio[(i - 1) * length(top_bins) + j] <- total_neg / total_posi
    
  }
  
} 

means_test_rsq_dt <- 
  lapply(means_test_rsq, rbindlist) |> rbindlist()

means_test_coeffs_dt <- 
  lapply(means_test_coeffs, rbindlist) |> rbindlist()

means_test_coeffs_dt[, qvalue := rep(bins, each = 40)]
means_test_coeffs_dt[, top_bins := rep(top_bins, each = 2, times = 16)]

means_test_rsq_dt[, qvalue := rep(bins, each = 20)]
means_test_rsq_dt[, top_bins := rep(top_bins, each = 1, times = 16)]

means_test_all_info <- 
merge.data.table(means_test_coeffs_dt, means_test_rsq_dt,
                 by = c("qvalue", "top_bins"))

means_test_all_info[, coeff := rep(c("Intercept", "strengthened"),
                                   times = 320)]

means_test_all_info <- 
  dcast.data.table(means_test_all_info,
                   qvalue + top_bins ~ coeff,
                   value.var = colnames(means_test_all_info)[3:7])

means_test_all_info <- 
means_test_all_info[, !c("rsq_Intercept", "Pr(>|t|)_Intercept"),
                    with = FALSE]

fwrite(means_test_all_info,
       file = xargs$meanSdiffCoef,
       sep = "\t")

means_test_all_info_plot <- 
  means_test_all_info |>
  ggplot(aes(x = qvalue, y = top_bins)) +
  geom_tile(aes(fill = Estimate_strengthened),
            stat = "identity") +
  geom_text(aes(label = round(Estimate_strengthened, digits = 3)),
            size = 2) +
  scale_fill_viridis(option = "C") +
  labs(x = "q-value", y = "Top most different slope",
       fill = "slope diff\nmean") +
  theme_minimal()


# diff_mean <- 
# lapply(diff_mean, rbindlist) |> rbindlist()
# 
# diff_mean[, qvalue := rep(bins, each = 40)]
# diff_mean[, top_bins := rep(top_bins, each = 2, times = 16)]
# 
# diff_mean <- 
# dcast.data.table(diff_mean, qvalue + top_bins ~ diff_sign,
#                  value.var = "mean")
# 
# diff_mean[, mean_diff := `+` - `-`] # positivos -> slope diminuiu
#                                     # negativ0s -> slope aumentou

#se a diff for +, média dos que diminuem é maior
#se a diff for -, a média dos que aumentam é maior

# diff_mean_plot <- 
# diff_mean |>
#   ggplot(aes(x = qvalue, y = top_bins)) +
#   geom_tile(aes(fill = mean_diff),
#             stat = "identity") +
#   geom_text(aes(label = round(mean_diff, digits = 3)),
#             size = 2) +
#   scale_fill_viridis(option = "C") +
#   labs(x = "q-value", y = "Top most different slope",
#        fill = "slope diff\nmean") +
#   theme_minimal()


jpeg(xargs$meanSdiff,
     width = 25,
     height = 23,
     units = 'cm',
     res = 1200)

means_test_all_info_plot

dev.off()

change_dir <- 
data.table(q_value = rep(bins, each= 20),
           top_diff = rep(top_bins, times = 16),
           change_ratio = change_ratio)

change_dir[, change_ratio := round(change_ratio, 2)]

simetry_plot <- 
change_dir |>
  ggplot(aes(x=q_value, y = top_diff)) +
  geom_tile(aes(fill = change_ratio),
             stat = "identity") +
  geom_text(aes(label = round(change_ratio, digits = 3)),
            size = 2) +
  scale_fill_viridis(option = "C") +
  labs(x = "q-value", y = "Top most different slope",
       fill = "Sign Ratio") +
  theme_minimal()


Slopes_lm_plot <- 
all_q[concor_0.1 == 1] |>
  ggplot(aes(x = slope_coefficients_97,
             y = slope_coefficients_0910)) +
  geom_point(aes(color = diff_slope_abs)) +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = 0.15, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -0.15, linetype = "dashed") +
  geom_smooth(method = 'lm') +
  labs(x = "Slope Coefficient in 1997",
       y = "Slope Coefficient in 2009/2010",
       color = "Abs. \nSlope Diff.") +
  scale_color_viridis(option = "A") +
  theme_minimal()


coeffs <- vector("list", length(concor))
for (i in seq_along(concor)) {
  model <- 
    lm(slope_coefficients_0910 ~ + slope_coefficients_97,
       all_q[get(concor[i]) == 1])
  coeffs[[i]] <- summary(model)$coefficients |> data.table()
}

coeffs <- 
rbindlist(coeffs)

coeffs[, qvalue_bin := rep(bins, each =2)]
coeffs[, coef_type := rep(c("intercept", "slope"), times = length(bins))]
setorder(coeffs, `Pr(>|t|)`)
coeffs[coef_type == "slope", n := .I]
coeffs[coef_type == "intercept", n := .I]
coeffs[,adj_critical_alpha := 0.05/(nrow(coeffs)+1 - n)]
setorder(coeffs, qvalue_bin)

fwrite(coeffs,
       file = xargs$coeff,
       sep = "\t")

# coeffs_plot <- 
# coeffs[coef_type == "slope"] |>
#   ggplot(aes(x=qvalue_bin, y = Estimate)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = Estimate - `Std. Error`,
#                     ymax =  Estimate + `Std. Error`)) +
#   labs(y = "slope coefficient") +
#   theme_minimal()

#main plot ---------------------------------------------------------------------

pdf(xargs$painel,
    width = 40/4,
    height = 29/4)

(Venn_gg + Prop_plot) / (Slopes_lm_plot + simetry_plot) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "plain")) 

dev.off()

# pdf(file = "concord_fig_test.pdf",
#     width = 40/4,
#     height = 29/4,)
# 
# (Venn_gg + Prop_plot + plot_layout(widths = c(1, 1.5))) / 
#   (Slopes_lm_plot + simetry_plot + plot_layout(widths = c(1, 1.5))) +
#   plot_layout(heights = c(1, 1.5)) +
#   plot_annotation(tag_levels = 'A') &
#   theme(plot.tag = element_text(face = "plain")) 
# 
# dev.off()
