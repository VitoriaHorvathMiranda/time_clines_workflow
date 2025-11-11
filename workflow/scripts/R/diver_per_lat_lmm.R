library(data.table)
library(tidyverse)
#library(lme4)
library(lmerTest)

Global_pi <- fread("/dados/time_clines/analysis/pop_stats/Global_stat_pi_per_pop.tsv")
Global_th <- fread("/dados/time_clines/analysis/pop_stats/Global_stat_th_per_pop.tsv")

lm_pi_lat <- lmer(data = Global_pi, global_stat ~ latitude + (1 | CHROM))
lm_th_lat <- lmer(data = Global_th, global_stat ~ latitude + (1 | CHROM))

lm_pi_y <- lmer(data = Global_pi, global_stat ~ collection_year + (1 | CHROM))
lm_th_y <- lmer(data = Global_th, global_stat ~ collection_year + (1 | CHROM))


lm_pi_lat_noCMD97B <- lmer(data = Global_pi[population != "CMD97B"], global_stat ~ latitude + (1 | CHROM))
lm_th_lat_noCMD97B <- lmer(data = Global_th[population != "CMD97B"], global_stat ~ latitude + (1 | CHROM))

lm_pi_y_noCMD97B <- lmer(data = Global_pi[population != "CMD97B"], global_stat ~ collection_year + (1 | CHROM))
lm_th_y_noCMD97B <- lmer(data = Global_th[population != "CMD97B"], global_stat ~ collection_year + (1 | CHROM))


sink("diver_per_lat_lmm.txt")

summary(lm_pi_lat)
print('#--------------------------------------------------------')
summary(lm_pi_lat_noCMD97B)
print('#--------------------------------------------------------')
summary(lm_pi_y)
print('#--------------------------------------------------------')
summary(lm_pi_y_noCMD97B)
print('#--------------------------------------------------------')
summary(lm_th_lat)
print('#--------------------------------------------------------')
summary(lm_th_lat_noCMD97B)
print('#--------------------------------------------------------')
summary(lm_th_y)
print('#--------------------------------------------------------')
summary(lm_th_y_noCMD97B)

sink()

Global_pi |>
  ggplot(aes(x = collection_year, y = global_stat)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(vars(CHROM))
