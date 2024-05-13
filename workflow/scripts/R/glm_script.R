#!/usr/bin/env Rscript

library(argparse)
library(purrr)
library(data.table)


#parse arguments
parser <- ArgumentParser(description= "computes a glm for each snp")
parser$add_argument('--timePops', '-tPops',
                    help = 'outputs from separate_time_pops_script.R')
parser$add_argument('--output', '-o',
                    help = 'table with inclination coefficient and 
                    p-values for each SNP, .tsv')

xargs<- parser$parse_args()


#reads data 
joined_pops <- fread(file = xargs$timePops, header = FALSE)
setnames(joined_pops, colnames(joined_pops), c("population", "latitude", "n_chrom", 
"CHROM", "POS", "REF", "ALT", "depth", "freq","NE"))

print("data read")

#finds SNPs that aren't present in that period or that are present in only
#one  or two population
joined_pops[, freq_status := 1*(freq!= 0)]
no_0freq <- joined_pops[, .(freq_status_sum = sum(freq_status)),
                        by = "POS"][freq_status_sum <= 2,]

#filter those SNPs out
joined_pops <- joined_pops[!no_0freq, on = "POS"]

#data.table function for nesting:
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}


#nest data by snp
nested_snps <- joined_pops[,group_nest_dt(.SD, POS)]

print("table nested")
#runs glm for each snp
nested_snps[, models := purrr::map(data, ~ glm(freq~latitude, 
                                               weights = NE,
                                               data = .x,
                                               family = binomial()))]

print("models done")
#drop data column 
nested_snps <- nested_snps[, !c("data")]

#gets lat coefficient
nested_snps[, slope_coefficients := purrr::map_dbl(models, ~coef(.x) |> pluck("latitude")) ]
nested_snps[, intercept_coefficients := purrr::map_dbl(models, ~coef(.x) |> pluck("(Intercept)")) ]


print("coefficients collected")
#gets p-value
nested_snps[, p_value := purrr::map_dbl(models, 
                                         ~summary(.x) %>% 
                                           pluck("coefficients") %>% pluck(8))]
print("p-values collected")
#drops models column (it is too heavy)
nested_snps <- nested_snps[, !c("models")]
print("models droped")

#saves
output_path <- xargs$output
fwrite(nested_snps, output_path, sep = "\t")

