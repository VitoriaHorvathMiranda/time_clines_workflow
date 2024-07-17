#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(data.table)
library(terra)

#parse arguments 
parser <- ArgumentParser(description= "Makes table with world_clim summary stats of tmax, tmin and prec from 6 moths pre collection month")
parser$add_argument('--meta', '-meta', help = "metadata table")
parser$add_argument('--year', '-year', help = "year")
parser$add_argument('--pathRaster', '-path', help= 'directory of raster files')
parser$add_argument('--output', '-o', help = "output")
xargs<- parser$parse_args()

meta <- fread(xargs$meta)
year <- xargs$year
raster_path <- xargs$pathRaster

meta <- meta[collection_year == year]
patt <- paste0("wc2.1_2.5m_[a-z]{4}_", year, "-[0-9]{2}.tif")
# meta <- fread('/dados/time_clines/data/meta/seq_metadata.tsv')
# meta <- meta[collection_year == "2009"]
meta[, ID := .I]

file_names <- 
list.files(path = raster_path,
           pattern = patt,
           full.names = TRUE)

rasters <- lapply(file_names, rast)

coords <- data.frame(long = meta$longitude, lat = meta$latitude)

worldclim <- lapply(rasters, function(x) terra::extract(x, coords)) |>
  lapply(as.data.table)

fill_variables <- 
function(x) {
  variable <- colnames(x)[2] |> str_split_i(pattern = "_", i = 3)
  month <- colnames(x)[2] |> str_split_i(pattern = "_", i = 4)
  x[, "variable" := variable]
  x[, "month" := month]
  
}

lapply(worldclim, fill_variables)

lapply(worldclim, function(x) setnames(x, colnames(x)[2], "value"))

all_var <- rbindlist(worldclim)
all_var[, month := as.integer(month)]

all_var <- merge.data.table(all_var, meta[, .(ID, population, collection_month)])

all_var[, six_month_pre := collection_month-6]
all_var <- all_var[six_month_pre <= month & collection_month >= month]

summary_var <- all_var[, .(min = min(value), max = max(value), mean = mean(value)),
        by = c("population", "variable")] |>
  dcast.data.table(formula = population ~ variable,
                   value.var = c("min", "max", "mean"))

summary_var <- summary_var[, !c("min_tmax", "max_tmin")]

fwrite(summary_var, file = xargs$output, sep = "\t")

