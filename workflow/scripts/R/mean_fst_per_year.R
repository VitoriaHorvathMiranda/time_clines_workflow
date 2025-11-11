#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(geosphere)

#parse arguments
parser <- ArgumentParser(description= "boxplots FST per year and lm fst ~year")
parser$add_argument('--FST', '-fst', help= 'FST pair list')
parser$add_argument('--outputLM', '-lm', help= 'output with lm summarys')
parser$add_argument('--outputBX', '-box', help= 'box plot pest per year, jpeg')
parser$add_argument('--distanceLM', '-distLM', help = 'distance lm summary')
parser$add_argument('--distFST', '-dist', help = 'scatter plot of distance per fst per collection year, jpeg')
parser$add_argument('--meta', '-meta')
xargs<- parser$parse_args()

# meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv",
#               select = c("population", "latitude", "longitude"))
# FST_auto <- fread("/dados/time_clines/analysis/fst/genome/Genome_FST_autosome_fst-list.csv")
FST_auto <- fread(xargs$FST)
meta <- fread(xargs$meta)

setnames(FST_auto, "FST", "fst")

FST_auto[, year_1 := fcase(first %like% "10" | first %like% "09", "2009/2010",
                           first %like% "97", "1997",
                           first %like% "17", "2017",
                           first %like% "22" | first %like% "23", "2022/2023")]

FST_auto[, year_2 := fcase(second %like% "10" | second %like% "09", "2009/2010",
                           second %like% "97", "1997",
                           second %like% "17", "2017",
                           second %like% "22" | second %like% "23", "2022/2023")]


FST_auto[, matching_year := fcase(year_1 == "2009/2010" & year_2 == "2009/2010", "2009/2010",
                                  year_1 == "1997" & year_2 == "1997", "1997",
                                  year_1 == "2017" & year_2 == "2017", "2017",
                                  year_1 == "2022/2023" & year_2 == "2022/2023", "2022/2023",
                                  default = "no_match")]

FST_with17 <- FST_auto[matching_year != "no_match"]
FST <- FST_auto[matching_year != "no_match" & year_1 != "2017"]

# test differences -------------------------------------------------------
#não usei esse modelo
FST_with17[, collection_year := case_when(matching_year == "2009/2010" ~ 2009.5,
                                        matching_year == "2022/2023" ~ 2022.5,
                                        .default = as.double(matching_year))]

lm <- lm(fst ~ collection_year,
         data = FST_with17)
mary <- summary(lm)

# table <- 
# mary$coefficients |> 
#   as.data.frame() |>
#   rownames_to_column() 

#remove CMD97B just to be sure
lm_noCMD97B <- lm(fst ~ collection_year,
   data = FST_with17[first != "CMD97B" & second != "CMD97B"])

sink(xargs$outputLM)
print(summary(lm))
print(summary(lm_noCMD97B))
sink()


# mary_noCMD97B <- 
# summary(lm_noCMD97B)
# 
# table_noCMD97B <- 
# mary_noCMD97B$coefficients |> 
#   as.data.frame() |>
#   rownames_to_column() 


#saves -------
# setnames(table_noCMD97B, colnames(table_noCMD97B)[-1],
#          paste0(colnames(table_noCMD97B)[-1], "_noCMD97B"))
# 
# save_table <- 
# merge(table, table_noCMD97B, by = ("rowname"))
# 
# fwrite(save_table, file = xargs$outputLM, sep = "\t")

# plot -------------------
MEAN_FST_YEAR <- 
FST[matching_year != "no_match"] |>
  ggplot() +
  geom_boxplot(aes(x = matching_year, y = fst, fill = matching_year)) +
  labs(x = "Year") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#E7298A")) +
  theme_light() +
  theme(legend.position = "none")

jpeg(filename = xargs$outputBX,
     width = 23,
     height = 19,
     units = "cm",
     res = 1200)

MEAN_FST_YEAR

dev.off()


#FST per dist -----------------------------------------------------

#compute dist--
FST <- 
  merge.data.table(FST, meta, by.x = "first", by.y = "population")

FST <- 
  merge.data.table(FST, meta, by.x = "second", by.y = "population")


calculate_distance <- function(lon1, lat1, lon2, lat2) {
  point1 <- c(as.numeric(lon1), as.numeric(lat1))
  point2 <- c(as.numeric(lon2), as.numeric(lat2))
  distHaversine(point1, point2) / 1000  # Convert to kilometers
}


FST[, distance := mapply(calculate_distance,
                         lon1 = longitude.x,
                         lat1 = latitude.x,
                         lon2 = longitude.y,
                         lat2 = latitude.y)]
FST_DIST <- 
FST[matching_year != "2022/2023"] |> 
  ggplot(aes(x = distance, y = fst, color = matching_year)) +
  geom_point(aes(), size = 5) +
  geom_smooth(method = 'lm') +
  labs(x = "Distance in km", y = expression("F"[ST]),
       color = "Collection Year") +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#E7298A")) +
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

jpeg(filename = xargs$distFST, #"fst_year_dist_autossome.jpeg", 
     width = 27,
     height = 15,
     units = "cm",
     res = 1200)

FST_DIST

dev.off()


distance_lm <- lm(fst~ distance + matching_year,
                  FST[matching_year != "2022/2023"])

distance_lm_int <- lm(fst~ distance * matching_year,
                  FST[matching_year != "2022/2023"])

distance_lm_noCMD97B <- lm(fst~ distance + matching_year,
                           FST[first != "CMD97B" & second != "CMD97B" &
                                 matching_year != "2022/2023"])

distance_lm_int_noCMD97B <- lm(fst~ distance * matching_year,
                               FST[first != "CMD97B" & second != "CMD97B" &
                                     matching_year != "2022/2023"])


sink(xargs$distanceLM)
print(summary(distance_lm))
print(summary(distance_lm_int))
print(anova(distance_lm_int, distance_lm))
sink()

sink("lm_fst_dist_year_noCMD97B.txt")
summary(distance_lm_noCMD97B)
summary(distance_lm_int_noCMD97B)
anova(distance_lm_noCMD97B, distance_lm_int_noCMD97B)
sink()



