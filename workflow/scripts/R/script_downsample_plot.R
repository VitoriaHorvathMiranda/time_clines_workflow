library(tidyverse)
library(data.table)
library(patchwork)
library(argparse)

#parse arguments
parser <- ArgumentParser(description= "Makes depth plot per coverage and per position")
parser$add_argument('--meta', '-m', help= 'metadados table')
parser$add_argument('--PreCovPath', '-precp', help= 'Pre downsample samtools coverage path')
parser$add_argument('--PosCovPath', '-poscp', help= 'Pos downsample samtools coverage path')
parser$add_argument('--DepthperPOS', '-dpos', help= 'path to folder with script_depth_per_pos output')
parser$add_argument('--output', '-o', help= 'Single Jpeg with all plots')
parser$add_argument('--boxplot', '-box', help = 'boxplot of mean depth per chrom per year')
xargs<- parser$parse_args()



#Coverage vs Mean Depth Plots ---------------------------------------------------
## Read data -------------------------------------------------------------------
#meta <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")
meta <- fread(xargs$meta, fill = TRUE)

meta[, population := ifelse(population == "HFL97_new", "ESC97_n", population)]

files_cov_pre_down <- list.files(path = "/dados/time_clines/data/seqs/processed/qltctrl/pre_downsample/samtools_coverage",
                                 pattern = "coverage.tsv", full.names = TRUE)

files_cov_pos_down1 <- list.files(path = "/dados/time_clines/data/seqs/processed/qltctrl/pos_downsample/samtools_coverage",
                                 pattern = "(1|2)_total_coverage.tsv", full.names = TRUE)

files_cov_pos_down2 <- list.files(path = "/dados/time_clines/data/seqs/processed/qltctrl/pos_downsample/samtools_coverage",
                                 pattern = "_all_chrom_total_coverage.tsv", full.names = TRUE)

male_ids <- c("dlSC10", "dlGA10", "dlFL10", "HFL97downto60mi")

files_cov_pos_down_all <- files_cov_pos_down1[!str_detect(files_cov_pos_down1, paste(male_ids, collapse = "|"))]

files_cov_pos_down_all <- c(files_cov_pos_down_all, files_cov_pos_down2)



# files_cov_pre_down <- list.files(path = xargs$PreCovPath,
#                                  pattern = "coverage.tsv", full.names = TRUE)
# files_cov_pos_down <- list.files(path = xargs$PosCovPath,
#                                  pattern = "coverage.tsv",  full.names = TRUE)


#gets the samples ids because files are out of order
sample_ids_pre <- lapply(files_cov_pre_down, str_split_i, pattern = "/", i = 10) |>
  lapply(str_split_i, pattern = "_", i = 1)

sample_ids_pos_down_all <- lapply(files_cov_pos_down_all, str_split_i, pattern = "/", i = 10) |>
  lapply(str_split_i, pattern = "_", i = 1)

cov_pre_down <- lapply(files_cov_pre_down, fread) # reads all tables and stores them in a list
lapply(cov_pre_down, setnames, "#rname", "chrom") # change the first collunm name of every list
chrom_pre_down <- lapply(cov_pre_down, # filters only relevant chrom arms
                         filter, chrom %in% c("2L", "2R", "3L", "3R", "X"))

cov_pos_down_all <- lapply(files_cov_pos_down_all, fread)
lapply(cov_pos_down_all, setnames, "#rname", "chrom")
chrom_pos_down1 <- lapply(cov_pos_down_all,
                         filter, chrom %in% c("2L", "2R", "3L", "3R", "X"))


print("coverage files:")
files_cov_pre_down
files_cov_pos_down
print("sample_ids:")
sample_ids
print("meta file:")
meta
## adds pops names and merge files ----------------------------------------------

get_names <- 
function(data, name) {
  output <- vector("list", length(data))
  for (i in seq_along(name)) {
    output[[i]] <- data[[i]][,seq_label := name[[i]]]
  }
  return(rbindlist(output))
}

all_cov_pre_down <- get_names(chrom_pre_down, sample_ids_pre)
all_cov_pos_down1 <- get_names(chrom_pos_down1, sample_ids_pos_down_all)


## plots -----------------------------------------------------------------------

boxplot_data <- 
function(all_cov, meta){
test <- 
all_cov |>
  mutate(seq_label = if_else(seq_label == "09", "9", seq_label)) |>
  left_join(meta, by = "seq_label") |>
  mutate(year = case_when(population %like% "97" ~ "1997",
                          population %like% "09" | population %like% "10" ~ "2009/2010",
                          population %like% "17" ~ "2017",
                          population %like% "22" | population %like% "23" ~ "2022/2023")) 
return(test)

}


box_table_pre <- boxplot_data(all_cov_pre_down, meta)
box_table_pos1 <- boxplot_data(all_cov_pos_down1, meta)

BOX_PRE <- 
box_table_pre |>
  ggplot(aes(y = meandepth)) +
  geom_boxplot(aes(color = year)) +
  facet_wrap(vars(chrom), nrow = 1) +
  labs(y = "Mean Read Depth", 
       title = "Before downsample") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 

BOX_POS <- 
box_table_pos1 |>
  ggplot(aes(y = meandepth)) +
  geom_boxplot(aes(color = year)) +
  facet_wrap(vars(chrom), nrow = 1) +
  labs(y = "Mean Read Depth", 
       title = "After downsample")  +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank())


jpeg(filename = xargs$boxplot,
     width = 25,
     height = 30,
     units = "cm",
     res = 1200)

BOX_PRE / BOX_POS

dev.off()



summarise_depths <- 
  function(data){
    plot_table <- data |>
      group_by(population, year, n_females, fly_sex) %>%
      summarise(total_pos = sum(endpos),
                total_covbases = sum(covbases),
                total_coverage = sum(covbases)/sum(endpos),
                total_meandepth = weighted.mean(meandepth, endpos)) |>
      ungroup()
      return(plot_table)
  }



plot_table_pre <- summarise_depths(box_table_pre)
plot_table_pos1 <- summarise_depths(box_table_pos1)

A <- plot_table_pre %>%
  ggplot(aes(x = total_meandepth, y = total_coverage)) +
  geom_point(aes(color = year), size = 2, show.legend = FALSE) +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.003) +
  labs(x = "Mean Read Depth", y = "Coverage",
       title = "Before downsample") +
  #theme(legend.position = "none") +
  theme_light()

B <- plot_table_pos1 %>%
  ggplot(aes(x = total_meandepth, y = total_coverage)) +
  geom_point(aes(color = year), size = 2) +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.003) +
  labs(x = "Mean Read Depth", y = "Coverage",
       title = "After downsample") +
  theme_light()


# Depths per Pos PLot ----------------------------------------------------------
## read window depth -----------------------------------------------------------

depth_per_win_paths <- 
  list.files(path = xargs$DepthperPOS,
           pattern = '_80kb_window_(2L|2R|3R|3L|X).tsv',  full.names = TRUE)

depths_per_window <- lapply(depth_per_win_paths, fread)

depths_per_window <- rbindlist(depths_per_window)
depths_per_window <- depths_per_window[CHROM %in% c("2L", "2R", "3L", "3R", "X")]
#depths_per_window[, window := as.factor(window)]

depths_per_window <- 
depths_per_window[meta[, .(seq_label, population)], on = "seq_label"]
depths_per_window <- depths_per_window[CHROM %in% c("2L", "2R", "3L", "3R", "X")]


#extracts smallest position of each window
depths_per_window[,min_pos := tstrsplit(window, ",", keep = 1)][
  , min_pos := sub(pattern = "\\(",
                   replacement = "",
                   x = min_pos)]


#computes mid position of each window
depths_per_window[, min_pos := as.double(min_pos)][
  , mid_pos := (min_pos+40000)
]


##plots
C <- depths_per_window %>%
  ggplot() +
  geom_line(aes(x = mid_pos, y = mean_depth_pre,
                color = population, group = population),
            #show.legend = FALSE,
            linewidth = 0.3) +
  facet_wrap(vars(CHROM), scales = "free_x") +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5),
        legend.position = "none") +
  labs(y = "Mean Depth", x = "Position",
       title = "Mean Depth Before Downsample")


D <- depths_per_window %>%
  ggplot() +
  geom_line(aes(x = mid_pos, y = mean_depth_pos,
                color = population, group = population),
            #show.legend = FALSE,
            linewidth = 0.3) +
  facet_wrap(vars(CHROM), scales = "free_x") +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 5)) +
  labs(y = "Mean Depth", x = "Position",
       title = "Mean Depth After Downsample")


## Join Plots
plot_all <- (A + B) / (C + D) +
  plot_annotation(tag_levels = "A")# + 



xargs$output

jpeg(filename = xargs$output,
     width = 25,
     height = 15,
     units = "cm",
     res = 1200)
plot_all

dev.off()


