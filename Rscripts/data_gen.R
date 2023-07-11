### data wrangling for Marginal Zone compartmnet
rm(list = ls()); gc();
## loading libraries
library(tidyverse)
library(readxl)

## import data
## BrdU pulse chase in developing B cells
scale_fun <- function(x) {x*1e6}

counts_df <- read_excel("datafiles/BrdU_data.xlsx", sheet =1) %>%
  mutate(across("Sample", str_replace_all, c("Stain 2 BM development" = "", " " = "_")),
         across("Experiment_Date", str_replace, "2014-", "")) %>%
  unite(sample_id, Experiment_Date, Sample) %>%
  mutate_at(vars(contains("_B")), scale_fun) %>%
  select(-Total_cells_BM) %>%
  select(sample_id, Time_h, Genotype, Pro_B, large_pre_B, small_pre_B, 
         BrdU_Pro_B, BrdU_large_pre_B, BrdU_small_pre_B) %>%
  gather(-c(sample_id, Time_h, Genotype), key = 'subpop', value = "cell_counts") %>%
  mutate(pop_id = ifelse(grepl("Pro_B", subpop), "Pro_B",
                         ifelse(grepl("large_pre_B", subpop), "large_pre_B", 
                                "small_pre_B")),
         BrdU_status = ifelse(grepl("BrdU", subpop), "BrdU_pos", 
                              "Total"))

fracs_df <- read_excel("datafiles/BrdU_data.xlsx", sheet =2) %>%
  mutate(across("Sample", str_replace_all, c("Stain 2 BM development" = "", " " = "_")),
         across("Experiment_Date", str_replace, "2014-", "")) %>%
  unite(sample_id, Experiment_Date, Sample) %>%
  select(sample_id, Time_h, Genotype, 
         BrdU_Pro_B, BrdU_large_pre_B, BrdU_small_pre_B)

fracs_wt <- fracs_df %>%
  filter(Genotype == "WT") 
fracs_dko <- fracs_df %>%
  filter(Genotype == "dKO") 

## Unique time points with indices to map
solve_time <- c(4, 18, 30)
time_index1 <- purrr::map_dbl(fracs_wt$Time_h, function(x) which(x == solve_time))     # keeping track of index of time point in relation to solve_time
time_index2 <- purrr::map_dbl(fracs_dko$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs1 <- nrow(fracs_wt)
numObs2 <- nrow(fracs_dko)
n_shards <- length(solve_time)
largePreB_wt <- fracs_wt$BrdU_large_pre_B/100
largePreB_dko <- fracs_dko$BrdU_large_pre_B/100
smallPreB_wt <- fracs_wt$BrdU_small_pre_B/100
smallPreB_dko <- fracs_dko$BrdU_small_pre_B/100
# time sequence for predictions specific to age bins within the data
ts_pred <- 10^seq(log10(0.1), log10(30), length.out = 300)
numPred <- length(ts_pred)


stan_rdump(c("numObs1",  "numObs2", "n_shards", "solve_time", "time_index1", "time_index2",
             "largePreB_wt",  "smallPreB_wt", "largePreB_dko", "smallPreB_dko",
             "ts_pred", "numPred"),
             file = file.path('datafiles', paste0('brdu_stanfit',".Rdump")))








