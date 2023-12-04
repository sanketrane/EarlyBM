rm(list = ls()); gc();

############################################################
### Preamble
############################################################

## loading libraries
library(tidyverse)
library(readxl)
library(deSolve)

#### plotting style
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks = as.numeric(1:10 %o% 10^(-4:8))


############################################################
############################################################

## import data
## BrdU pulse chase in developing B cells
scale_fun <- function(x) {x*1e6}

counts_df <- read_excel("BrdU_data.xlsx", sheet =1) %>%
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

ggplot(counts_df, aes(x= Time_h, y=cell_counts, col = Genotype)) +
  geom_point()+ xlim(0, 35)+
  scale_y_log10() + 
  #geom_hline(yintercept = 510111, col="darkblue")+
  #geom_hline(yintercept = 131619, col="darkred")+
  labs(x='Days post trabsfer', y=NULL)+
  facet_grid(BrdU_status ~ factor(pop_id, levels = c("Pro_B", "large_pre_B", "small_pre_B",
                                                  "BrdU_Pro_B", "BrdU_large_pre_B", "BrdU_small_pre_B"))) +
  myTheme 



fracs_df <- read_excel("datafiles/BrdU_data.xlsx", sheet =2) %>%
  mutate(across("Sample", str_replace_all, c("Stain 2 BM development" = "", " " = "_")),
         across("Experiment_Date", str_replace, "2014-", "")) %>%
  unite(sample_id, Experiment_Date, Sample) %>%
  select(sample_id, Time_h, Genotype, 
         BrdU_Pro_B, BrdU_large_pre_B, BrdU_small_pre_B) %>%
  gather(-c(sample_id, Time_h, Genotype), key = 'subpop', value = "prop_brdu")


ggplot(fracs_df, aes(x= Time_h, y=prop_brdu, col = Genotype)) +
  geom_point()+ xlim(0, 35) + ylim(0, 100)+ #scale_y_continuous(trans = "log10", limits = c(1, 100))+
  labs(x="Hours post BrdU injection", y="% BrdU+ cells") +
  facet_grid(Genotype ~ factor(subpop, levels = c("BrdU_Pro_B", "BrdU_large_pre_B", "BrdU_small_pre_B"))) +
  myTheme + guides(col="none")

