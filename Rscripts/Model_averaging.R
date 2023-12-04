rm(list = ls())  
gc()    
library(loo)
library(tidyverse)
library(gridExtra)
library(formattable)

## directories for saving outputs
OutputDir <- file.path('output_fit')

customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customOrange = "#FF8C00"
customRed = "#FF6347"

## function for exporting a data table of delta loo ic values and akaike weights
#takes 2 separate lists of the name of the models and loo-ic values
model_compare <- function(looiclist){
  # delta loo-ic
  deltaloo_list <- looiclist - min(looiclist)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('deltaloo' = round(deltaloo_list, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  colnames(export_table)[1:2] <-  c(paste0('\u0394', 'LooIC'),  paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}

LooDir <- file.path("loo_fit")

### reading the loo objects for each model
M2_base = paste0('loosave_M2_base.rds')
M2_base_NoMu = paste0('loosave_M2_base_NoMu.rds')
M2_base_diffLoss = paste0('loosave_M2_base_diffLoss.rds')
M2_base_diffMu = paste0('loosave_M2_base_diffMu.rds')
M2_diffLoss = paste0('loosave_M2_base_diffLoss.rds')
M2_diffMu= paste0('loosave_M2_diffMu.rds')

model_list <- list("M2_base" = readRDS(file.path(LooDir, M2_base)), 
                   "M2_base_NoMu" = readRDS(file.path(LooDir, M2_base_NoMu)), 
                   "M2_base_diffLoss" = readRDS(file.path(LooDir, M2_base_diffLoss)), 
                   "M2_base_diffMu" = readRDS(file.path(LooDir, M2_base_diffMu)), 
                   "M2_diffLoss" = readRDS(file.path(LooDir, M2_diffLoss)))
compare_mods <- loo_compare(model_list)
print(compare_mods, simplify = F)

compare_mods
loo_model_weights(model_list, method = 'pseudobma') * 100


