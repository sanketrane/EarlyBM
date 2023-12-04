## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(tidyverse)
####################################################################################

#### plotting style
myTheme <- theme(text = element_text(size = 15), axis.text = element_text(size = 15),
                 axis.title =  element_text(size = 15, face = "bold"),
                 plot.title = element_text(size=15,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())


####### plotting
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

log10minorbreaks=as.numeric(1:10 %o% 10^(3:8))

## Setting all the directories for opeartions
projectDir <- getwd()
dataDir <- file.path(projectDir, "datafiles")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')


## model specific details that needs to be change for every run
modelName <- "M2_base"

stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_5",".csv")))

fit1 <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars_M1 <- which(fit1@model_pars %in% "sigma1") -1      
params_M1 <- fit1@model_pars[1:num_pars_M1]


# compiling multiple stan objects together that ran on different nodes

modelName <- "M2_base_diffMu"

stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_5",".csv")))

fit2 <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5))


# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars_M2 <- which(fit2@model_pars %in% "sigma1") -1      
params_M2 <- fit2@model_pars[1:num_pars_M2]



### plots

expit_func <- function(x){exp(x)/(1+exp(x))}

Param_mat <- as.data.frame(fit2, pars = params_M2) %>%
  mutate_at(vars(matches("Log")), exp) %>%
  mutate(rho_dko_Log = rho_Log * expit_func(kappa),
         lambda_dko_Log = lambda_Log,
         convert_Log = 1/mu_Log,
         convert_dko_Log = 1/mu_dko_Log,
         lifespan_Log = 1/lambda_Log,
         lifespan_dko_Log = 1/lambda_Log,
         intdiv_Log = 1/rho_Log,
         intdiv_dko_Log = 1/rho_dko_Log) %>%
  select(contains("Log")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            estimate = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>% 
  mutate(Genotype = ifelse(grepl("dko", key), "dKO", "WT"),
         param = ifelse(grepl("rho", key), "rho",
                        ifelse(grepl("lambda", key), "lambda",
                               ifelse(grepl("mu", key), "mu", 
                                      ifelse(grepl("intdiv", key), "intdiv", 
                                             ifelse(grepl("convert", key), "convert", "lifespan"))))))

fac_labels <- c(`rho`= paste0("Rate of division (/h)"), 
                `lambda`="Rate of loss (/h)",
                `mu`= paste0("Rate of conversion \n from transient to resting (/h)"),
                `lifespan`= paste0("Mean residence time \n in Pre B (h)"),
                `convert`= paste0("Mean residence time \n in transient state (h)"),
                `intdiv`= paste0("Mean time between \n two division events (h)"))

blank_data <- data.frame(param = Param_mat$param,
                         Genotype = Param_mat$Genotype,
                         estimate = c(0.01, 0.001, 0.001, 0.001,0.01, 0.001, 3, 0.3, 3, 3, 3, 3,
                                      0.3, 0.1, 0.2, 1.5, 0.3, 0.1, 400, 50, 300, 300, 30, 100))

fac_order <- c("rho",  "mu", "lambda", "intdiv", "convert", "lifespan")


ggplot(Param_mat, aes(y=estimate, x=(Genotype), col=Genotype))+
  labs(y=NULL) +
  geom_errorbar(aes(y=estimate, ymin=lb, ymax=ub, x=Genotype),
                width=0.2, linetype=1,  position=position_dodge(0.4)) +
  geom_blank(data = blank_data)+
  geom_point(position=position_dodge(width=0.4), stat = "identity", size=4) + 
  facet_wrap(~ factor(param, levels = fac_order), scales = "free", ncol=3,
             labeller = as_labeller(fac_labels)) + 
  #expand_limits(y = 0.01)  +
  scale_y_log10()+
  myTheme + theme(axis.text.x=element_blank(),
                  axis.title.x = element_blank())

ggsave(file.path("output_fit", paste0("param_plot", modelName, ".pdf")), last_plot(), device = "pdf", width = 10, height = 8, units = "in")














