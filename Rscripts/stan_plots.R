## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
library(readxl)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "M2_diffMu"

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "Rscripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "datafiles")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv/rsquare')
LooDir <- file.path('loo_fit') 

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_5",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma4")      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burnin samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

### data munging
brdu_df <- read_excel("datafiles/BrdU_data.xlsx", sheet =2) %>%
  mutate(across("Sample", str_replace_all, c("Stain 2 BM development" = "", " " = "_")),
         across("Experiment_Date", str_replace, "2014-", "")) %>%
  unite(sample_id, Experiment_Date, Sample) %>%
  select(sample_id, Time_h, Genotype, 
         BrdU_Pro_B, BrdU_large_pre_B, BrdU_small_pre_B) %>%
  mutate(prop_large_PreB = BrdU_large_pre_B/100,
         prop_small_PreB = BrdU_small_pre_B/100)%>%
  filter(sample_id != "03-19__BM_13.fcs")

fracs_wt <- brdu_df %>%
  filter(Genotype == "WT") 

fracs_dko <- brdu_df %>%
  filter(Genotype == "dKO") 


# ################################################################################################
# calculating PSIS-L00-CV for the fit
largePrBwt_loglik <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
largePrBdko_loglik <- extract_log_lik(fit, parameter_name = "log_lik3", merge_chains = TRUE)
smallPrBwt_loglik <- extract_log_lik(fit, parameter_name = "log_lik2", merge_chains = TRUE)
smallPrBdko_loglik <- extract_log_lik(fit, parameter_name = "log_lik4", merge_chains = TRUE)

log_lik_comb <- cbind(largePrBwt_loglik, largePrBdko_loglik, smallPrBwt_loglik, smallPrBdko_loglik)
# optional but recommended
ll_array <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)
loofilename <- paste0("loosave_", modelName, ".rds")
write_rds(loo_loglik, file  = file.path(LooDir, loofilename))

ploocv <- data.frame("Model" = modelName,
                     "LooIC" = loo_loglik$estimates[3],
                     "SE" = loo_loglik$estimates[6], 
                     "PLoo" = loo_loglik$estimates[2])

ploocv

#write.table(ploocv, file = file.path(outputDir, 'timeinfluxfit', "stat_table_MZB1.csv"),
#            sep = ",", append = T, quote = FALSE,
#            col.names = F, row.names = FALSE)


################################################################################################
################################################################################################
## posterior predictive distributions
### parameters table
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
exp(out_table)
write.csv(out_table, file = file.path(outputDir, paste0('params_', modelName, ".csv")))


# time sequence for predictions 
ts_pred <- 10^seq(log10(0.1), log10(30), length.out = 300)
numPred <- length(ts_pred)


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


Y1pred <- as.data.frame(fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("Time_h" = ts_pred,
            "subpop" = "LargePreB",
            "Genotype" = "WT")


Y2pred <- as.data.frame(fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("Time_h" = ts_pred,
            "subpop" = "SmallPreB",
            "Genotype" = "WT")

Y3pred <- as.data.frame(fit, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("Time_h" = ts_pred,
            "subpop" = "LargePreB",
            "Genotype" = "dKO")


Y4pred <- as.data.frame(fit, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("Time_h" = ts_pred,
            "subpop" = "SmallPreB",
            "Genotype" = "dKO")


Ypred <- rbind(Y1pred, Y2pred, Y3pred, Y4pred)

brdu_plot <- brdu_df %>%
  select(sample_id, Time_h, Genotype, BrdU_large_pre_B, BrdU_small_pre_B) %>%
  gather(-c(sample_id, Time_h, Genotype), key = "SubPop", value = "prop_brdu") %>%
  mutate(subpop = ifelse(grepl("large", SubPop), "LargePreB", "SmallPreB"))

brdu_plot$Genotype = factor(brdu_plot$Genotype, levels = c("WT", "dKO"))

ggplot()+
  geom_point(data = brdu_plot, aes(x=Time_h, y=prop_brdu, col=subpop)) +
  geom_line(data = Ypred, aes(x=Time_h, y=median*100, col=subpop)) + 
  geom_ribbon(data = Ypred, aes(x = Time_h, ymin = lb*100, ymax = ub*100, fill=subpop), alpha = 0.25) +
  scale_color_discrete(name = NULL)+
  facet_grid(factor(Genotype) ~ subpop) +
  labs(x = "Time post BrdU injection (hours)", y= paste0("% BrdU+ cells")) +
  #scale_y_continuous(limits=c(1, 3e2), trans = "log10")+
  #scale_x_log10(limits=c(1, 30))
  xlim(0, 31) + ylim(0, 100) + guides(col="none", fill="none")



ggsave(file.path("output_fit/rsquare", paste0("P1_", modelName, ".pdf")), last_plot(), device = "pdf", width = 10, height = 7)


eps_function <- function(Time, r_eps){
  #value = Time * exp(-r_eps * Time^2);
  value = 1.0/(1+((Time-0)/r_eps)^0.5);
  return(value)
}
eps_vec_pred <- sapply(ts_pred, eps_function, r_eps = 3.5)

ggplot()+
  geom_line(aes(x = ts_pred, y=eps_vec_pred), col=4) +
  #scale_x_log10(limits=c(0.1, 30), breaks=c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) +
  labs(x = "Time post BrdU injection (hours)", y =NULL, title = paste0("Labelling efficiency of BrdU")) 

ggsave(file.path("output_fit", paste0("P2_", modelName, ".pdf")), last_plot(), device = "pdf", width = 6, height = 4.5)


### Residual plots
resid_d1  <- t(as.data.frame(fit, pars = "resid_d1"))[,1]
resid_d2  <- t(as.data.frame(fit, pars = "resid_d2"))[,1]
resid_d3  <- t(as.data.frame(fit, pars = "resid_d3"))[,1]
resid_d4  <- t(as.data.frame(fit, pars = "resid_d4"))[,1]


p1.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = fracs_wt, aes(x = Time_h, y = resid_d1), col=6) +
  labs(title=paste("Residuals CAR GC counts WT"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = fracs_wt, aes(x = Time_h, y = resid_d2), col=2) +
  labs(title=paste("Residuals CAR MZ counts WT"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p3.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = fracs_dko, aes(x = Time_h, y = resid_d3), col=4) +
  labs(title=paste("Residuals GC counts N2KO"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p4.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = fracs_dko, aes(x = Time_h, y = resid_d4), col=4) +
  labs(title=paste("Residuals MZ counts N2KO"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

cowplot::plot_grid(p1.1, p2.1, p3.1, p4.1, nrow  = 2)



################################################################################################
## open graphics device 
## saving  plots for quality control 
library(bayesplot)

pdf(file = file.path(outputDir, paste(modelName,"StanPlots%03d.pdf", sep = "")),
    width = 8, height = 5, onefile = FALSE, useDingbats = FALSE)
pairs(fit, pars = parametersToPlot)
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
bayesplot::color_scheme_set(scheme = "viridis")

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

dev.off()


