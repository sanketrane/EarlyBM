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

ggplot(counts_df, aes(x= Time_h, y=cell_counts, col = Genotype)) +
  geom_point()+ xlim(0, 35)+
  scale_y_log10() + 
  labs(x='Days post trabsfer', y=NULL)+
  facet_wrap(~ factor(subpop, levels = c("Pro_B", "large_pre_B", "small_pre_B",
                                   "BrdU_Pro_B", "BrdU_large_pre_B", "BrdU_small_pre_B"))) +
  myTheme + guides(col="none")

ggplot(counts_df, aes(x= Time_h, y=cell_counts, col = BrdU_status)) +
  geom_point()+ xlim(0, 35)+
  scale_y_log10() + 
  labs(x='Days post trabsfer', y=NULL)+
  facet_grid(Genotype ~ factor(pop_id, levels = c("Pro_B", "large_pre_B", "small_pre_B",
                                         "BrdU_Pro_B", "BrdU_large_pre_B", "BrdU_small_pre_B"))) +
  myTheme 

ggplot(counts_df, aes(x= Time_h, y=cell_counts, col = Genotype)) +
  geom_point()+ xlim(0, 35)+
  scale_y_log10() + 
  geom_hline(yintercept = 510111, col="darkblue")+
  geom_hline(yintercept = 131619, col="darkred")+
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
  geom_point()+ xlim(0, 35) + scale_y_continuous(trans = "log10", limits = c(1, 100))+
  labs(x='Days post trabsfer', y=NULL) +
  facet_wrap(~ factor(subpop, levels = c("BrdU_Pro_B", "BrdU_large_pre_B", "BrdU_small_pre_B"))) +
  myTheme


counts_df %>%
  group_by(Genotype, pop_id, BrdU_status) %>%
  summarise(mean_c = mean(cell_counts))

# ##stan model
# library(rstan)
# stanmodel <- "stan_models/multiplex_model.stan"
# expose_stan_functions(stanmodel)
# 
# solve_time <- c(4, 18, 30)
# init_cond <- c(0, 0, 372151, 0, 0, 83027, 0, 0, 4318182, 0, 0, 607546)
# params <- c(6, 0.2, 0.01, 0.02, 0.1, 0.01, 0.05)
# ypred <- solve_ODE_sys(solve_time, init_cond, parms = params)



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


solve_time <- c(4, 18, 30)
time_index1 <- purrr::map_dbl(fracs_wt$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index2 <- purrr::map_dbl(fracs_dko$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
# 
# ## Data to import in Stan
# numObs1 <- nrow(large_preB_wt)
# numObs2 <- nrow(large_preB_dko)
# BrdU_wt <- large_preB_wt$prop_brdu
# BrdU_dko <- large_preB_dko$prop_brdu
# n_shards <- length(solve_time)
# 
# # time sequence for predictions specific to age bins within the data
# ts_pred <- seq(0, 30, length.out = 300)
# numPred <- length(ts_pred)
# 
# 
# stan_rdump(c("numObs1",  "numObs2", "solve_time", "n_shards", "time_index1", "time_index2",
#              "BrdU_wt",  "BrdU_dko", "ts_pred", "numPred"),
#            file = file.path('datafiles', paste0('BrdU_stanfit',".Rdump")))


##############################################################################################################################
### model
eps_func_exp <- function(Time, r_eps){
  exp(-r_eps * Time)
}

eps_func_hill <- function(Time, r_eps){
  1/(1+ (Time/r_eps)^1)
}

eps_func <- function(Time, r_eps){
  #eps_func_exp(Time, r_eps)
  eps_func_hill(Time, r_eps)
}

asinsq_inv <- function(x){sin(x)^2}

ode_func <- function (Time, y, parms) {
  pro_pos = (131619) #counts_df %>% filter(subpop == "BrdU_Pro_B") %>%
    #summarise(mean(cell_counts))
  pro_neg = (378491.4) #counts_df %>% filter(subpop == "Pro_B") %>%
    #summarise(mean(cell_counts)) - pro_pos
  with(as.list(c(y, parms)),{
    rho = 0.192868923# exp(rho_log)
    delta = 0.078966709 #exp(delta_log)
    rho_dko = 0.090248636 #exp(rho_dko_log)
    delta_dko = 0.078966709
    lambda = exp(lambda_log)
    lambda_dko = exp(lambda_log)
    alpha = 0# exp(alpha_log)
    mu = exp(mu_log)
    mu_dko = exp(mu_log)
    r_eps = 1.64255675#1.64255675
    
    theta = (rho - delta) * (372151) / (pro_neg+pro_pos);
    theta_dko = (rho_dko - delta) * (83027) / (pro_neg+pro_pos);
    #lambda = (alpha - mu) * (4318182) / (372151);
    #lambda_dko = (alpha - mu_dko) * (607546) / (83027);
    
    del = delta - mu
    del_dko = delta_dko - mu_dko
      
    #theta = (rho - delta) * (y1+y2+y3)/(pro_neg+pro_pos)
    #theta_dko = (rho_dko - delta_dko) * (y4+y5+y6)/(pro_neg+pro_pos)
    
    ## in WT
    # L2 pop in large Pre B
    dy1 <- theta * (pro_pos) + rho * eps_func(Time, r_eps) * (2*y2 + y1) - rho * (1-eps_func(Time, r_eps)) * y1 - (mu + del) * y1
    # L1 pop in large Pre B
    dy2 <- theta * (pro_neg) + rho * eps_func(Time, r_eps) * (2*y3 - y2) + 2 * rho * (1-eps_func(Time, r_eps)) * y1 - (mu + del) * y2
    # U pop in large Pre B
    dy3 <- - rho * eps_func(Time, r_eps) * y3 + rho * (1-eps_func(Time, r_eps)) * (y2 + y3) - (mu + del) * y3
    
    ## in DKO
    # L2 pop in large Pre B
    dy4 <-  theta_dko * (pro_pos) + rho_dko * eps_func(Time, r_eps) * (2*y5+ y4) - rho_dko * (1-eps_func(Time, r_eps)) * y4 - (mu_dko + del_dko) * y4
    # L1 pop in large Pre B
    dy5 <- theta_dko * (pro_neg) + rho_dko * eps_func(Time, r_eps) * (2*y6 - y4) + 2 * rho_dko * (1-eps_func(Time, r_eps)) * y4 - (mu_dko + del_dko) * y5
    # U pop in large Pre B
    dy6 <- - rho_dko * eps_func(Time, r_eps) * y6 + rho_dko * (1-eps_func(Time, r_eps)) * (y5 + y6) - (mu_dko + del_dko) * y6
    
    ## in WT
    # L2 pop in small Pre B
    dy7 <- mu * y1 + alpha * eps_func(Time, r_eps) * (2*y8 + y7) - alpha * (1-eps_func(Time, r_eps)) * y7 - lambda * y7
    # L1 pop in small Pre B
    dy8 <- mu * y2 + alpha * eps_func(Time, r_eps) * (2*y9 - y8) + 2 * alpha * (1-eps_func(Time, r_eps)) * y7 - lambda * y8
    # U pop in small Pre B
    dy9 <- mu * y3 - alpha * eps_func(Time, r_eps) * y9 + alpha * (1-eps_func(Time, r_eps)) * (y8 + y9) - lambda * y9
    
    ##in DKO
    # L2 pop in small PreB
    dy10 <- mu_dko * y4 + alpha * eps_func(Time, r_eps) * (2*y11 + y10) - alpha * (1-eps_func(Time, r_eps)) * y10 - lambda_dko * y10
    # L1 pop in small Pre B
    dy11 <- mu_dko * y5 + alpha * eps_func(Time, r_eps) * (2*y12 - y10) + 2 * alpha * (1-eps_func(Time, r_eps)) * y10 - lambda_dko * y11
    # U pop in small Pre B
    dy12 <- mu_dko * y6 - alpha * eps_func(Time, r_eps) * y12 + alpha * (1-eps_func(Time, r_eps)) * (y11 + y12) - lambda_dko * y12

    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10, dy11, dy12))
  })
}

init_cond_wt <- c("y1" = 0.0, "y2" = 0, "y3" = 372151)
init_cond_dko <- c("y4" = 0.0, "y5" = 0, "y6" = 83027)
init_cond_wt2 <- c("y7" = 0.0, "y8" = 0, "y9" = 4318182)
init_cond_dko2 <- c("y10" = 0.0, "y11" = 0, "y12" = 607546)
init_cond <- c(init_cond_wt, init_cond_dko, init_cond_wt2, init_cond_dko2)
params <- c("lambda_log" = -4,  "mu_log" = -4)
asin_transf <- function(x){asin(sqrt(x))}


#model optimization 
LL_func <- function(param, boot_data1, boot_data2, init_conds) {
  k <- length(param)        #number of unknown parameters 
  n1 <-  nrow(boot_data1)
  n2 <-  nrow(boot_data2)
  
  pred_df <- data.frame(ode(y=init_conds, times=c(0, 4, 18, 30), func=ode_func, parms=param))%>%
    mutate(#wt_large_Total = y1+y2+y3,
           #wt_large_prop_brdu = (y1+y2)/wt_large_Total,
           wt_small_Total = y7+y8+y9,
           wt_small_prop_brdu = (y7+y8)/wt_small_Total,
           #dko_large_Total = y4+y5+y6,
           #dko_large_prop_brdu = (y4+y5)/dko_large_Total,
           dko_small_Total = y10+y11+y12,
           dko_small_prop_brdu = (y10+y11)/dko_small_Total
    ) %>%
    filter(time != 0) %>% rename(Time_h = time) %>% select(Time_h, contains('brdu'))
  
  #wt_large <-  boot_data1$prop_large_PreB 
  wt_small <-  boot_data1$prop_small_PreB
  
  #dko_large <-  boot_data2$prop_large_PreB 
  dko_small <-  boot_data2$prop_small_PreB
  
  
  #R1 <- sum((asin_transf(pred_df$wt_large_prop_brdu[time_index1]) - asin_transf(wt_large))^2)
  R2 <- sum((asin_transf(pred_df$wt_small_prop_brdu[time_index1]) - asin_transf(wt_small))^2)
  #R3 <- sum((asin_transf(pred_df$dko_large_prop_brdu[time_index2]) - asin_transf(dko_large))^2)
  R4 <- sum((asin_transf(pred_df$dko_small_prop_brdu[time_index2]) - asin_transf(dko_small))^2)
  
  logl <- -(n1/2)*(log(R2))  - (n2/2)*(log(R4))  
  
  return(-logl)
} 

fit_LL <- optim(par=params, fn=LL_func, 
                boot_data1 = fracs_wt, boot_data2 = fracs_dko, init_conds = init_cond,
                method = "Nelder-Mead",
                control = list(trace = 6, maxit=2000))

fit_LL
par_est <- fit_LL$par
aic_val <- 2* length(fit_LL$par) + 2*fit_LL$value
aic_val
## predictions
pred_df <- data.frame(ode(y=init_cond, times=seq(0, 30, 0.1), func=ode_func, parms=par_est))%>%
  mutate(
    #wt_large_Total = y1+y2+y3,
    #wt_large_prop_brdu = (y1+y2)/wt_large_Total,
    wt_small_Total = y7+y8+y9,
    wt_small_prop_brdu = (y7+y8)/wt_small_Total,
    #dko_large_Total = y4+y5+y6,
    #dko_large_prop_brdu = (y4+y5)/dko_large_Total,
    dko_small_Total = y10+y11+y12,
    dko_small_prop_brdu = (y10+y11)/dko_small_Total
  ) %>%
  filter(time != 0) %>% rename(Time_h = time) %>%
  select(Time_h, contains('prop')) %>%
  gather(-Time_h, key = "genotype", value = "prop_brdu") %>%
  mutate(Genotype = ifelse(grepl("wt", genotype), "WT", "dKO"),
         subpop = ifelse(grepl("large", genotype), "Large_Pre_B", "Small_Pre_B"))

brdu_plot <- brdu_df %>%
  select(sample_id, Time_h, Genotype, BrdU_small_pre_B) %>%
  gather(-c(sample_id, Time_h, Genotype), key = "SubPop", value = "prop_brdu") %>%
  mutate(subpop = ifelse(grepl("large",SubPop), "Large_Pre_B", "Small_Pre_B"))

# ggplot()+
#   geom_point(data = brdu_plot, aes(x=Time_h, y=prop_brdu, col=Genotype))+
#   geom_line(data = pred_df, aes(x=Time_h, y=prop_brdu*100, col=Genotype)) +
#   facet_wrap(.~ subpop)+
#   xlim(0, 31) + ylim(0, 80)

small_preB_df <- brdu_df %>%
  select(sample_id, Time_h, Genotype, BrdU_small_pre_B) %>%
  rename(prop_brdu = BrdU_small_pre_B)


exp(fit_LL$par)
ggplot()+
  geom_point(data = small_preB_df, aes(x=Time_h, y=prop_brdu, col=Genotype))+
  geom_line(data = pred_df, aes(x=Time_h, y=prop_brdu*100, col=Genotype)) +
  facet_wrap(.~ Genotype)+
  labs(x = "Time post BrdU injection (hours)", y= paste0("% BrdU+ cells"))+
  xlim(0, 31) + ylim(0, 60)

modelName <- "M0v2.1"
ggsave(file.path("output_fit", paste0("P3_", modelName, ".pdf")), last_plot(), device = "pdf", width = 9, height = 4.5)

ts_pred <- seq(0.1, 30, 0.05)
eps_vec_pred <- sapply(ts_pred, eps_func, r_eps = 2.5)
ggplot()+
  geom_line(aes(x = ts_pred, y=eps_vec_pred), col=4)+
  scale_x_log10(limits=c(0.1, 30), breaks=c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30))

