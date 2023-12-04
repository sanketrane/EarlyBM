rm(list = ls()); gc();

############################################################
### Preamble
############################################################

## loading libraries
library(tidyverse)
library(readxl)
library(deSolve)
library(boot)

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
 library(rstan)
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



proB_frac <- function(Time, b0, r){
  #b0 * exp(r*Time)
  b0/(1 + exp(-Time * r))
}
# 
# nls_proB_wt <- nls((BrdU_Pro_B/100) ~ (proB_frac(Time_h, b0, r)), 
#                    data = fracs_wt,
#                    start = list("b0" = 0.2, "r" = 0.05))
# parwt <- coef(nls_proB_wt)
# 
# nls_proB_dko <- nls((BrdU_Pro_B/100) ~ (proB_frac(Time_h, b0, r)), 
#                    data = fracs_dko,
#                    start = list("b0" = 0.2, "r" = 0.05))
# parko <- coef(nls_proB_dko)
# 
# pred_ts <- seq(0, 30, 0.1)
# proBWT_vec <- sapply(pred_ts, proB_frac, b0 = parwt[1], r=parwt[2])
# proBdko_vec <- sapply(pred_ts, proB_frac,b0 = parko[1], r=parko[2])
# 
# ggplot() +
#   geom_point(data = fracs_wt, aes(x=Time_h, y=BrdU_Pro_B)) +
#   geom_line(aes(x=pred_ts, y= proBWT_vec * 100)) +
#   ylim(0,50)
# 
# ggplot() +
#   geom_point(data = fracs_dko, aes(x=Time_h, y=BrdU_Pro_B), col=2) +
#   geom_line(aes(x=pred_ts, y= proBdko_vec * 100), col=2) +
#   ylim(0,50)

solve_time <- c(4, 18, 30)
time_index1 <- purrr::map_dbl(fracs_wt$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index2 <- purrr::map_dbl(fracs_dko$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
# 
## Data to import in Stan
numObs1 <- nrow(fracs_wt)
numObs2 <- nrow(fracs_dko)
largePreB_wt <- fracs_wt$BrdU_large_pre_B/100
largePreB_dko <- fracs_dko$BrdU_large_pre_B/100
smallPreB_wt <- fracs_wt$BrdU_small_pre_B/100
smallPreB_dko <- fracs_dko$BrdU_small_pre_B/100
n_shards = length(solve_time)

# time sequence for predictions specific to age bins within the data
ts_pred <- 10^seq(log10(0.1), log10(30), length.out = 300)
numPred <- length(ts_pred)


stan_rdump(c("numObs1",  "numObs2", "solve_time", "n_shards", "time_index1", "time_index2",
             "largePreB_wt",  "largePreB_dko", "smallPreB_wt",  "smallPreB_dko", "ts_pred", "numPred"),
           file = file.path('datafiles', paste0('BrdU_stanfit',".Rdump")))

##############################################################################################################################
### model
eps_func_exp <- function(Time, r_eps){
  exp(-r_eps * Time)
}

eps_func_hill <- function(Time, r_eps){
  1/(1+ (Time/r_eps)^0.5)
}

eps_func <- function(Time, r_eps){
  #eps_func_exp(Time, r_eps)
 eps_func_hill(Time, r_eps)
}


asinsq_inv <- function(x){sin(x)^2}
expit_fun <- function(x){exp(x)/(1+exp(x))}
# 
 # counts_df %>% filter(subpop == "Pro_B") %>%
 #   group_by(Genotype) %>%
 #   summarise(mean(cell_counts))
 # 
 # counts_df %>% filter(subpop == "BrdU_Pro_B") %>%
 #   group_by(Genotype) %>%
 #   summarise(mean(cell_counts))

counts_df %>% filter(subpop == "small_pre_B") %>%
     group_by(Genotype) %>%
     summarise(mean(cell_counts))

ode_func <- function (Time, y, parms) {
  pro_pos_wt = proB_frac(Time, 0.28, 0.17)  * (472334) 
  pro_neg_wt =  (1-proB_frac(Time, 0.28, 0.17))  * (472334) 
  pro_pos_dko = proB_frac(Time, 0.434, 0.027)  * (549776) 
  pro_neg_dko =  (1-proB_frac(Time, 0.434, 0.027))  * (549776) 
  
  with(as.list(c(y, parms)),{
    rho = exp(rho_log) 
    rho_dko = rho * (expit_fun(kappa)) 
    lambda = exp(lambda_log)
    lambda_dko = exp(lambda_log)
    mu = exp(mu_log)
    mu_dko = exp(mu_dko_log)
    r_eps =  4
    
    theta = (rho - lambda) * (y7+y8) / (pro_neg_wt+pro_pos_wt);
    theta_dko = (rho_dko - lambda_dko) * (y9+y10) / (pro_neg_dko+pro_pos_dko);
    
    ## in WT
    # L2 pop in large Pre B
    dy1 <-  rho * eps_func(Time, r_eps) * (2*y7 + 2*y2 + y1) - rho * (1-eps_func(Time, r_eps)) * y1 - (lambda + mu) * y1
    # L1 pop in large Pre B
    dy2 <-  rho * eps_func(Time, r_eps) * (2*y8 + 2*y3 - y2) + rho * (1-eps_func(Time, r_eps)) * 2*y1 - (lambda + mu) * y2
    # U pop in large Pre B
    dy3 <-  - rho * eps_func(Time, r_eps) * y3 + rho * (1-eps_func(Time, r_eps)) * (2*y7 + 2*y8 + y2 + y3) - (lambda + mu) * y3
    
    ## in DKO
    # L2 pop in large Pre B
    dy4 <-  rho_dko * eps_func(Time, r_eps) * (2*y9 + 2*y5+ y4) - rho_dko * (1-eps_func(Time, r_eps)) * y4 - (lambda_dko + mu_dko) * y4
    # L1 pop in large Pre B
    dy5 <-  rho_dko * eps_func(Time, r_eps) * (2*y10 + 2*y6 - y5) + rho_dko * (1-eps_func(Time, r_eps)) * 2*y4 - (lambda_dko + mu_dko) * y5
    # U pop in large Pre B
    dy6 <-  - rho_dko * eps_func(Time, r_eps) * y6 + rho_dko * (1-eps_func(Time, r_eps)) * (2*y9 + 2*y10 + y5 + y6) - (lambda_dko + mu_dko) * y6
    
    ## WT
    ## Labelled in small preB assuming small PreB dont deivide -- alpha = 0
    dy7 <- theta * (pro_pos_wt) + mu * (y1 + y2) - (rho + lambda) * y7
    ## Unlablled in small preB assuming small PreB dont deivide -- alpha = 0
    dy8 <- theta * (pro_neg_wt)  + mu * (y3) - (rho + lambda) * y8
    
    ## dKO
    ## Labelled in small preB assuming small PreB dont deivide -- alpha = 0
    dy9 <- theta_dko * (pro_pos_dko)  + mu_dko * (y4 + y5) - (rho_dko + lambda_dko) * y9
    ## Unlablled in small preB assuming small PreB dont deivide -- alpha = 0
    dy10 <- theta_dko * (pro_neg_dko) + mu_dko * (y6) - (rho_dko + lambda_dko) * y10
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10))
  })
}

init_cond_wt <- c("y1" = 0.0, "y2" = 0, "y3" = 372151)
init_cond_dko <- c("y4" = 0.0, "y5" = 0, "y6" = 83027)
init_cond_wt2 <- c("y7" = 0.0, "y8" = 4318182)
init_cond_dko2 <- c("y9" = 0.0, "y10" = 607546)
init_cond <- c(init_cond_wt, init_cond_dko, init_cond_wt2, init_cond_dko2)
params <- c("rho_log" = -2,  "kappa" = -2, "lambda_log" = -5, "mu_log"  = -5, "mu_dko_log" = -5)
asin_transf <- function(x){asin(sqrt(x))}

boot_data <- rbind(fracs_wt, fracs_dko)

#model optimization 
LL_func <- function(param, boot_data1, boot_data2) {
  k <- length(param)        #number of unknown parameters 
  n1 <- 21
  n2 <- 19
  
  #boot_data1 <- data_tab[1:n1, ]
  #boot_data2 <- data_tab[n1+1:n2, ]
  
  pred_df <- data.frame(ode(y=init_cond, times=c(0, 4, 18, 30), func=ode_func, parms=param))%>%
    mutate(wt_large_Total = y1+y2+y3,
           wt_large_prop_brdu = (y1+y2)/wt_large_Total,
           wt_small_Total = y7+y8,
           wt_small_prop_brdu = (y7)/wt_small_Total,
           dko_large_Total = y4+y5+y6,
           dko_large_prop_brdu = (y4+y5)/dko_large_Total,
           dko_small_Total = y9+y10,
           dko_small_prop_brdu = (y9)/dko_small_Total
    ) %>%
    filter(time != 0) %>% rename(Time_h = time) %>% select(Time_h, contains('brdu'))
  
  wt_large <-  boot_data1$prop_large_PreB 
  wt_small <-  boot_data1$prop_small_PreB
  
  dko_large <-  boot_data2$prop_large_PreB 
  dko_small <-  boot_data2$prop_small_PreB
  
  
  R1 <- sum((asin_transf(pred_df$wt_large_prop_brdu[time_index1]) - asin_transf(wt_large))^2)
  R2 <- sum((asin_transf(pred_df$wt_small_prop_brdu[time_index1]) - asin_transf(wt_small))^2)
  R3 <- sum((asin_transf(pred_df$dko_large_prop_brdu[time_index2]) - asin_transf(dko_large))^2)
  R4 <- sum((asin_transf(pred_df$dko_small_prop_brdu[time_index2]) - asin_transf(dko_small))^2)
  
  logl <- -(n1/2)*(log(R1))  - (n2/2)*(log(R2)) - (n1/2)*(log(R3))  - (n2/2)*(log(R4))  
  
  return(-logl)
} 

fit_LL <- optim(par=params, fn=LL_func, 
                boot_data1 =  fracs_wt, boot_data2 = fracs_dko,
                method = "Nelder-Mead",
                control = list(trace = 6, maxit=2000))

fit_LL
par_est <- fit_LL$par
aic_val <- 2* length(fit_LL$par) + 2*fit_LL$value
aic_val


## predictions
par_est <- c("rho_log" = -2.16, "kappa" = -0.58, "lambda_log" = -3.01, "mu_log" = -1.96, "mu_dko_log" = -0.8)
pred_df <- data.frame(ode(y=init_cond, times=seq(0, 30, 0.1), func=ode_func, parms=par_est))%>%
  mutate(
    wt_large_Total = y1+y2+y3,
    wt_large_prop_brdu = (y1+y2)/wt_large_Total,
    wt_small_Total = y7+y8,
    wt_small_prop_brdu = (y7)/wt_small_Total,
    dko_large_Total = y4+y5+y6,
    dko_large_prop_brdu = (y4+y5)/dko_large_Total,
    dko_small_Total = y10+y9,
    dko_small_prop_brdu = (y9)/dko_small_Total
  ) %>%
  filter(time != 0) %>% rename(Time_h = time) %>%
  select(Time_h, contains('prop')) %>%
  gather(-Time_h, key = "genotype", value = "prop_brdu") %>%
  mutate(Genotype = ifelse(grepl("wt", genotype), "WT", "dKO"),
         subpop = ifelse(grepl("large", genotype), "Large Pre B", "Small Pre B"))

brdu_plot <- brdu_df %>%
  select(sample_id, Time_h, Genotype, BrdU_large_pre_B, BrdU_small_pre_B) %>%
  gather(-c(sample_id, Time_h, Genotype), key = "SubPop", value = "prop_brdu") %>%
  mutate(subpop = ifelse(grepl("large", SubPop), "Large Pre B", "Small Pre B"))


1/exp(fit_LL$par)
aic_val
par_est
exp(par_est["lambda_log"])
exp(par_est["mu_log"])
exp(par_est["mu_dko_log"])
exp(par_est["rho_log"])
expit_fun(par_est["kappa"])

ggplot()+
  geom_point(data = brdu_plot, aes(x=Time_h, y=prop_brdu, col=subpop)) +
  geom_line(data = pred_df, aes(x=Time_h, y=prop_brdu*100, col=subpop)) +
  scale_color_discrete(name = NULL)+
  facet_grid(subpop ~ Genotype) +
  labs(x = "Time post BrdU injection (hours)", y= paste0("% BrdU+ cells")) +
  #scale_y_log10()
  xlim(0, 31) + ylim(0, 100)

modelName <- "MV2.05_diffMu"
ggsave(file.path("output_fit/combined", paste0("P1_", modelName, ".pdf")), last_plot(), device = "pdf", width = 9, height = 4.5)

ts_pred <- seq(0.1, 30, 0.05)
eps_vec_pred <- sapply(ts_pred, eps_func, r_eps = 4)
ggplot()+
  geom_line(aes(x = ts_pred, y=eps_vec_pred), col=4) +
  #scale_x_log10(limits=c(0.1, 30), breaks=c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) +
  labs(x = "Time post BrdU injection (hours)", y =NULL, title = paste0("Labelling efficiency of BrdU")) 

ggsave(file.path("output_fit/combined", paste0("P2_", modelName, ".pdf")), last_plot(), device = "pdf", width = 6, height = 4.5)
