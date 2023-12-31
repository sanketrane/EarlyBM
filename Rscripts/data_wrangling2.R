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

# stan model
# library(rstan)
# stanmodel <- "stan_models/multiplex_model.stan"
# expose_stan_functions(stanmodel)
# 
# params <- c(0.3, 0.2, 0.1, 0.1)
# ypred <- solve_ODE_sys(solve_time, init_cond, parms = params)



### data munging
large_preB_df <- counts_df %>% filter(pop_id == "large_pre_B") %>% 
  select(sample_id, Time_h, cell_counts, BrdU_status, Genotype) %>%
  spread(key = BrdU_status, value = cell_counts) %>%
  mutate(prop_brdu = BrdU_pos/Total) 

large_preB_wt <- large_preB_df %>% filter(Genotype == "WT") %>% arrange(Time_h)
large_preB_dko <- large_preB_df %>% filter(Genotype == "dKO") %>% arrange(Time_h)

logit_func <- function(x){log(x/(1-x))}
logit_func(large_preB_dko$prop_brdu)

solve_time <- c(4, 18, 30)
time_index1 <- purrr::map_dbl(large_preB_wt$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index2 <- purrr::map_dbl(large_preB_dko$Time_h, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
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
eps_func <- function(Time, r_eps){
  exp(-r_eps * Time)
}
ts_pred <- seq(0, 35, by=0.1)
eps_vec <- sapply(ts_pred, eps_func, r_eps = 0.05)
qplot(x = ts_pred, y=eps_vec)

ode_func <- function (Time, y, parms) {
  pro_pos = (131619) #counts_df %>% filter(subpop == "BrdU_Pro_B") %>%
    #summarise(mean(cell_counts))
  pro_neg = (378491.4) #counts_df %>% filter(subpop == "Pro_B") %>%
    #summarise(mean(cell_counts)) - pro_pos
  
  with(as.list(c(y, parms)),{
    delta = exp(delta_log)
    rho = exp(rho_log)
    #rho_dko = exp(rho_dko_log)
    
    theta = (rho - delta) * (y1+y2+y3)/(pro_neg+pro_pos)
    #theta_dko = (rho_dko - delta) * (y4+y5+y6)/(pro_neg+pro_pos)
    alpha = 1.0
    
    ## in WT
    # L2 pop
    dy1 <- theta * (1 - alpha) * pro_pos + rho * eps_func(Time, r_eps) * (2*y2 + y1) - rho * (1-eps_func(Time, r_eps)) * y1 - delta * y1
    # L1 pop
    dy2 <- theta * alpha * pro_neg + rho * eps_func(Time, r_eps) * (2*y3 - y1) + 2 * rho * (1-eps_func(Time, r_eps)) * y1 - delta * y2
    # U pop
    dy3 <- theta * pro_neg - rho * eps_func(Time, r_eps) * y3 + rho * (1-eps_func(Time, r_eps)) * (y2 + y3) - delta * y3
    
    # ## in DKO
    # # L2 pop
    # dy4 <- theta_dko * (1 - alpha) * pro_pos + rho_dko * eps_func(Time) * (2*y5+ y4) - rho_dko * (1-eps_func(Time)) * y4 - delta * y4
    # # L1 pop
    # dy5 <- theta_dko * alpha * pro_neg + rho_dko * eps_func(Time) * (2*y6 - y4) + 2 * rho_dko * (1-eps_func(Time)) * y4 - delta * y5
    # # U pop
    # dy6 <- theta_dko * pro_neg - rho_dko * eps_func(Time) * y6 + rho_dko * (1-eps_func(Time)) * (y5 + y6) - delta * y6
    # 
    list(c(dy1, dy2, dy3))
  })
}

init_cond_wt <- c("y1" = 0.0, "y2" = 0, "y3" = 372151)
init_cond_dko <- c("y1" = 0.0, "y2" = 0, "y3" = 83027)
params <- c("rho_log" = -2, "delta_log" = -3, "r_eps" = 0.1)
asin_transf <- function(x){asin(sqrt(x))}

pred_df <- data.frame(ode(y=init_cond_wt, times=c(0, 4, 18, 30), func=ode_func, parms=params))%>%
  mutate(wt_Total = y1+y2+y3,
         wt_prop_brdu = (y1+y2)/wt_Total) %>%
  filter(time != 0) %>% rename(Time_h = time)


#model optimization 
LL_func <- function(param, boot_data, init_conds) {
  rho    <- param[1]             #parametrs to be estimated as part of a vector
  delta  <- param[2]
  r_eps  <- param[3]
  
  k <- length(param)        #number of unknown parameters 
  n1 <-  nrow(boot_data)
  
  pred_df <- data.frame(ode(y=init_conds, times=c(0, 4, 18, 30), func=ode_func, parms=param))%>%
    mutate(wt_Total = y1+y2+y3,
           wt_prop_brdu = (y1+y2)/wt_Total) %>%
    filter(time != 0) %>% rename(Time_h = time)
  
  pred_wt <- pred_df %>%
    select(Time_h, contains('wt'))
  
  wt_total <-  boot_data$Total 
  wt_frac <-  boot_data$prop_brdu
  
  R1 <- sum((asin_transf(pred_wt$wt_prop_brdu[time_index1]) - asin_transf(wt_frac))^2)
  
  logl <- -(n1/2)*(log(R1)) 
  
  return(-logl)
} 

fit_LL <- optim(par=params, fn=LL_func, 
                boot_data = large_preB_wt, init_conds = init_cond_wt,
                method = "Nelder-Mead",
                control = list(trace = 6))

fit_LL_dko <- optim(par=params, fn=LL_func, 
                boot_data = large_preB_dko, init_conds = init_cond_dko,
                method = "Nelder-Mead",
                control = list(trace = 6))

fit_LL
par_est <- fit_LL$par
exp(fit_LL$par)

fit_LL_dko
par_est_dko <- fit_LL_dko$par
exp(fit_LL_dko$par)


## predictions
pred_wt <- data.frame(ode(y=init_cond_wt, times=seq(0, 30, by=0.01), func=ode_func, parms=par_est))%>%
  mutate(Total = y1+y2+y3,
         prop_brdu = (y1+y2)/Total) %>% rename(Time_h = time) %>%
  select(Time_h, contains('prop'))

pred_dko <- data.frame(ode(y=init_cond_dko, times=seq(0, 30, by=0.01), func=ode_func, parms=par_est_dko))%>%
  mutate(Total = y1+y2+y3,
         prop_brdu = (y1+y2)/Total) %>% rename(Time_h = time) %>%
  select(Time_h, contains('prop'))

pred_df <- full_join(pred_wt, pred_dko, by = "Time_h", suffix = c(".WT", ".dKO")) %>%
  gather(-Time_h, key = "genotype", value = "prop_brdu") %>%
  mutate(Genotype = ifelse(grepl("WT", genotype), "WT", "dKO"))


ggplot()+
  geom_point(data = large_preB_df, aes(x=Time_h, y=prop_brdu *100, col=Genotype))+
  geom_line(data = pred_df, aes(x=Time_h, y=prop_brdu*100)) +
  facet_wrap(.~ Genotype)+
  xlim(0, 31) + ylim(0, 100)

