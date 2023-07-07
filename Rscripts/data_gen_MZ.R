### data wrangling for Marginal Zone compartmnet
rm(list = ls()); gc();

library(tidyverse)

## importing data to be fitted 
counts_file <- file.path("datafiles", "counts_MZ.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K)

Nfd_file <- file.path("datafiles", "Nfd_MZ.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2)

ki_file <- file.path("datafiles", "Ki67_MZ.csv")
ki_data <- read.csv(ki_file) %>% 
  arrange(age.at.S1K) %>%
  rename(host_ki = host_ki67_MZ,
         donor_ki = donor_ki67_MZ)

## pooled data
chimera_data <- full_join(counts_data, Nfd_data, 
                          by = c("Lamis.ID", "age.at.S1K", "age.at.bmt", "days.post.bmt")) %>%
  full_join(ki_data, 
            by = c("Lamis.ID", "age.at.S1K", "age.at.bmt", "days.post.bmt")) %>%
  ### Binning the data for easy predictions
  mutate(age_bins = ifelse(age.at.bmt <= 56, "age_bin1",
                           ifelse(age.at.bmt <= 77, "age_bin2", "age_bin3")))

## defining the function to calculate mode of a vector series
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## mean of age at BMT within each age at BMT bin
chimera_data %>% arrange(age.at.bmt) %>%
  group_by(age_bins) %>%
  summarise("Mode"= getmode(age.at.bmt),
            "Mean"= mean(age.at.bmt))


## Unique time points with indices to map
unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K  ## unique time points in the data
## Map of the unique time points on all the timepoints
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs <- length(data_time_chi)
n_shards <- length(solve_time_chi)
solve_time <- solve_time_chi
time_index <- time_index_chi
dpBMT <- chimera_data$age.at.S1K -chimera_data$age.at.bmt
ageAtBMT <- unique_times_chi$age.at.bmt
counts <- chimera_data$total_counts
Nfd <- chimera_data$Nfd
ki_donor <- chimera_data$donor_ki
ki_host <- chimera_data$host_ki

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(45), log10(750), length.out = 300)
ts_pred2 <- 10^seq(log10(67), log10(750), length.out = 300)
ts_pred3 <- 10^seq(log10(89), log10(750), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(67, 300)
tb_pred3 <- rep(89, 300)
numPred <- length(ts_pred1)


stan_rdump(c("numObs",  "n_shards", "solve_time", "dpBMT", "ageAtBMT", "time_index",
             "counts",  "Nfd", "ki_donor", "ki_host",
             "ts_pred1", "ts_pred2", "ts_pred3",
             "tb_pred1", "tb_pred2", "tb_pred3", "numPred"),
             file = file.path('datafiles', paste0('MZ_data',".Rdump")))


theta_spline <- function(Time){
  exp(14.36) * exp(-0.00186 * (Time))
}

Chi_T1 <- function(Time){
  chiEst = 0.76; qEst = 0.094;
 # if (Time - 10 < 0){
 #   chi = 0;                       
 # } else {
    chi = chiEst * (1 - exp(-qEst * (Time-10)));
 # }
}

eps_donor=0.98
eps_host = 0.94
neutral_mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    tb = Pars[3]
    
    dY1 <- psi * theta_spline(Time) * Chi_T1(Time - tb) - lambda * Y1
    
    dY2 <- psi * theta_spline(Time) * (1 - Chi_T1(Time - tb)) - lambda * Y2
    
    return(list(c(dY1, dY2)))
  })
}

yini <- c(Y1 = 0, Y2= exp(14.7))

out_df <- data.frame()
for (i in 1:length(solve_time)){
  pars <- c(psi = 0.3, lambda = 0.04, ageAtBMT[i])
  
  out_df[i, 1:3] <- data.frame(deSolve::ode(yini, c(40, solve_time[i]), neutral_mod, pars))[2,1:3]
}

logit_transf <- function(x){asin(sqrt(x))}
#fitting log N_total & logit N_fd to the NIMR data
LL_neutral <- function(param) { 
  psi <- param[1]
  lambda  <- param[2]
  
  k  <- length(param)             #number of unknown parameters 
  n <- length(solve_time)           #number of observations in dataset1
  
  
  sol_df <- data.frame()
  pred_counts <- c()
  pred_Nfd <- c()
  sol_Nfd <- c()
  sol_counts <- c()
  for (i in 1:length(solve_time)){
    pars <- c(psi = param[1], lambda = param[2], ageAtBMT[i])
    
    sol_df[i, 1:3] <- data.frame(deSolve::ode(yini, c(40, solve_time[i]), neutral_mod, pars))[2,1:3]
    
    sol_counts[i] = sol_df[time_index[i], 2] + sol_df[time_index[i], 3]
    sol_Nfd[i] = sol_df[time_index[i], 2]/((sol_df[time_index[i], 2] + sol_df[time_index[i], 3]) * Chi_T1(solve_time[time_index[i]] - ageAtBMT[time_index[i]])) 
  }
  
  
  for (i in 1: length(data_time_chi)){
    
    pred_counts[i] <- sol_counts[time_index[i]]
    pred_Nfd[i] <- sol_Nfd[time_index[i]] 
  }

  R1 <- sum((log(counts) - log(pred_counts))^2)              #SSR for dataset1
  R2 <- sum((logit_transf(Nfd) - logit_transf(pred_Nfd))^2)  #SSR for dataset2
  
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <- -(n/2)*log(R1) - (n/2)*log(R2)
  
  aiccd4NIMR <<- -2*logl + 2*k
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

fit_neutral_mod <- optim(par=c(0.3, 0.05), fn=LL_neutral, control = list(trace = 6))
fit_neutral_mod

pred_df1 <- data.frame()
for (i in 1:length(ts_pred1)){
  pars <- c(psi = fit_neutral_mod$par[1], lambda = fit_neutral_mod$par[2], tb_pred1[i])
  
  pred_df1[i, 1:3] <- data.frame(deSolve::ode(yini, c(40, ts_pred1[i]), neutral_mod, pars))[2,1:3]
}

pred_counts1 <- c()
pred_Nfd1 <- c()
for (i in 1:300) {
  pred_counts1[i] <- pred_df1[i, 2] + pred_df1[i, 3] 
  pred_Nfd1[i] <- pred_df1[i, 2]/(( pred_df1[i, 2] + pred_df1[i, 3]) * Chi_T1(ts_pred1[i] - tb_pred1[i]))
}

pred_df2 <- data.frame()
for (i in 1:length(ts_pred2)){
  pars <- c(psi = fit_neutral_mod$par[1], lambda = fit_neutral_mod$par[2], tb_pred2[i])
  
  pred_df2[i, 1:3] <- data.frame(deSolve::ode(yini, c(40, ts_pred2[i]), neutral_mod, pars))[2,1:3]
}

pred_counts2 <- c()
pred_Nfd2 <- c()
for (i in 1:300) {
  pred_counts2[i] <- pred_df2[i, 2] + pred_df2[i, 3] 
  pred_Nfd2[i] <- pred_df2[i, 2]/(( pred_df2[i, 2] + pred_df2[i, 3]) * Chi_T1(ts_pred2[i] - tb_pred2[i]))
}

pred_df3 <- data.frame()
for (i in 1:length(ts_pred3)){
  pars <- c(psi = fit_neutral_mod$par[1], lambda = fit_neutral_mod$par[2], tb_pred3[i])
  
  pred_df3[i, 1:3] <- data.frame(deSolve::ode(yini, c(40, ts_pred3[i]), neutral_mod, pars))[2,1:3]
}

pred_counts3 <- c()
pred_Nfd3 <- c()
for (i in 1:300) {
  pred_counts3[i] <- pred_df3[i, 2] + pred_df3[i, 3] 
  pred_Nfd3[i] <- pred_df3[i, 2]/(( pred_df3[i, 2] + pred_df3[i, 3]) * Chi_T1(ts_pred3[i] - tb_pred3[i]))
}


ggplot() +
  geom_line(aes(x = ts_pred1, y = pred_Nfd1), col=2, size=1.2) +
  geom_line(aes(x = ts_pred2, y = pred_Nfd2), col=4, size=1.2) +
  geom_line(aes(x = ts_pred3, y = pred_Nfd3), col=6, size=1.2) +
  geom_point(aes(x = data_time_chi, y = Nfd), size=2) 

















