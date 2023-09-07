rm(list = ls()); gc()

library(rstan)
library(tidyverse)

## model
modelName <- "Simple_multiplex_model"
rstan::expose_stan_functions(file.path("stan_models", paste0(modelName, ".stan")))
params <- c(0.1, 0.01, 0.04)
params1 <- c(0.21, 0.34, 0.18)
set.seed(5454)
timeseq <- seq(0.1, 30, 0.1)
time_set <- round(runif(25, 1, 30), 0) %>% sort() 
time_obs <- time_set %>% unique()
time_set
time_obs
time_index <- purrr::map_dbl(time_set, function(x) which(x == time_obs))
init_cond <- c("y1" = 0.0, "y2" = 0, "y3" = 372151, "y4" = 0.0, "y5" = 0, "y6" = 83027)
ode_df <- solve_ODE_sys(time_obs, init_cond, params)
ode_pred <- solve_ODE_sys(time_obs, init_cond, params1)


noise1 <- rnorm(length(time_set), 0.1, 0.05)
noise2 <- rnorm(length(time_set), 0.15, 0.05)

stan_pred_df <- data.frame("time_h" = time_obs,
                            "y_pred" = matrix(unlist(ode_df), nrow = length(time_obs), byrow = TRUE)) %>%
  mutate(wt_brdu = (y_pred.1 + y_pred.2)/(y_pred.1 + y_pred.2 + y_pred.3),
         dko_brdu = (y_pred.4 + y_pred.5)/(y_pred.4 + y_pred.5 + y_pred.6)) %>%
  select(time_h, wt_brdu, dko_brdu)

fit_pred <- data.frame("time_h" = time_obs,
                           "y_pred" = matrix(unlist(ode_pred), nrow = length(time_obs), byrow = TRUE)) %>%
  mutate(wt_brdu = (y_pred.1 + y_pred.2)/(y_pred.1 + y_pred.2 + y_pred.3),
         dko_brdu = (y_pred.4 + y_pred.5)/(y_pred.4 + y_pred.5 + y_pred.6)) %>%
  select(time_h, wt_brdu, dko_brdu) %>%
  gather(-time_h, key = "Genotype", value = "prop_brdu")

artf_df <- data.frame("time_h" = time_set,
                      wt_brdu = stan_pred_df$wt_brdu[time_index] + noise1,
                      dko_brdu = stan_pred_df$dko_brdu[time_index]+ noise2)

artf_df %>%
  gather(-time_h, key = "Genotype", value = "prop_brdu") %>%
  ggplot(aes(x=time_h, y=prop_brdu, col=Genotype)) +
  geom_point(size=1.2) + facet_wrap(.~ Genotype)+
  geom_line(data = fit_pred) + facet_wrap(.~ Genotype)

write.csv(artf_df, "artf_df.csv", row.names = F)


## Data to import in Stan
numObs1 <- length(time_set)
numObs2 <- length(time_set)
n_shards <- length(time_obs)
largePreB_wt <- artf_df$wt_brdu
largePreB_dko <- artf_df$dko_brdu
time_index1 = time_index
time_index2 = time_index
solve_time = time_obs
# time sequence for predictions specific to age bins within the data
ts_pred <- 10^seq(log10(0.1), log10(30), length.out = 300)
numPred <- length(ts_pred)


stan_rdump(c("numObs1",  "numObs2", "n_shards", "solve_time", "time_index1", "time_index2",
             "largePreB_wt", "largePreB_dko",
             "ts_pred", "numPred"),
           file = file.path('datafiles', paste0('artf_data',".Rdump")))



read_rdump(file.path('datafiles', paste0('artf_data',".Rdump")))


