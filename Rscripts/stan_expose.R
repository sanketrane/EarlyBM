rm(list = ls()); gc()

library(rstan)
library(tidyverse)

## model
modelName <- "Simple_multiplex_model"
rstan::expose_stan_functions(file.path("stan_models", paste0(modelName, ".stan")))
par_inc <- c(1.6, 0.2, 0.08, 0.1)

timeseq <- seq(0.1, 30, 0.1)
init_cond <- c("y1" = 0.0, "y2" = 0, "y3" = 372151, "y4" = 0.0, "y5" = 0, "y6" = 83027)
sol_vec <- solve_ODE_sys(timeseq, init_cond, par_inc)




