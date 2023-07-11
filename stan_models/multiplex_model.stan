functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real eps_function(real time, real[] params){
    real r_eps = params[1];
    //real value = exp(-r_eps * time);
    real value = 1.0/(1+(time/r_eps)^1);
    return value;
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real rho = parms[2];
     real delta = parms[3];
     real rho_dko = parms[4];
     real alpha = parms[5];
     real lambda = parms[6];
     real alpha_dko = parms[7];

     real x_pos= 131619;
     real x_neg = 378491.4;

    //  real phi = (rho - delta) * (y[1]+y[2]+y[3]) / (x_pos + x_neg);
    //  real phi_dko = (rho_dko - delta) * (y[4]+y[5]+y[6]) / (x_pos + x_neg);
    //  real mu = (alpha - lambda) * (y[7]+y[8]+y[9]) / (y[1]+y[2]+y[3]);
    //  real mu_dko = (alpha_dko - lambda) * (y[10]+y[11]+y[12]) / (y[4]+y[5]+y[6]);
     real phi = (rho - delta) * (372151) / (x_pos + x_neg);
     real phi_dko = (rho_dko - delta) * (83027) / (x_pos + x_neg);
     real mu = (alpha - lambda) * (4318182) / (372151);
     real mu_dko = (alpha_dko - lambda) * (607546) / (83027);
     real psi = 0.0;

     // the system of ODEs
     real dydt[12];
     // L2 large pre B cells in WT
     dydt[1] = phi * psi * x_pos + rho * eps_function(time, parms) * (2 * y[2] + y[1]) - rho * (1 - eps_function(time, parms)) * y[1] - (mu + delta) * y[1];
     // L1 large pre B cells in WT
     dydt[2] = phi * (1 - psi) * x_pos + rho * eps_function(time, parms) * (2 * y[3] - y[2]) + rho * (1 - eps_function(time, parms)) * (2 *  y[1]) - (mu + delta) * y[2];
     // U large pre B cells in WT
     dydt[3] = phi * x_neg - rho * eps_function(time, parms) * y[3] + rho * (1 - eps_function(time, parms)) * (y[2] + y[3]) - (mu + delta) * y[3];

     // L2 large pre B cells in dko
     dydt[4] = phi_dko * psi * x_pos + rho_dko * eps_function(time, parms) * (2 * y[5] + y[4]) - rho_dko * (1 - eps_function(time, parms)) * y[4] - (mu + delta) * y[4];
     // L1 large pre B cells in dko
     dydt[5] = phi_dko * (1 - psi) * x_pos + rho_dko * eps_function(time, parms) * (2 * y[6] - y[5]) + rho_dko * (1 - eps_function(time, parms)) * (2 * y[4]) - (mu + delta) * y[5];
     // U large pre B cells in dko
     dydt[6] = phi_dko * x_neg - rho_dko * eps_function(time, parms) * y[6] + rho_dko * (1 - eps_function(time, parms)) * (y[5] + y[6]) - (mu + delta) * y[6];
     
     // L2 small pre B cells in WT
     dydt[7] = mu * y[1] + alpha * eps_function(time, parms) * (2 * y[8] + y[7]) - alpha * (1 - eps_function(time, parms)) * y[7] - lambda * y[7];
     // L1 small pre B cells in WT
     dydt[8] = mu * y[2] + alpha * eps_function(time, parms) * (2 * y[9] - y[8]) + alpha * (1 - eps_function(time, parms)) * (2 *  y[7]) - lambda * y[8];
     // U small pre B cells in WT
     dydt[9] = mu * y[3] - alpha * eps_function(time, parms) * y[9] + alpha * (1 - eps_function(time, parms)) * (y[8] + y[9]) - lambda * y[9];

     // L2 small pre B cells in dko
     dydt[10] = mu_dko * y[4] + alpha_dko * eps_function(time, parms) * (2 * y[11] + y[10]) - alpha_dko * (1 - eps_function(time, parms)) * y[10] - lambda * y[10];
     // L1 small pre B cells in dko
     dydt[11] = mu_dko * y[5] + alpha_dko * eps_function(time, parms) * (2 * y[12] - y[11]) + alpha_dko * (1 - eps_function(time, parms)) * (2 * y[10]) - lambda * y[11];
     // U small pre B cells in dko
     dydt[12] = mu_dko * y[6] - alpha_dko * eps_function(time, parms) * y[12] + alpha_dko * (1 - eps_function(time, parms)) * (y[11] + y[12]) - lambda * y[12];
     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 12];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 0.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}

data{
  int<lower  = 1> numObs1;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> numObs2;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> n_shards;
  int<lower  = 1> numPred;
  int<lower  = 1> time_index1[numObs1];
  int<lower  = 1> time_index2[numObs2];
  real<lower = 0> solve_time[n_shards];
  real<lower = 0> largePreB_wt[numObs1];
  real<lower = 0> largePreB_dko[numObs2];
  real<lower = 0> smallPreB_wt[numObs1];
  real<lower = 0> smallPreB_dko[numObs2];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> rho;
  real<lower = 0, upper=rho> rho_dko;
  real<lower = 0, upper=rho_dko> delta;
  real<lower = 0> r_eps;
  real<lower = 0> alpha;
  real<lower = 0, upper=alpha> alpha_dko;
  real<lower = 0, upper=alpha_dko> lambda;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }


transformed parameters{
  real y_hat[n_shards, 12];     // declaring the array for ODE solution
  real largePreB_wt_mean[numObs1];
  real largePreB_dko_mean[numObs2];
  real smallPreB_wt_mean[numObs1];
  real smallPreB_dko_mean[numObs2];
  real parms[7];                  // declaring the array for parameters
  real init_cond[12];              // declaring the array for state variables

  // initial conditions and parameters
  init_cond[1] = 0.0;     // L2 cells in WT at t0
  init_cond[2] = 0.0;     // L1 cells in WT at t0
  init_cond[3] = 372151;  // U cell in WT at t0
  init_cond[4] = 0.0;     // L2 cells in WT at t0
  init_cond[5] = 0.0;     // L1 cells in WT at t0
  init_cond[6] = 83027;   // U cell in WT at t0
  init_cond[7] = 0.0;     // L2 cells in WT at t0
  init_cond[8] = 0.0;     // L1 cells in WT at t0
  init_cond[9] = 4318182;  // U cell in WT at t0
  init_cond[10] = 0.0;     // L2 cells in WT at t0
  init_cond[11] = 0.0;     // L1 cells in WT at t0
  init_cond[12] = 607546;   // U cell in WT at t0

  parms[1] = r_eps;
  parms[2] = rho;
  parms[3] = delta;
  parms[4] = rho_dko;
  parms[5] = alpha;
  parms[6] = lambda;
  parms[7] = alpha_dko;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);

  for (i in 1:numObs1){
    largePreB_wt_mean[i] = (y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2])/(y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2] + y_hat[time_index1[i], 3]);
    smallPreB_wt_mean[i] = (y_hat[time_index1[i], 7] + y_hat[time_index1[i], 8])/(y_hat[time_index1[i], 7] + y_hat[time_index1[i], 8] + y_hat[time_index1[i], 9]);
  }
  for (i in 1:numObs2){
    largePreB_dko_mean[i] = (y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5])/(y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5] + y_hat[time_index2[i], 6]);
    smallPreB_dko_mean[i] = (y_hat[time_index2[i], 10] + y_hat[time_index2[i], 11])/(y_hat[time_index2[i], 10] + y_hat[time_index2[i], 11] + y_hat[time_index2[i], 12]);
  }
}

model{
  // prior distribution for model parameters
  r_eps ~ normal(5, 1);
  rho ~ normal(0.2, 0.5);
  delta ~ normal(0.1, 0.5);
  rho_dko ~ normal(0.1, 0.5);
  alpha ~ normal(0.2, 0.5);
  lambda ~ normal(0.1, 0.5);
  alpha_dko ~ normal(0.1, 0.5);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);
  sigma4 ~ normal(0, 2.5);

  // model fitting on to data
  (largePreB_wt) ~ normal((largePreB_wt_mean), sigma1);
  (smallPreB_wt) ~ normal((smallPreB_wt_mean), sigma2);
  (largePreB_dko) ~ normal((largePreB_dko_mean), sigma3);
  (smallPreB_dko) ~ normal((smallPreB_dko_mean), sigma4);
}

generated quantities{
   // ODE predictions
   real y_hat_pred[numPred, 12];
   // variables for model predictions
   real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred]; real y4_mean_pred[numPred];
   // variables for model predictions with stdev
   real largePreB_wt_pred[numPred]; real smallPreB_wt_pred[numPred]; real largePreB_dko_pred[numPred]; real smallPreB_dko_pred[numPred];
   // Residuals
   vector[numObs1] resid_d1; vector[numObs1] resid_d2; vector[numObs2] resid_d3; vector[numObs2] resid_d4; 
   // log likelihoods
   vector[numObs1] log_lik1; vector[numObs1] log_lik2; vector[numObs2] log_lik3; vector[numObs2] log_lik4;

   //ODE solution
   y_hat_pred[1] = init_cond;
   y_hat_pred[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);

   // model predictions with stdev
   for (i in 1:numPred){
     //WT predictions
     y1_mean_pred[i] = (y_hat_pred[i, 1] + y_hat_pred[i, 2])/(y_hat_pred[i, 1] + y_hat_pred[i, 2] + y_hat_pred[i, 3]); // large Pre B
     y2_mean_pred[i] = (y_hat_pred[i, 7] + y_hat_pred[i, 8])/(y_hat_pred[i, 7] + y_hat_pred[i, 8] + y_hat_pred[i, 9]); // small Pre B
     largePreB_wt_pred[i] = inv_logit(normal_rng(logit(y1_mean_pred[i]), sigma1));
     smallPreB_wt_pred[i] = inv_logit(normal_rng(logit(y2_mean_pred[i]), sigma2));
     //dKO predictions
     y3_mean_pred[i] = (y_hat_pred[i, 4] + y_hat_pred[i, 5])/(y_hat_pred[i, 4] + y_hat_pred[i, 5] + y_hat_pred[i, 6]); // large Pre B
     y4_mean_pred[i] = (y_hat_pred[i, 10] + y_hat_pred[i, 11])/(y_hat_pred[i, 10] + y_hat_pred[i, 5] + y_hat_pred[i, 12]); // small Pre B
     largePreB_dko_pred[i] = inv_logit(normal_rng(logit(y3_mean_pred[i]), sigma3));
     smallPreB_dko_pred[i] = inv_logit(normal_rng(logit(y4_mean_pred[i]), sigma4));
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs1) {
     resid_d1[n] = logit(largePreB_wt[n]) - logit(largePreB_wt_mean[n]);
     resid_d1[n] = logit(smallPreB_wt[n]) - logit(smallPreB_wt_mean[n]);
     log_lik1[n] = normal_lpdf(logit(largePreB_wt[n]) | logit(largePreB_wt_mean[n]), sigma1);
     log_lik1[n] = normal_lpdf(logit(smallPreB_wt[n]) | logit(smallPreB_wt_mean[n]), sigma2);
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs2) {
     resid_d1[n] = logit(largePreB_dko[n]) - logit(largePreB_dko_mean[n]);
     log_lik1[n] = normal_lpdf(logit(largePreB_dko[n]) | logit(largePreB_dko_mean[n]), sigma3);
     resid_d1[n] = logit(smallPreB_dko[n]) - logit(smallPreB_dko_mean[n]);
     log_lik1[n] = normal_lpdf(logit(smallPreB_dko[n]) | logit(smallPreB_dko_mean[n]), sigma4);
   }
}
