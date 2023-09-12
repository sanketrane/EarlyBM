functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real eps_function(real time, real[] params){
    real r_eps = params[5];
    //real value = exp(-r_eps * time);
    real value = 1.0/(1+(time/r_eps)^2);
    return value;
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real rho = parms[1];
     real delta = parms[2];
     real rho_dko = parms[3];
     real delta_dko = parms[4];

     real x_pos= 131619;
     real x_neg = 378491;

     real phi = (rho - delta) * (y[1] + y[2] + y[3]) / (x_pos + x_neg);
     real phi_dko = (rho_dko - delta_dko) * (y[4] + y[5] + y[6]) / (x_pos + x_neg);

     // the system of ODEs
     real dydt[6];
     // L2 large pre B cells in WT
     dydt[1] =  rho * eps_function(time, parms) * (2 * y[2] + y[1]) - rho * (1 - eps_function(time, parms)) * y[1] - delta * y[1];
     // L1 large pre B cells in WT
     dydt[2] = phi * x_pos + rho * eps_function(time, parms) * (2 * y[3] - y[2]) + rho * (1 - eps_function(time, parms)) * (2 *  y[1]) - delta * y[2];
     // U large pre B cells in WT
     dydt[3] = phi * x_neg - rho * eps_function(time, parms) * y[3] + rho * (1 - eps_function(time, parms)) * (y[2] + y[3]) - delta * y[3];

     // L2 large pre B cells in dko
     dydt[4] =  rho_dko * eps_function(time, parms) * (2 * y[5] + y[4]) - rho_dko * (1 - eps_function(time, parms)) * y[4] - delta_dko * y[4];
     // L1 large pre B cells in dko
     dydt[5] = phi_dko * x_pos + rho_dko * eps_function(time, parms) * (2 * y[6] - y[5]) + rho_dko * (1 - eps_function(time, parms)) * (2 * y[4]) - delta_dko * y[5];
     // U large pre B cells in dko
     dydt[6] = phi_dko * x_neg - rho_dko * eps_function(time, parms) * y[6] + rho_dko * (1 - eps_function(time, parms)) * (y[5] + y[6]) - delta_dko * y[6];

     //// L2 small pre B cells in WT
     //dydt[1] = mu * y[1] + alpha * eps_function(time, parms) * (2 * y[8] + y[7]) - alpha * (1 - eps_function(time, parms)) * y[7] - lambda * y[7];
     //// L1 small pre B cells in WT
     //dydt[2] = mu * y[2] + alpha * eps_function(time, parms) * (2 * y[9] - y[8]) + rho * (1 - eps_function(time, parms)) * (2 *  y[1]) - delta * y[2];
     //// U small pre B cells in WT
     //dydt[3] = phi * x_neg - rho * eps_function(time, parms) * y[3] + rho * (1 - eps_function(time, parms)) * (y[2] + y[3]) - delta * y[3];
//
     //// L2 small pre B cells in dko
     //dydt[4] =  rho_dko * eps_function(time, parms) * (2 * y[5] + y[4]) - rho_dko * (1 - eps_function(time, parms)) * y[4] - delta_dko * y[4];
     //// L1 small pre B cells in dko
     //dydt[5] = phi_dko * x_pos + rho_dko * eps_function(time, parms) * (2 * y[6] - y[5]) + rho_dko * (1 - eps_function(time, parms)) * (2 * y[4]) - delta_dko * y[5];
     //// U small pre B cells in dko
     //dydt[6] = phi_dko * x_neg - rho_dko * eps_function(time, parms) * y[6] + rho_dko * (1 - eps_function(time, parms)) * (y[5] + y[6]) - delta_dko * y[6];
     
     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 6];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 0.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }

   // functions for transformation of fractions in (0,a), where a >=1
   real[] asinsqrt_array(real[] x){
     int ndims = size(x);
     real answer[ndims];
     real a = 1.2;

     for (i in 1: ndims){
       answer[i] = asin(sqrt(x[i])/sqrt(a));
     }
     return answer;
   }

   real asinsqrt_real(real x){
     real a = 1.2;

     real answer = asin(sqrt(x)/sqrt(a));
     return answer;
   }

   real asinsqrt_inv(real x){
     real a = 1.2;

     real answer = a * (sin(x))^2;
     return answer;
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
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real rho_dko_Log;
  real<lower=rho_dko_Log>rho_Log;
  real delta_Log;
  real delta_dko_Log;
  real<lower = 0> r_eps;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  }


transformed parameters{
  real y_hat[n_shards, 6];     // declaring the array for ODE solution
  real largePreB_wt_mean[numObs1];
  real largePreB_dko_mean[numObs2];
  //real smallPreB_wt_mean[numObs1];
  //real smallPreB_dko_mean[numObs2];
  real parms[5];                  // declaring the array for parameters
  real init_cond[6];              // declaring the array for state variables

  // initial conditions and parameters
  init_cond[1] = 0.0;     // L2 cells in WT at t0
  init_cond[2] = 0.0;     // L1 cells in WT at t0
  init_cond[3] = 372151;  // U cell in WT at t0
  init_cond[4] = 0.0;     // L2 cells in dko at t0
  init_cond[5] = 0.0;     // L1 cells in dko at t0
  init_cond[6] = 83027;   // U cell in dko at t0

  
  parms[1] = exp(rho_Log);
  parms[2] = exp(delta_Log);
  parms[3] = exp(rho_dko_Log);
  parms[4] = exp(delta_dko_Log);
  parms[5] = r_eps;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);
  

  for (i in 1:numObs1){
    largePreB_wt_mean[i] = (y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2])/(y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2] + y_hat[time_index1[i], 3]);
  }
  for (i in 1:numObs2){
    largePreB_dko_mean[i] = (y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5])/(y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5] + y_hat[time_index2[i], 6]);
  }
}

model{
  // prior distribution for model parameters
  rho_Log ~ normal(-2, 0.5);
  delta_Log ~ normal(-3, 0.5);
  rho_dko_Log ~ normal(-4, 0.5);
  delta_dko_Log ~ normal(-3, 0.5);
  r_eps ~ normal(5, 1);

  sigma1 ~ normal(0.1, 0.5);
  sigma2 ~ normal(0.1, 0.5);

  // model fitting on to data
  (largePreB_wt) ~ normal((largePreB_wt_mean), sigma1);
  (largePreB_dko) ~ normal((largePreB_dko_mean), sigma2);
}

generated quantities{
   // ODE predictions
   real y_hat_pred[numPred, 6];
   // variables for model predictions
   real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; 
   // variables for model predictions with stdev
   real largePreB_wt_pred[numPred]; real largePreB_dko_pred[numPred]; 
   //real smallPreB_wt_pred[numPred]; real smallPreB_dko_pred[numPred];
   // Residuals
   //vector[numObs1] resid_d1; vector[numObs2] resid_d2; 
   // log likelihoods
   vector[numObs1] log_lik1; vector[numObs2] log_lik2; 

   //ODE solution
   y_hat_pred[1] = init_cond;
   y_hat_pred[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);

   // model predictions with stdev
   for (i in 1:numPred){
     //WT predictions
     y1_mean_pred[i] = (y_hat_pred[i, 1] + y_hat_pred[i, 2])/(y_hat_pred[i, 1] + y_hat_pred[i, 2] + y_hat_pred[i, 3]); // large Pre B
     largePreB_wt_pred[i] = (normal_rng((y1_mean_pred[i]), sigma1));
     //smallPreB_wt_pred[i] = inv_logit(normal_rng(logit(y2_mean_pred[i]), sigma2));
     //dKO predictions
     y2_mean_pred[i] = (y_hat_pred[i, 4] + y_hat_pred[i, 5])/(y_hat_pred[i, 4] + y_hat_pred[i, 5] + y_hat_pred[i, 6]); // large Pre B
     largePreB_dko_pred[i] = (normal_rng((y2_mean_pred[i]), sigma2));
     //smallPreB_dko_pred[i] = inv_logit(normal_rng(logit(y4_mean_pred[i]), sigma4));
   }

    // calculating the log predictive accuracy for each point
    for (n in 1:numObs1) {
    //  resid_d1[n] = logit(largePreB_wt[n]) - logit(largePreB_wt_mean[n]);
      log_lik1[n] = normal_lpdf((largePreB_wt[n]) | (largePreB_wt_mean[n]), sigma1);
    }
    //
    //// calculating the log predictive accuracy for each point
    for (n in 1:numObs2) {
    //  resid_d2[n] = logit(largePreB_dko[n]) - logit(largePreB_dko_mean[n]);
      log_lik2[n] = normal_lpdf((largePreB_dko[n]) | (largePreB_dko_mean[n]), sigma2);
    }
}
