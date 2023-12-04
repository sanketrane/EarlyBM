functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real eps_function(real time, real[] params){
    real r_eps = params[1];
    //real value = exp(-r_eps * time);
    real value = 1/(1.0 + (time/r_eps)^0.5);
    return value;
   }
   
   real proB_frac(real time, real b0, real r1){
    return b0/(1 + exp(-time * r1));
   }

   real expit_func(real z){
    return exp(z)/(1 + exp(z));
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real rho = exp(parms[2]);
     real kappa = parms[3];
     real lambda = exp(parms[4]);
     real lambda_dko = exp(parms[4]);
     real mu = exp(parms[5]);
     real mu_dko = exp(parms[5]);

     real rho_dko = expit_func(kappa) * rho;

     //real x_pos_wt = proB_frac(time, 0.28, 0.17) * (472334);
     //real x_neg_wt = (1 - proB_frac(time, 0.28, 0.17)) * (472334);
     //real x_pos_dko = proB_frac(time, 0.434, 0.027) * (549776);
     //real x_neg_dko = (1 - proB_frac(time, 0.434, 0.027)) * (549776);
     real x_pos = proB_frac(time, 0.28, 0.27) * (510111);
     real x_neg= (1 - proB_frac(time, 0.28, 0.27)) * (510111);
     
     //real theta = ((rho + lambda) * (y[7] + y[8]) - mu * (y[1] + y[2] + y[3]))/(x_pos_wt + x_neg_wt);
     //real theta_dko = ((rho_dko + lambda_dko) * (y[9] + y[10]) - mu * (y[4] + y[5] + y[6]))/(x_pos_dko + x_neg_dko);
     real theta = (lambda - rho) * (y[1] + y[2] + y[3] + y[7] + y[8])/(x_pos + x_neg);
     real theta_dko = (lambda_dko - rho_dko) * (y[4] + y[5] + y[6] + y[9] + y[10])/(x_pos + x_neg);


     // the system of ODEs
     real dydt[10];
     // L2 large pre B cells in WT
     dydt[1] = rho * eps_function(time, parms) * (2 * y[7] + y[2] + y[1]) - rho * (1 - eps_function(time, parms)) * y[1] - (mu + lambda) * y[1];
     // L1 large pre B cells in WT
     dydt[2] = rho * eps_function(time, parms) * (2 * y[8] +  2 * y[3]) + rho * (1 - eps_function(time, parms)) * (2 *  y[1]) - (mu + lambda) * y[2];
     // U large pre B cells in WT
     dydt[3] = - rho * eps_function(time, parms) * y[3] + rho * (1 - eps_function(time, parms)) * (2* y[7] + 2* y[8] + y[2] + y[3]) - (mu + lambda) * y[3];

     // L2 large pre B cells in dko
     dydt[4] = rho_dko * eps_function(time, parms) * (2 * y[9] + y[5] + y[4]) - rho_dko * (1 - eps_function(time, parms)) * y[4] - (mu_dko + lambda_dko) * y[4];
     // L1 large pre B cells in dko
     dydt[5] = rho_dko * eps_function(time, parms) * (2 * y[10] + 2 * y[6]) + rho_dko * (1 - eps_function(time, parms)) * (2 * y[4]) - (mu_dko + lambda_dko) * y[5];
     // U large pre B cells in dko
     dydt[6] = - rho_dko * eps_function(time, parms) * y[6] + rho_dko * (1 - eps_function(time, parms)) * (2* y[9] + 2* y[10] + y[5] + y[6]) - (mu_dko + lambda_dko) * y[6];
     
     // L2 small pre B cells in WT
     dydt[7] = theta * x_pos + mu * (y[1] + y[2]) - (rho + lambda) * y[7];
     // L1 small pre B cells in WT
     dydt[8] = theta * x_neg + mu * (y[3]) - (rho + lambda) * y[8];

     // L2 small pre B cells in dko
     dydt[9] = theta_dko * x_pos + mu_dko * (y[4] + y[5]) - (rho_dko + lambda_dko) * y[9];
     // L1 small pre B cells in dko
     dydt[10] = theta_dko * x_neg + mu_dko * (y[6]) - (rho_dko + lambda_dko) * y[10];

     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 10];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 0.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }

   // functions for transformation of fractions in (0,a), where a >=1
    real[] asinsqrt_array(real[] x){
      int ndims = size(x);
      real answer[ndims];
      real a = 1.0;

      for (i in 1: ndims){
        answer[i] = asin(sqrt(x[i])/sqrt(a));
      }
      return answer;
    }

    real asinsqrt_real(real x){
      real a = 1.0;

      real answer = asin(sqrt(x)/sqrt(a));
      return answer;
    }

    real asinsqrt_inv(real x){
      real a = 1.0;

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
  real<lower = 0> largePreB_fracs_wt[numObs1];
  real<lower = 0> largePreB_fracs_dko[numObs2];
  real<lower = 0> smallPreB_fracs_wt[numObs1];
  real<lower = 0> smallPreB_fracs_dko[numObs2];
  real<lower = 0> largePreB_counts_wt[numObs1];
  real<lower = 0> largePreB_counts_dko[numObs2];
  real<lower = 0> smallPreB_counts_wt[numObs1];
  real<lower = 0> smallPreB_counts_dko[numObs2];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real rho_Log;
  real kappa;
  real lambda_Log;
  real mu_Log;
  //real<lower=0> r_eps;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  //real<lower = 0> sigma5;
  //real<lower = 0> sigma6;
  //real<lower = 0> sigma7;
  //real<lower = 0> sigma8;
  }


transformed parameters{
  real y_hat[n_shards, 10];     // declaring the array for ODE solution
  real largePreB_fracs_wt_mean[numObs1];
  real largePreB_fracs_dko_mean[numObs2];
  real smallPreB_fracs_wt_mean[numObs1];
  real smallPreB_fracs_dko_mean[numObs2];
  real largePreB_counts_wt_mean[numObs1];
  real largePreB_counts_dko_mean[numObs2];
  real smallPreB_counts_wt_mean[numObs1];
  real smallPreB_counts_dko_mean[numObs2];

  real parms[5];                  // declaring the array for parameters
  real init_cond[10];              // declaring the array for state variables

  // initial conditions and parameters
  init_cond[1] = 0.0;     // L2 largePreB cells in WT at t0
  init_cond[2] = 0.0;     // L1 largePreB cells in WT at t0
  init_cond[3] = 372151;  // U largePreB cells in WT at t0
  init_cond[4] = 0.0;     // L2 largePreB cells in WT at t0
  init_cond[5] = 0.0;     // L1 largePreB cells in WT at t0
  init_cond[6] = 83027;   // U largePreB cells in WT at t0
  init_cond[7] = 0.0;     // Labelled smallPreB cells in WT at t0
  init_cond[8] = 4318182;  // Unlabelled smallPreB cells in WT at t0
  init_cond[9] = 0.0;     // Labelled smallPreB cells in dko at t0
  init_cond[10] = 607546;   // Unlabelled smallPreB cells in dko at t0

  parms[1] = 4.0;
  parms[2] = rho_Log;
  parms[3] = kappa;
  parms[4] = lambda_Log;
  parms[5] = mu_Log;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);

  for (i in 1:numObs1){
    largePreB_counts_wt_mean[i] = y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2] + y_hat[time_index1[i], 3];
    smallPreB_counts_wt_mean[i] = y_hat[time_index1[i], 7] + y_hat[time_index1[i], 8];
    largePreB_fracs_wt_mean[i] = (y_hat[time_index1[i], 1] + y_hat[time_index1[i], 2])/largePreB_counts_wt_mean[i];
    smallPreB_fracs_wt_mean[i] = y_hat[time_index1[i], 7]/smallPreB_counts_wt_mean[i];
  }
  for (i in 1:numObs2){
    largePreB_counts_dko_mean[i] = y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5] + y_hat[time_index2[i], 6];
    smallPreB_counts_dko_mean[i] = y_hat[time_index2[i], 9] + y_hat[time_index2[i], 10];
    largePreB_fracs_dko_mean[i] = (y_hat[time_index2[i], 4] + y_hat[time_index2[i], 5])/largePreB_counts_dko_mean[i];
    smallPreB_fracs_dko_mean[i] = y_hat[time_index2[i], 9]/smallPreB_counts_dko_mean[i];
  }
}

model{
  // prior distribution for model parameters
  //r_eps ~ normal(4, 1);
  rho_Log ~ normal(-2, 2);
  kappa ~ normal(-0.5, 1);
  lambda_Log ~ normal(-3, 2);
  mu_Log ~ normal(-4, 2);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);
  sigma4 ~ normal(0, 2.5);
  //sigma5 ~ normal(0, 2.5);
  //sigma6 ~ normal(0, 2.5);
  //sigma7 ~ normal(0, 2.5);
  //sigma8 ~ normal(0, 2.5);

  // model fitting on to data
  asinsqrt_array(largePreB_fracs_wt)  ~ normal(asinsqrt_array(largePreB_fracs_wt_mean), sigma1);
  asinsqrt_array(smallPreB_fracs_wt)  ~ normal(asinsqrt_array(smallPreB_fracs_wt_mean), sigma2);
  asinsqrt_array(largePreB_fracs_dko) ~ normal(asinsqrt_array(largePreB_fracs_dko_mean), sigma3);
  asinsqrt_array(smallPreB_fracs_dko) ~ normal(asinsqrt_array(smallPreB_fracs_dko_mean), sigma4);

  //log(largePreB_counts_wt)  ~ normal(log(largePreB_counts_wt_mean), sigma5);
  //log(smallPreB_counts_wt)  ~ normal(log(smallPreB_counts_wt_mean), sigma6);
  //log(largePreB_counts_dko) ~ normal(log(largePreB_counts_dko_mean), sigma7);
  //log(smallPreB_counts_dko) ~ normal(log(smallPreB_counts_dko_mean), sigma8);
}

generated quantities{
    // ODE predictions
    real y_hat_pred[numPred, 10];
    
    // variables for model predictions
    real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred]; real y4_mean_pred[numPred];
    real y5_mean_pred[numPred]; real y6_mean_pred[numPred]; real y7_mean_pred[numPred]; real y8_mean_pred[numPred];
    
    // variables for model predictions with stdev
    //real largePreB_wt_pred[numPred]; real smallPreB_wt_pred[numPred]; real largePreB_dko_pred[numPred]; real smallPreB_dko_pred[numPred];
    
    // Residuals
    //vector[numObs1] resid_d1; vector[numObs1] resid_d2; vector[numObs2] resid_d3; vector[numObs2] resid_d4; 
    // log likelihoods
    vector[numObs1] log_lik1; vector[numObs1] log_lik2; vector[numObs2] log_lik3; vector[numObs2] log_lik4;
    //vector[numObs1] log_lik5; vector[numObs1] log_lik6; vector[numObs2] log_lik7; vector[numObs2] log_lik8;
    
    //ODE solution
    y_hat_pred[1] = init_cond;
    y_hat_pred[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);
    
    // model predictions with stdev
    for (i in 1:numPred){
      //WT predictions
      //counts
      y5_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2] + y_hat_pred[i, 3]; // large Pre B counts
      y6_mean_pred[i] = y_hat_pred[i, 7] + y_hat_pred[i, 8];  // small Pre B counts
     
      // fracs
      y1_mean_pred[i] = (y_hat_pred[i, 1] + y_hat_pred[i, 2])/y5_mean_pred[i]; // large Pre B fracs
      y2_mean_pred[i] = (y_hat_pred[i, 7])/y6_mean_pred[i]; // small Pre B fracs
      
      //largePreB_wt_pred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y1_mean_pred[i]), sigma1));
      //smallPreB_wt_pred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred[i]), sigma2));
      
      //dKO predictions
      //counts
      y7_mean_pred[i] = y_hat_pred[i, 4] + y_hat_pred[i, 5] + y_hat_pred[i, 6]; // large Pre B counts
      y8_mean_pred[i] = y_hat_pred[i, 9] + y_hat_pred[i, 10];  // small Pre B counts
     
      y3_mean_pred[i] = (y_hat_pred[i, 4] + y_hat_pred[i, 5])/ y7_mean_pred[i]; // large Pre B fracs
      y4_mean_pred[i] = y_hat_pred[i, 9]/ y8_mean_pred[i]; // large Pre B fracs
      
      //largePreB_dko_pred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y3_mean_pred[i]), sigma3));
      //smallPreB_dko_pred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y4_mean_pred[i]), sigma4));
    }
    
    // calculating the log predictive accuracy for each point
    for (n in 1:numObs1) {
      log_lik1[n] = normal_lpdf(asinsqrt_real(largePreB_fracs_wt[n]) | asinsqrt_real(largePreB_fracs_wt_mean[n]), sigma1);
      log_lik2[n] = normal_lpdf(asinsqrt_real(smallPreB_fracs_wt[n]) | asinsqrt_real(smallPreB_fracs_wt_mean[n]), sigma2);
      //log_lik5[n] = normal_lpdf(log(largePreB_counts_wt[n]) | log(largePreB_counts_wt_mean[n]), sigma5);
      //log_lik6[n] = normal_lpdf(log(smallPreB_counts_wt[n]) | log(smallPreB_counts_wt_mean[n]), sigma6);
    }
    
    // calculating the log predictive accuracy for each point
    for (n in 1:numObs2) {
      //resid_d3[n] = asinsqrt_real(largePreB_dko[n]) - asinsqrt_real(largePreB_dko_mean[n]);
      //resid_d4[n] = asinsqrt_real(smallPreB_dko[n]) - asinsqrt_real(smallPreB_dko_mean[n]);
      log_lik3[n] = normal_lpdf(asinsqrt_real(largePreB_fracs_dko[n]) | asinsqrt_real(largePreB_fracs_dko_mean[n]), sigma3);
      log_lik4[n] = normal_lpdf(asinsqrt_real(smallPreB_fracs_dko[n]) | asinsqrt_real(smallPreB_fracs_dko_mean[n]), sigma4);
      //log_lik7[n] = normal_lpdf(log(largePreB_counts_dko[n]) | log(largePreB_counts_dko_mean[n]), sigma7);
      //log_lik8[n] = normal_lpdf(log(smallPreB_counts_dko[n]) | log(smallPreB_counts_dko_mean[n]), sigma8);
    }
}
 