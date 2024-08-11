functions {
#include /include/SR.stan
#include /include/gnormal_lpdf.stan
}

data {
  int<lower=1> SR_fun;          // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;               // total number of cases in all pops and years
  array[N] int<lower=1,upper=N> pop;  // population index
  array[N] int<lower=1,upper=N> year; // brood year index
  int<lower=1,upper=N> N_fit;   // number of cases used in fitting (non-missing S_obs and R_obs)
  array[N_fit] int<lower=1,upper=N> which_fit; // cases used in fitting
  vector<lower=0>[N] S_obs;     // observed annual total spawner abundance (not density)
  vector<lower=0>[N] R_obs;     // total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  // habitat area associated with each spawner abundance obs
  array[N] int<lower=0,upper=1> S_NA; // logical indicating whether S_obs is missing and should be simulated
  array[N] int<lower=0,upper=1> R_NA; // logical indicating whether R_obs is missing and should be simulated
  int<lower=2> N_age;           // number of adult age classes
  int<lower=2> max_age;         // maximum adult age
  matrix<lower=0,upper=1>[max(pop),N_age] p_pop_obs; // average recruit age distributions by pop 
}

transformed data {
  int<lower=1,upper=N> N_pop;      // number of populations
  int<lower=1,upper=N> N_year;     // number of years
  array[N_age] int<lower=2> ages;  // adult ages
  real mu_mu_Rmax = max(log(R_obs[which_fit] ./ A[which_fit]));   // prior mean of mu_Rmax
  real sigma_mu_Rmax = sd(log(R_obs[which_fit] ./ A[which_fit])); // prior SD of mu_Rmax
  
  N_pop = max(pop);
  N_year = max(year);
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
}

parameters {
  real mu_alpha;                        // hyper-mean log intrinsic productivity of wild spawners
  real<lower=0> sigma_alpha;            // hyper-SD log intrinsic productivity
  vector[N_pop] zeta_alpha;             // log intrinsic prod of wild spawners (Z-scores)
  real mu_Rmax;                         // hyper-mean log maximum recruitment
  real<lower=0> sigma_Rmax;             // hyper-SD log maximum recruitment
  vector[N_pop] zeta_Rmax;              // log maximum recruitment (Z-scores)
  real<lower=-1,upper=1> rho_alphaRmax; // correlation between log(alpha) and log(Rmax)
  real<lower=-1,upper=1> rho_R;         // AR(1) coef for brood year log productivity anomalies
  real<lower=0> sigma_year_R;           // hyper-SD of brood year log productivity anomalies
  vector[max(year)] zeta_year;          // log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma_R;                // residual error SD
}

transformed parameters {
  vector<lower=0>[N_pop] alpha;       // intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;        // maximum recruitment 
  vector[N_year] eta_year_R;          // log brood year productivity anomalies
  vector<lower=0>[N] R_hat;           // expected recruit abundance (not density) by brood year
  
  // Multivariate Matt trick for [log(alpha), log(Rmax)]
  {
    matrix[2,2] L_alphaRmax;         // Cholesky factor of corr matrix of log(alpha), log(Rmax)
    matrix[N_pop,2] zeta_alphaRmax;  // [log(alpha), log(Rmax)] random effects (z-scored)
    matrix[N_pop,2] eta_alphaRmax;   // [log(alpha), log(Rmax)] random effects
    vector[2] sigma_alphaRmax;       // SD vector of [log(alpha), log(Rmax)]
    
    L_alphaRmax[1,1] = 1;
    L_alphaRmax[2,1] = rho_alphaRmax;
    L_alphaRmax[1,2] = 0;
    L_alphaRmax[2,2] = sqrt(1 - rho_alphaRmax^2);
    sigma_alphaRmax[1] = sigma_alpha;
    sigma_alphaRmax[2] = sigma_Rmax;
    zeta_alphaRmax = append_col(zeta_alpha, zeta_Rmax);
    eta_alphaRmax = diag_pre_multiply(sigma_alphaRmax, L_alphaRmax * zeta_alphaRmax')';
    alpha = exp(mu_alpha + eta_alphaRmax[,1]);
    Rmax = exp(mu_Rmax + eta_alphaRmax[,2]);
  }
  
  // AR(1) model for eta_year_R
  eta_year_R[1] = zeta_year[1]*sigma_year_R/sqrt(1 - rho_R^2); // initial anomaly
  for(i in 2:N_year)
    eta_year_R[i] = rho_R*eta_year_R[i-1] + zeta_year[i]*sigma_year_R;
  eta_year_R = eta_year_R - mean(eta_year_R);  // constrain log anomalies to sum to zero

  // Predict recruitment
  R_hat = rep_vector(0,N);
  for(i in 1:N_fit)
    R_hat[which_fit[i]] = SR(SR_fun, {0,0}, alpha[pop[which_fit[i]]], 0, 0,
                             Rmax[pop[which_fit[i]]], 0, 0,
                             S_obs[which_fit[i]], 0, 0, A[which_fit[i]]);
}

model {
  // Priors
  mu_alpha ~ normal(0,5);
  sigma_alpha ~ gnormal(0,3,10);
  mu_Rmax ~ normal(mu_mu_Rmax, sigma_mu_Rmax);
  sigma_Rmax ~ normal(0,3);
  rho_R ~ gnormal(0,0.85,50);  // mildly regularize to ensure stationarity
  sigma_year_R ~ normal(0,3);
  sigma_R ~ normal(0,2);
  
  // Hierarchical priors
  // [log(alpha), log(Rmax)] ~ MVN(0, D*R_log_aRmax*D), where D = diag_matrix(sigma_alpha, sigma_Rmax)
  zeta_alpha ~ std_normal();
  zeta_Rmax ~ std_normal();
  zeta_year ~ std_normal(); // eta_year_R ~ N(0, sigma_year_R)
  
  // Likelihood
  R_obs[which_fit] ~ lognormal(log(R_hat[which_fit]) + eta_year_R[year[which_fit]], sigma_R);
}

generated quantities {
  vector[N] S_sim;    // simulated spawners
  vector[N] R_sim;    // simulated recruits

  S_sim = S_obs;
  R_sim = R_obs;

  for(i in 1:N)
  {
    if(S_NA[i] == 1)
    {
      if(i >= max_age && pop[i-max_age] == pop[i])
      {
        S_sim[i] = 0;
        for(a in 1:N_age)
          S_sim[i] = S_sim[i] + R_sim[i-ages[a]]*p_pop_obs[pop[i],a];
      }
    }

    if(R_NA[i] == 1)
      R_sim[i] = SR(SR_fun, {0,0}, alpha[pop[i]], 0, 0, Rmax[pop[i]], 0, 0, 
                    S_sim[i], 0, 0, A[i]) * lognormal_rng(eta_year_R[year[i]], sigma_R);
  }
}
