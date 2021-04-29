functions {
  // spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      // discrete exponential
      R = alpha*S;
    else if(SR_fun == 2) // Beverton-Holt
      R = alpha*S/(1 + alpha*S/(A*Rmax));
    else if(SR_fun == 3) // Ricker
      R = alpha*S*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  // Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
  
  // Left multiply vector by matrix
  // works even if size is zero
  vector mat_lmult(matrix X, vector v) {
    vector[rows(X)] Xv;
    Xv = rows_dot_product(X, rep_matrix(to_row_vector(v), rows(X)));
    return(Xv); 
  }
  
  // Quantiles of a vector
  real quantile(vector v, real p) {
    int N = num_elements(v);
    real Np = round(N*p);
    real q;
    
    for(i in 1:N) {
      if(i - Np == 0.0) q = v[i];
    }
    return(q);
  }
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // calendar year identifier
  // fecundity
  int<lower=0> N_E;                    // number of fecundity observations
  int<lower=1> age_E[N_E];             // female age (years) for fecundity observations
  vector<lower=0>[N_E] E_obs;          // observed fecundity
  // spawner-smolt productivity
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  int<lower=0> N_X_M;                  // number of spawner-smolt productivity covariates
  matrix[max(year),N_X_M] X_M;         // spawner-smolt productivity covariates
  int<lower=1> smolt_age;              // smolt age
  // smolt abundance and observation error
  int<lower=1,upper=N> N_M_obs;        // number of cases with non-missing smolt abundance obs 
  int<lower=1,upper=N> which_M_obs[N_M_obs]; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            // observed annual smolt abundance (not density)
  int<lower=1,upper=N> N_tau_M_obs;    // number of cases with known smolt observation error SD
  int<lower=1,upper=N> which_tau_M_obs[N_tau_M_obs]; // cases with known smolt observation error SD
  vector<lower=0>[N] tau_M_obs;        // known smolt observation error SDs
  int<lower=0,upper=N> N_upstream;     // number of smolt obs upstream of a trap in another pop
  int<lower=1,upper=N> which_upstream[N_upstream]; // smolt obs upstream of a trap in another pop
  int<lower=1,upper=N> downstream_trap[N_upstream]; // row index for the downstream trap
  // SAR (smolt-spawner survival)
  int<lower=0> N_X_MS;                 // number of SAR productivity covariates
  matrix[max(year),N_X_MS] X_MS;       // SAR covariates
  // fishery and hatchery removals and translocations
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  vector<lower=0>[N] S_add_obs;        // number of translocated spawners added to population
  // spawner abundance and observation error
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  int<lower=1,upper=N> N_tau_S_obs;    // number of cases with known spawner observation error SD
  int<lower=1,upper=N> which_tau_S_obs[N_tau_S_obs]; // cases with known spawner observation error SD
  vector<lower=0>[N] tau_S_obs;        // known spawner observation error SDs
  // spawner age structure and sex ratio
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  int<lower=0> n_M_obs[N];             // observed count of male spawners
  int<lower=0> n_F_obs[N];             // observed count of female spawners
  vector<lower=0,upper=1>[N] p_G_obs;  // proportion green (fully fecund) females
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // observed count of wild spawners
  int<lower=0> n_H_obs[N_H];           // observed count of hatchery spawners
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=0> ocean_ages[N_age];          // ocean ages
  vector<lower=1,upper=1>[N_age] ones_N_age = rep_vector(1,N_age); // for rowsums of p matrix 
  vector<lower=1,upper=1>[N] ones_N = rep_vector(1,N); // for elementwise inverse of rowsums 
  int<lower=1> max_ocean_age = max_age - smolt_age; // maximum ocean age
  int<lower=1> min_ocean_age = max_ocean_age - N_age + 1; // minimum ocean age
  int<lower=1> min_age = max_age - N_age + 1; // minimum adult age
  int<lower=1> a_E[N_E];                   // female age index for fecundity observations
  int<lower=1> pop_year_indx[N];           // index of years within each pop, starting at 1
  int<lower=0> n_MF_obs[N];                // total sample sizes for sex frequencies
  vector<lower=0,upper=1>[N] p_NG_obs = 1 - p_G_obs; // proportion non-green females
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  real mu_mu_Mmax = quantile(log(M_obs[which_M_obs]/1000), 0.9); // prior mean of mu_Mmax (thousands of smolts)
  real mu_M_init = mean(log(M_obs[which_M_obs])); // prior log-mean of smolt abundance in years 1:smolt_age
  real sigma_M_init = 2*sd(log(M_obs[which_M_obs])); // prior log-SD of smolt abundance in years 1:smolt_age
  vector[max_ocean_age*N_pop] mu_S_init;   // prior mean of total spawner abundance in years 1:max_ocean_age
  real sigma_S_init = 2*sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_ocean_age
  matrix[N_age,max_ocean_age*N_pop] mu_q_init; // prior counts of wild spawner age distns in years 1:max_ocean_age
  
  for(a in 1:N_age) ocean_ages[a] = min_ocean_age - 1 + a;
  for(i in 1:N) n_MF_obs[i] = n_M_obs[i] + n_F_obs[i];
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  for(i in 1:N_E) a_E[i] = age_E[i] - min_age + 1;
  
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  
  for(i in 1:max_ocean_age)
  {
    int N_orphan_age = N_age - max(i - min_ocean_age, 0); // number of orphan age classes
    int N_amalg_age = N_age - N_orphan_age + 1;           // number of amalgamated age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_ocean_age + i; // index into S_init, q_init
      
      // S_init prior mean that scales observed log-mean by fraction of orphan age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_age) - log(N_age);
      
      // prior on q_init that implies q_orphan ~ Dir(1)
      mu_q_init[,ii] = append_row(rep_vector(1.0/N_amalg_age, N_amalg_age), 
                                  rep_vector(1, N_orphan_age - 1));
    }
  }
}

parameters {
  // fecundity
  vector<lower=0>[N_age] mu_E;           // mean fecundity at each female age
  vector<lower=0>[N_age] sigma_E;        // among-female SD of fecundity at each age
  real<lower=0,upper=1> delta_NG;        // fecundity of non-green females as proportion of green
  // spawner-smolt productivity
  real<lower=0,upper=1> mu_psi;          // hyper-mean density-independent egg-smolt survival
  real<lower=0> sigma_psi;               // hyper-SD logit density-independent egg-smolt survival
  vector[N_pop] zeta_psi;                // logit density-independent egg-smolt survival (Z-scores)
  real mu_Mmax;                          // hyper-mean log maximum smolt recruitment (thousands)
  real<lower=0> sigma_Mmax;              // among-pop hyper-SD of log maximum smolt recruitment
  vector[N_pop] zeta_Mmax;               // log maximum smolt recruitment (Z-scores)
  vector[N_X_M] beta_M;                  // regression coefs for annual spawner-smolt productivity anomalies
  real<lower=-1,upper=1> rho_M;          // AR(1) coef of annual spawner-smolt productivity anomalies
  real<lower=0> sigma_year_M;            // innovation SD of annual spawner-smolt productivity anomalies
  vector[N_year] zeta_year_M;            // annual spawner-smolt productivity anomalies (Z-scores)
  real<lower=0> sigma_M;                 // SD of unique spawner-smolt productivity process errors
  vector[N] zeta_M;                      // unique spawner-smolt productivity process errors (Z-scores)
  // SAR
  real<lower=0,upper=1> mu_MS;           // hyper-mean SAR
  vector[N_X_MS] beta_MS;                // regression coefs for annual logit SAR anomalies
  real<lower=-1,upper=1> rho_MS;         // AR(1) coef of annual logit SAR anomalies
  real<lower=0> sigma_year_MS;           // innovation SD of annual logit SAR anomalies
  vector[N_year] zeta_year_MS;           // annual logit SAR anomalies (Z-scores)
  real<lower=0> sigma_MS;                // SD of unique logit SAR process errors
  vector[N] zeta_MS;                     // unique logit SAR process errors (Z-scores)
  // spawner age structure
  simplex[N_age] mu_p;                   // among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_pop_p;  // among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_pop_p; // Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] zeta_pop_p;      // population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_p;      // SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_p;     // Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;              // log-ratio cohort age distributions (Z-scores)
  // spawner sex ratio
  real<lower=0,upper=1> mu_F;            // hyper-mean proportion female spawners 
  real<lower=0> sigma_pop_F;             // among-pop SD of mean logit proportion female
  vector[N_pop] zeta_pop_F;              // population mean logit proportion female (Z-scores)
  real<lower=0> sigma_F;                 // annual SD of logit proportion female
  vector[N] zeta_F;                      // logit proportion females by outmigration year (Z-scores)
  // H/W composition and removals
  vector<lower=0,upper=1>[N_H] p_HOS;    // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;   // true broodstock take rate when B_take > 0
  // initial states and observation error
  vector<lower=0>[smolt_age*N_pop] M_init; // true smolt abundance in years 1:smolt_age
  vector<lower=0>[max_ocean_age*N_pop] S_init; // true total spawner abundance in years 1:max_ocean_age
  simplex[N_age] q_init[max_ocean_age*N_pop];  // true wild spawner age distributions in years 1:max_ocean_age
  vector<lower=0,upper=1>[max_ocean_age*N_pop] q_F_init; // true proportion females in years 1:max_ocean_age
  real<lower=0> mu_tau_M;                // median smolt observation error SD
  real<lower=0> sigma_tau_M;             // log-SD of smolt observation error SDs
  vector[N] tau_M_z;                     // log of smolt observation error SDs (Z-scores)
  real<lower=0> mu_tau_S;                // median spawner observation error SD
  real<lower=0> sigma_tau_S;             // log-SD of spawner observation error SDs
  vector[N] tau_S_z;                     // log of spawner observation error SDs (Z-scores)
}

transformed parameters {
  // fecundity
  vector<lower=0>[N] E_hat;              // expected egg deposition
  // spawner-smolt productivity
  vector<lower=0,upper=1>[N_pop] psi;    // pop-specific density-independent egg-smolt survival 
  vector<lower=0>[N_pop] Mmax;           // pop-specific maximum smolt recruitment (thousands)
  vector[N_year] eta_year_M;             // annual spawner-smolt productivity anomalies
  vector<lower=0>[N] M_hat;              // expected smolt abundance by brood year
  vector<lower=0>[N] M0;                 // true smolt abundance by brood year
  vector<lower=0>[N] M;                  // true smolt abundance by outmigration year
  // SAR
  vector[N_year] eta_year_MS;            // annual logit SAR anomalies
  vector<lower=0,upper=1>[N] s_MS;       // true SAR
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  row_vector[N_age-1] mu_alr_p;          // mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] mu_pop_alr_p;    // population mean log-ratio age distributions
  matrix[N,N_age] alr_p;                 // alr(p[i,])
  matrix<lower=0>[N,N_age] exp_alr_p;    // exp(alr(p[i,]))
  matrix<lower=0,upper=1>[N,N_age] p;    // true adult age distributions by outmigration year
  matrix<lower=0,upper=1>[N,N_age] q;    // true spawner age distributions
  // spawner sex ratio
  vector[N] p_F;                         // proportion females by outmigration year
  vector[N] q_F;                         // proportion females by return year
  // observation error SDs
  vector<lower=0>[N] tau_M;              // smolt observation error SDs
  vector<lower=0>[N] tau_S;              // spawner observation error SDs

  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;

  // Population-specific egg-smolt survival and capacity terms (independent)
  psi = inv_logit(logit(mu_psi) + sigma_psi*zeta_psi);
  Mmax = exp(mu_Mmax + sigma_Mmax*zeta_Mmax);
  
  // AR(1) model for spawner-smolt productivity and SAR anomalies
  eta_year_M[1] = zeta_year_M[1]*sigma_year_M/sqrt(1 - rho_M^2);     // initial anomaly
  eta_year_MS[1] = zeta_year_MS[1]*sigma_year_MS/sqrt(1 - rho_MS^2); // initial anomaly
  for(i in 2:N_year)
  {
    eta_year_M[i] = rho_M*eta_year_M[i-1] + zeta_year_M[i]*sigma_year_M;
    eta_year_MS[i] = rho_MS*eta_year_MS[i-1] + zeta_year_MS[i]*sigma_year_MS;
  }
  // constrain "fitted" anomalies to sum to 0 (columns of X should be centered)
  eta_year_M = eta_year_M - mean(eta_year_M[1:N_year]) + mat_lmult(X_M, beta_M);
  eta_year_MS = eta_year_MS - mean(eta_year_MS[1:N_year]) + mat_lmult(X_MS, beta_MS);
  // annual population-specific SAR
  s_MS = inv_logit(logit(mu_MS) + eta_year_MS[year] + sigma_MS*zeta_MS);
  
  // Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_alr_p = to_row_vector(log(mu_p[1:(N_age-1)]) - log(mu_p[N_age]));
  mu_pop_alr_p = rep_matrix(mu_alr_p,N_pop) + diag_pre_multiply(sigma_pop_p, L_pop_p * zeta_pop_p')';
  alr_p = append_col(mu_pop_alr_p[pop,] + diag_pre_multiply(sigma_p, L_p * zeta_p')', rep_vector(0,N));
  // Inverse log-ratio (softmax) transform of cohort age distn
  exp_alr_p = exp(alr_p);
  p = diag_pre_multiply(ones_N ./ (exp_alr_p * ones_N_age), exp_alr_p);
  
  // Annual population-specific proportion females
  p_F = inv_logit(logit(mu_F) + sigma_pop_F*zeta_pop_F[pop] + sigma_F*zeta_F);
  
  // Calculate true total wild and hatchery spawners, spawner age distribution, 
  // eggs and smolts, and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] S_W_a;       // wild spawners by age
    int ii;                        // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_ocean_age, 0), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution
    vector[N_age] q_F_a;           // proportion of each age that are female
    
    // Smolt recruitment
    if(pop_year_indx[i] <= smolt_age)
      M[i] = M_init[(pop[i]-1)*smolt_age + pop_year_indx[i]];  // use initial values
    else
      M[i] = M0[i-smolt_age];  // smolts from appropriate brood year
    
    // Spawners and age structure
    // Use initial values for orphan age classes, otherwise use process model
    if(pop_year_indx[i] <= max_ocean_age)
    {
      ii = (pop[i] - 1)*max_ocean_age + pop_year_indx[i];
      q_orphan = append_row(sum(head(q_init[ii], N_age - N_orphan_age + 1)), 
                                tail(q_init[ii], N_orphan_age - 1));
    }
    
    for(a in 1:N_age)
    {
      if(ocean_ages[a] < pop_year_indx[i])
      {
        // Use recruitment process model
        S_W_a[a] = M[i-ocean_ages[a]]*s_MS[i-ocean_ages[a]]*p[i-ocean_ages[a],a];
        q_F_a[a] = p_F[i-ocean_ages[a]];
      }
      else 
      {
        // Use initial values
        S_W_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
        q_F_a[a] = q_F_init[ii];
      }
    }
    
    // catch and broodstock removal and translocations
    S_W_a = S_W_a*(1 - F_rate[i])*(1 - B_rate_all[i]);
    S_W[i] = sum(S_W_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i] + S_add_obs[i];
    q[i,] = S_W_a/S_W[i];
    q_F[i] = q[i,]*q_F_a;
    
    // Density-independent egg production from brood year i
    // weighted by age structure and sex ratio 
    // discounted for proportion of non-green (not fully fecund) females
    if(is_nan(q[i,]*mu_E*q_F[i]*S[i]))
      print("i = ", i, "  q[i,] = ", q[i,], "  q_F[i] = ", q_F[i], "  S[i] = ", S[i], "  E_hat[i] = ", q[i,]*mu_E*p_F[i]*S[i]);
    E_hat[i] = q[i,] * mu_E * q_F[i] * (p_G_obs[i] + delta_NG * p_NG_obs[i]) * S[i];

    // Smolt production from brood year i
    // Density-dependent egg-to-smolt survival
    M_hat[i] = SR(SR_fun, psi[pop[i]], Mmax[pop[i]]*1e3, E_hat[i], A[i]);
    M0[i] = M_hat[i] * exp(eta_year_M[year[i]] + sigma_M*zeta_M[i]);
  }
  
  // Impute missing observation error SDs from lognormal hyperdistribution
  // (use Matt trick)
  tau_M = exp(log(mu_tau_M) + sigma_tau_M*tau_M_z);
  tau_M[which_tau_M_obs] = tau_M_obs[which_tau_M_obs];
  tau_S = exp(log(mu_tau_S) + sigma_tau_S*tau_S_z);
  tau_S[which_tau_S_obs] = tau_S_obs[which_tau_S_obs];
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  vector[N] M_downstream; // total smolts including upstream populations 
  
  // Priors
  
  // fecundity and sex ratio
  mu_E ~ normal(2500,500);
  sigma_E ~ normal(500,1000);

  // spawner-smolt productivity
  sigma_psi ~ normal(0,1);        // mildly regularize long right tail
  mu_Mmax ~ normal(mu_mu_Mmax,3); // units of Mmax: thousands of smolts
  sigma_Mmax ~ normal(0,3);
  zeta_psi ~ std_normal();        // logit(psi) ~ N(logit(mu_psi), sigma_psi)
  zeta_Mmax ~ std_normal();       // log(Mmax) ~ N(mu_Mmax, sigma_Mmax)
  beta_M ~ normal(0,2);
  rho_M ~ pexp(0,0.85,20);        // regularize away from 1 to ensure stationarity
  sigma_year_M ~ normal(0,2);
  zeta_year_M ~ std_normal();     // eta_year_M[i] ~ N(rho_M*phi_M[i-1], sigma_year_M)
  sigma_M ~ normal(0,2);
  zeta_M ~ std_normal();          // M0 ~ lognormal(log(M_hat) + eta_year_M, sigma_M)
  
  // SAR
  beta_MS ~ normal(0,2);
  rho_MS ~ pexp(0,0.85,20);       // regularize away from 1 to ensure stationarity
  sigma_year_MS ~ normal(0,2);
  zeta_year_MS ~ std_normal();    // phi_MS[i] ~ N(rho_MS*phi_MS[i-1], sigma_year_MS)
  sigma_MS ~ normal(0,2);
  zeta_MS ~ std_normal();         // SAR: logit(s_MS) ~ N(logit(s_MS_hat), sigma_MS)

  // spawner age structure
  to_vector(sigma_pop_p) ~ normal(0,2);
  to_vector(sigma_p) ~ normal(0,2);
  L_pop_p ~ lkj_corr_cholesky(1);
  L_p ~ lkj_corr_cholesky(1);
  // pop mean age probs logistic MVN: 
  // mu_pop_alr_p[i,] ~ MVN(mu_alr_p, D*R_pop_p*D), where D = diag_matrix(sigma_pop_p)
  to_vector(zeta_pop_p) ~ std_normal();
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_pop_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();
  
  // spawner sex ratio
  sigma_pop_F ~ normal(0,2);
  zeta_pop_F ~ std_normal();   // eta_pop_F ~ N(0, sigma_pop_F)
  sigma_F ~ normal(0,2);
  zeta_F ~ std_normal();       // logit(p_F[i]) ~ N(eta_pop_F[pop[i]], sigma_F)
  q_F_init ~ beta(3,3);        // mildly regularize initial states toward 0.5

  // removals
  log_B_take = log(S_W[which_B]) + logit(B_rate); // B_take = S_W*B_rate/(1 - B_rate)
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match 

  // initial states
  // (accounting for amalgamation of q_init to q_orphan)
  M_init ~ lognormal(mu_M_init, sigma_M_init);
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_ocean_age*N_pop] q_init_mat;
    
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // observation error SDs
  // hyperparameters
  mu_tau_M ~ normal(0,1); // half-normal
  sigma_tau_M ~ normal(0,5);
  mu_tau_S ~ normal(0,1); // half-normal
  sigma_tau_S ~ normal(0,5);
  // unknown SDs
  tau_M_z ~ std_normal(); // tau_M[-which_tau_M_obs] ~ lognormal(log(mu_tau_M), sigma_tau_M);
  tau_S_z ~ std_normal(); // tau_S[-which_tau_S_obs] ~ lognormal(log(mu_tau_S), sigma_tau_S);
  // known SDs
  tau_M_obs[which_tau_M_obs] ~ lognormal(log(mu_tau_M), sigma_tau_M);
  tau_S_obs[which_tau_S_obs] ~ lognormal(log(mu_tau_S), sigma_tau_S);

  // Observation model
  E_obs ~ normal(mu_E[a_E], sigma_E[a_E]);              // observed fecundity
  target += -normal_lccdf(0 | mu_E[a_E], sigma_E[a_E]); // zero-truncated normal (T[0, ] not vectorized)
  M_downstream = M;
  M_downstream[downstream_trap] += M[which_upstream];
  M_downstream[which_M_obs] ~ lognormal(log(M_obs[which_M_obs]), tau_M[which_M_obs]); // prior on smolts
  S[which_S_obs] ~ lognormal(log(S_obs[which_S_obs]), tau_S[which_S_obs]); // prior on spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
  n_F_obs ~ binomial(n_MF_obs, q_F);   // observed counts of female vs. male spawners
}

generated quantities {
  corr_matrix[N_age-1] R_pop_p; // among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_p;     // correlation matrix of within-pop cohort log-ratio age distns 
  // ?? vector[N_E] LL_E_obs; ??
  vector[N] M_downstream;       // total smolts including upstream populations 
  vector[N] LL_M_obs;           // pointwise log-likelihood of smolts
  vector[N] LL_S_obs;           // pointwise log-likelihood of spawners
  vector[N] LL_n_F_obs;         // pointwise log-likelihood of female vs. male frequencies
  vector[N_H] LL_n_H_obs;       // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;       // pointwise log-likelihood of wild age frequencies
  vector[N] LL;                 // total pointwise log-likelihood                              
  
  R_pop_p = multiply_lower_tri_self_transpose(L_pop_p);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  M_downstream = M;
  M_downstream[downstream_trap] += M[which_upstream];
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M[which_M_obs[i]] | log(M_obs[which_M_obs[i]]), tau_M[which_M_obs[i]]); 
    // LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M[which_M_obs[i]]); 
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S[which_S_obs[i]] | log(S_obs[which_S_obs[i]]), tau_S[which_S_obs[i]]); 
    // LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S[which_S_obs[i]]); 
  for(i in 1:N)
    LL_n_F_obs[i] = binomial_lpmf(n_F_obs[i] | n_MF_obs[i], q_F[i]);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL = LL_M_obs + LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
