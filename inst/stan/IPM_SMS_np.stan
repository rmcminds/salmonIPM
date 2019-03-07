functions {
  # spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      # discrete exponential
      R = alpha*S/A;
    else if(SR_fun == 2) # Beverton-Holt
      R = alpha*S/(A + alpha*S/Rmax);
    else if(SR_fun == 3) # Ricker
      R = alpha*(S/A)*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  # Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
}

data {
  int<lower=1> SR_fun;                 # S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;                      # total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         # population identifier
  int<lower=1,upper=N> year[N];        # calendar year identifier
  int<lower=1> N_X_M;                  # number of spawner-smolt productivity covariates
  matrix[max(year),N_X_M] X_M;         # spawner-smolt covariates (if none, use vector of zeros)
  int<lower=1> N_X_MS;                 # number of SAR productivity covariates
  matrix[max(year),N_X_MS] X_MS;       # SAR covariates (if none, use vector of zeros)
  int<lower=0,upper=max(pop)> N_pop_H; # number of populations with hatchery input
  int<lower=1,upper=max(pop)> which_pop_H[max(N_pop_H,1)]; # populations with hatchery input
  int<lower=1,upper=N> N_S_obs;        # number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; # cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            # observed annual total spawner abundance (not density)
  int<lower=1,upper=N> N_M_obs;        # number of cases with non-missing smolt abundance obs 
  int<lower=1,upper=N> which_M_obs[N_M_obs]; # cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            # observed annual smolt abundance (not density)
  int<lower=1> smolt_age;              # smolt age
  int<lower=2> N_age;                  # number of adult age classes
  int<lower=2> max_age;                # maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  # observed wild spawner age frequencies (all zero row = NA)  
  int<lower=0,upper=N> N_H;            # number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[max(N_H,1)]; # years with p_HOS > 0
  int<lower=0> n_W_obs[max(N_H,1)];    # count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[max(N_H,1)];    # count of hatchery spawners in samples (assumes no NAs)
  vector[N] A;                         # habitat area associated with each spawner abundance obs
  vector[N] F_rate;                    # fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            # number of years with B_take > 0
  int<lower=1,upper=N> which_B[max(N_B,1)]; # years with B_take > 0
  vector[max(N_B,1)] B_take_obs;       # observed broodstock take of wild adults
}

transformed data {
  int<lower=1,upper=N> N_pop;         # number of populations
  int<lower=1,upper=N> N_year;        # number of years
  int<lower=1> ocean_ages[N_age];     # ocean ages
  int<lower=1> pop_year_indx[N];      # index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[max(N_H,1)];  # total sample sizes for H/W frequencies
  
  N_pop = max(pop);
  N_year = max(year);
  for(a in 1:N_age)
    ocean_ages[a] = max_age - smolt_age - N_age + a;
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  for(i in 1:max(N_H,1)) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
}

parameters {
  vector<lower=0>[N_pop] alpha;               # intrinsic spawner-smolt productivity
  vector<lower=0>[N_pop] Rmax;                # asymptotic smolt recruitment
  matrix[N_pop,N_X_M] beta_M;                 # regression coefs for spawner-smolt productivity 
  vector<lower=-1,upper=1>[N_pop] rho_M;      # AR(1) coefs for spawner-smolt productivity
  vector<lower=0>[N_pop] sigma_M;             # spawner-smolt process error SDs
  vector[N] zeta_M;                           # smolt recruitment process errors (z-scored)
  vector<lower=0,upper=1>[N_pop] mu_MS;       # mean SAR
  matrix[N_pop,N_X_MS] beta_MS;               # regression coefs for SAR 
  vector<lower=-1,upper=1>[N_pop] rho_MS;     # AR(1) coefs for SAR
  vector<lower=0>[N_pop] sigma_MS;            # SAR process error SDs
  vector[N] zeta_MS;                          # SAR process errors (z-scored)
  simplex[N_age] mu_p[N_pop];                 # population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p;     # log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop];   # Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;                   # log-ratio cohort age distribution errors (Z-scores)
  vector<lower=0>[smolt_age*N_pop] M_init;    # true smolt abundance in years 1:smolt_age
  vector<lower=0>[max_age*N_pop] S_init;      # true total spawner abundance in years 1:max_age
  simplex[N_age] q_init[max_age*N_pop];       # true wild spawner age distributions in years 1:max_age
  vector<lower=0,upper=1>[max(N_H,1)] p_HOS;  # true p_HOS in years which_H
  vector<lower=0,upper=1>[max(N_B,1)] B_rate; # true broodstock take rate when B_take > 0
  vector<lower=0>[N_pop] tau_M;               # smolt observation error SDs
  vector<lower=0>[N_pop] tau_S;               # spawner observation error SDs
}

transformed parameters {
  vector<lower=0>[N] S_W;                # true total wild spawner abundance
  vector[N] S_H;                         # true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  # true total spawner abundance
  matrix[N_pop,N_age-1] gamma;           # population mean log ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    # true adult age distributions by outmigration year
  matrix<lower=0,upper=1>[N,N_age] q;    # true spawner age distributions
  vector[N] p_HOS_all;                   # true p_HOS in all years (can == 0)
  vector<lower=0>[N] M_hat;              # expected smolt abundance (not density) by brood year
  vector[N] epsilon_M;                   # process error in smolt abundance by brood year 
  vector<lower=0>[N] M0;                 # true smolt abundance (not density) by brood year
  vector<lower=0>[N] M;                  # true smolt abundance (not density) by outmigration year
  vector[N] epsilon_MS;                  # process error in SAR by outmigration year 
  vector<lower=0>[N] s_MS;               # true SAR by outmigration year
  vector<lower=0,upper=1>[N] B_rate_all; # true broodstock take rate in all years
  
  # Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  if(N_H > 0)
    p_HOS_all[which_H] = p_HOS;
  
  B_rate_all = rep_vector(0,N);
  if(N_B > 0)
    B_rate_all[which_B] = B_rate;
  
  # Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    gamma[j,] = to_row_vector(log(mu_p[j,1:(N_age-1)]) - log(mu_p[j,N_age]));
  
  # Calculate true total wild and hatchery spawners, spawner age distribution, and smolts,
  # and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] alr_p; # temp variable: alr(p[i,])
    row_vector[N_age] S_W_a; # temp variable: true wild spawners by age
    
    # Within-pop, time-varying IID age vectors
    # (multivariate Matt trick)
    alr_p = rep_row_vector(0,N_age);
    alr_p[1:(N_age-1)] = gamma[pop[i],] + sigma_p[pop[i],] .* (L_p[pop[i]] * zeta_p[i,]')';
    alr_p = exp(alr_p);
    p[i,] = alr_p/sum(alr_p);
    
    # AR(1) smolt recruitment process errors 
    # MAR(1) SAR process errors  
    if(pop_year_indx[i] == 1) # initial process error
    {
      epsilon_M[i] = zeta_M[i]*sigma_M[pop[i]]/sqrt(1 - rho_M[pop[i]]^2);
      epsilon_MS[i] = zeta_MS[i]*sigma_MS[pop[i]]/sqrt(1 - rho_MS[pop[i]]^2);
    }
    else
    {
      epsilon_M[i] = rho_M[pop[i]]*epsilon_M[i-1] + zeta_M[i]*sigma_M[pop[i]];
      epsilon_MS[i] = rho_MS[pop[i]]*epsilon_MS[i-1] + zeta_MS[i]*sigma_MS[pop[i]];
    }
    # SAR for outmigration year i
    s_MS[i] = inv_logit(logit(mu_MS[pop[i]]) + dot_product(X_MS[year[i],], beta_MS[pop[i],]) + epsilon_MS[i]); 
    
    # Smolt recruitment
    if(pop_year_indx[i] <= smolt_age)
      M[i] = M_init[(pop[i]-1)*smolt_age + pop_year_indx[i]];  # use initial values
    else
      M[i] = M0[i-smolt_age];  # smolts from appropriate brood year
    
    # Spawners and age structure
    if(pop_year_indx[i] <= max_age)
    {
      # use initial values
      S_W[i] = S_init[(pop[i]-1)*max_age + pop_year_indx[i]]*(1 - p_HOS_all[i]);        
      S_H[i] = S_init[(pop[i]-1)*max_age + pop_year_indx[i]]*p_HOS_all[i];
      q[i,1:N_age] = to_row_vector(q_init[(pop[i]-1)*max_age + pop_year_indx[i],1:N_age]);
      S_W_a = S_W[i]*q[i,];
    }
    else
    {
      for(a in 1:N_age)
        S_W_a[a] = M[i - ocean_ages[a]]*s_MS[i - ocean_ages[a]]*p[i - ocean_ages[a],a];
      # catch and broodstock removal (assumes no take of age 1)
      S_W_a[2:N_age] = S_W_a[2:N_age]*(1 - F_rate[i])*(1 - B_rate_all[i]);
      S_W[i] = sum(S_W_a);
      S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
      q[i,] = S_W_a/S_W[i];
    }
    
    S[i] = S_W[i] + S_H[i];
    
    # Smolt production
    M_hat[i] = A[i] * SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    # smolts from brood year i
    M0[i] = M_hat[i]*exp(dot_product(X_M[year[i],], beta_M[pop[i],]) + epsilon_M[i]); # smolts from brood year i
  }
}

model {
  vector[max(N_B,1)] B_take; # true broodstock take when B_take_obs > 0
  
  # Priors
  alpha ~ lognormal(2,2);
  Rmax ~ lognormal(2,3);
  to_vector(beta_M) ~ normal(0,5);
  to_vector(beta_MS) ~ normal(0,5);
  for(j in 1:N_pop)
  {
    rho_M[j] ~ pexp(0,0.85,20);   # mildly regularize rho to ensure stationarity
    sigma_M[j] ~ pexp(0,1,10);
    tau_M[j] ~ pexp(1,0.85,30);   # rule out tau < 0.1 to avoid divergences 
    rho_MS[j] ~ pexp(0,0.85,20);  # mildly regularize rho to ensure stationarity
    sigma_MS[j] ~ pexp(0,1,10);
    tau_S[j] ~ pexp(1,0.85,30);   # rule out tau < 0.1 to avoid divergences 
    L_p[j] ~ lkj_corr_cholesky(3);
  }
  to_vector(sigma_p) ~ normal(0,5);
  M_init ~ lognormal(0,5);
  S_init ~ lognormal(0,5);
  if(N_B > 0)
  {
    B_take = B_rate .* S_W[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
    B_take_obs ~ lognormal(log(B_take), 0.1); # penalty to force pred and obs broodstock take to match 
  }
  
  # Hierarchical priors
  # age probs logistic MVN: alr_p[i,] ~ MVN(gamma[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ normal(0,1);
  
  # Process model
  zeta_M ~ normal(0,1);  # total smolts: log(M) ~ normal(log(M_hat), sigma_M)
  zeta_MS ~ normal(0,1); # SAR: logit(s_MS) ~ normal(logit(s_MS_hat), sigma_MS)
  
  # Observation model
  M_obs[which_M_obs] ~ lognormal(log(M[which_M_obs]), tau_M[pop[which_M_obs]]);  # observed smolts
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S[pop[which_S_obs]]);  # observed spawners
  if(N_H > 0) n_H_obs ~ binomial(n_HW_obs, p_HOS);  # observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  # obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_p[N_pop]; # correlation matrices of within-pop cohort log-ratio age distns
  vector[N] LL_M_obs;              # pointwise log-likelihood of smolts
  vector[N] LL_S_obs;              # pointwise log-likelihood of spawners
  vector[max(N_H,1)] LL_n_H_obs;   # pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;          # pointwise log-likelihood of wild age frequencies
  vector[N] LL;                    # total pointwise log-likelihood                              
  
  for(j in 1:N_pop)
    R_p[j] = multiply_lower_tri_self_transpose(L_p[j]);
  
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M); 
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S); 
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_H_obs = rep_vector(0,max(N_H,1));
  if(N_H > 0)
  {
    for(i in 1:N_H)
      LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  }
  LL = LL_M_obs + LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
