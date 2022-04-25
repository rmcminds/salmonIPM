functions {
#include /include/SR.stan
#include /include/pexp_lpdf_vec.stan
#include /include/mat_lmult.stan
#include /include/quantile.stan
#include /include/posdiff.stan
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // calendar year identifier
  // smolt production
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> smolt_age;              // smolt age
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  int<lower=0> K_alpha;                // number of intrinsic productivity covariates
  matrix[N,K_alpha] X_alpha;           // intrinsic productivity covariates
  int<lower=0> K_Mmax;                 // number of maximum smolt recruitment covariates
  matrix[N,K_Mmax] X_Mmax;             // maximum smolt recruitment covariates
  int<lower=0> K_M;                    // number of smolt recruitment covariates
  row_vector[K_M] X_M[N];              // smolt recruitment covariates
  // smolt abundance
  int<lower=1,upper=N> N_M_obs;        // number of cases with non-missing smolt abundance obs 
  int<lower=1,upper=N> which_M_obs[N_M_obs]; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            // observed annual smolt abundance (not density)
  // SAR (sMolt-Spawner survival)
  int<lower=0> K_MS;                   // number of SAR covariates
  row_vector[K_MS] X_MS[N];            // SAR covariates
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  // spawner age structure
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=0> ocean_ages[N_age];          // ocean ages
  int<lower=1> max_ocean_age = max_age - smolt_age; // maximum ocean age
  int<lower=1> min_ocean_age = max_ocean_age - N_age + 1; // minimum ocean age
  int<lower=1> pop_year_indx[N];           // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  vector<lower=0>[N] B_take_all;           // broodstock take of wild adults in all cases
  int<lower=0,upper=1> use_B[N];           // binary indicator of B_take_obs > 0
  real mu_Mmax = quantile(log(M_obs[which_M_obs]), 0.9); // prior log-mean of Mmax
  real sigma_Mmax = sd(log(M_obs[which_M_obs])); // prior log-SD of Mmax
  real mu_M_init = mean(log(M_obs[which_M_obs])); // prior log-mean of smolt abundance in years 1:smolt_age
  real sigma_M_init = sd(log(M_obs[which_M_obs])); // prior log-SD of smolt abundance in years 1:smolt_age
  real sigma_S_init = sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_ocean_age
  vector[max_ocean_age*N_pop] mu_S_init;   // prior mean of total spawner abundance in years 1:max_ocean_age
  matrix[N_age,max_ocean_age*N_pop] mu_q_init; // prior counts of wild spawner age distns in years 1:max_ocean_age
  
  for(a in 1:N_age) ocean_ages[a] = min_ocean_age - 1 + a;
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  
  B_take_all = rep_vector(0,N);
  B_take_all[which_B] = B_take_obs;
  for(i in 1:N) use_B[i] = B_take_all[i] > 0;

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
  // smolt recruitment
  vector<lower=0>[N_pop] alpha;             // intrinsic spawner-smolt productivity
  matrix[N_pop,K_alpha] beta_alpha;         // regression coefs for log alpha
  vector<lower=0>[N_pop] Mmax;              // asymptotic smolt recruitment
  matrix[N_pop,K_Mmax] beta_Mmax;           // regression coefs for log Mmax
  matrix[N_pop,K_M] beta_M;                 // regression coefs for smolt recruitment
  vector<lower=-1,upper=1>[N_pop] rho_M;    // AR(1) coefs for spawner-smolt productivity
  vector<lower=0>[N_pop] sigma_M;           // spawner-smolt process error SDs
  vector[N] zeta_M;                         // smolt recruitment process errors (z-scored)
  // SAR
  vector<lower=0,upper=1>[N_pop] mu_MS;     // mean SAR
  matrix[N_pop,K_MS] beta_MS;               // regression coefs for SAR
  vector<lower=-1,upper=1>[N_pop] rho_MS;   // AR(1) coefs for SAR
  vector<lower=0>[N_pop] sigma_MS;          // SAR process error SDs
  vector[N] zeta_MS;                        // SAR process errors (z-scored)
  // spawner age structure
  simplex[N_age] mu_p[N_pop];               // population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p;   // log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop]; // Cholesky factors of correlation matrices of cohort log-ratio age distns
  matrix[N,N_age-1] zeta_p;                 // log-ratio cohort age distribution errors (Z-scores)
  // H/W composition
  vector<lower=0,upper=1>[N_H] p_HOS;       // true p_HOS in years which_H
  // initial states, observation error
  vector<lower=0>[smolt_age*N_pop] M_init;  // true smolt abundance in years 1:smolt_age
  vector<lower=0>[max_ocean_age*N_pop] S_init; // true total spawner abundance in years 1:max_ocean_age
  simplex[N_age] q_init[max_ocean_age*N_pop];  // true wild spawner age distributions in years 1:max_ocean_age
  vector<lower=0>[N_pop] tau_M;             // smolt observation error SDs
  vector<lower=0>[N_pop] tau_S;             // spawner observation error SDs
}

transformed parameters {
  // smolt recruitment
  vector<lower=0>[N] alpha_Xbeta;      // intrinsic productivity including covariate effects
  vector<lower=0>[N] Mmax_Xbeta;       // maximum recruitment including covariate effects
  vector<lower=0>[N] M_hat;            // expected smolt abundance (not density) by brood year
  vector[N] epsilon_M;                 // process error in smolt abundance by brood year 
  vector<lower=0>[N] M0;               // true smolt abundance (not density) by brood year
  vector<lower=0>[N] M;                // true smolt abundance (not density) by outmigration year
  // SAR
  vector[N] epsilon_MS;                // process error in SAR by outmigration year 
  vector<lower=0>[N] s_MS;             // true SAR by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                 // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;              // true total wild spawner abundance
  vector[N] S_H;                       // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                // true total spawner abundance
  vector<lower=0>[N] B_take_adj;       // adjusted broodstock take
  // spawner age structure
  vector[N_age-1] mu_alr_p[N_pop];     // population mean log ratio age distributions
  simplex[N_age] p[N];                 // true adult age distributions by outmigration year
  matrix<lower=0,upper=1>[N,N_age] q;  // true spawner age distributions
  
  // Pad p_HOS and B_take
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_take_adj = rep_vector(0,N);
    
  // S-R parameters including covariate effects
  alpha_Xbeta = alpha[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
  Mmax_Xbeta = Mmax[pop] .* exp(rows_dot_product(X_Mmax, beta_Mmax[pop,]));

  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    mu_alr_p[j] = log(head(mu_p[j], N_age-1)) - log(mu_p[j,N_age]);
  
  // Calculate true total wild and hatchery spawners, spawner age distribution, and smolts,
  // and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    int ii;                  // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_ocean_age, 0), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution (amalgamated simplex)
    vector[N_age] alr_p;     // alr(p[i,])
    row_vector[N_age] S_W_a; // true wild spawners by age
    
    // Within-pop, time-varying IID age vectors
    // (multivariate Matt trick)
    alr_p = rep_vector(0,N_age);
    alr_p[1:(N_age-1)] = mu_alr_p[pop[i]] + sigma_p[pop[i],]' .* (L_p[pop[i]] * zeta_p[i,]');
    alr_p = exp(alr_p);
    p[i] = alr_p / sum(alr_p);
    
    // AR(1) smolt recruitment and SAR process errors  
    if(pop_year_indx[i] == 1) // initial process error
    {
      epsilon_M[i] = zeta_M[i] * sigma_M[pop[i]] / sqrt(1 - rho_M[pop[i]]^2);
      epsilon_MS[i] = zeta_MS[i] * sigma_MS[pop[i]] / sqrt(1 - rho_MS[pop[i]]^2);
    }
    else
    {
      epsilon_M[i] = rho_M[pop[i]] * epsilon_M[i-1] + zeta_M[i] * sigma_M[pop[i]];
      epsilon_MS[i] = rho_MS[pop[i]] * epsilon_MS[i-1] + zeta_MS[i] * sigma_MS[pop[i]];
    }
    // SAR for outmigration year i
    s_MS[i] = inv_logit(logit(mu_MS[pop[i]]) + dot_product(X_MS[i], beta_MS[pop[i],]) + epsilon_MS[i]);
    
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
        // Use recruitment process model
        S_W_a[a] = M[i-ocean_ages[a]]*s_MS[i-ocean_ages[a]]*p[i-ocean_ages[a],a];
      else
        // Use initial values
        S_W_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
    }
    
    // catch and broodstock removal (assumes no take of age 1)
    S_W_a[2:N_age] = S_W_a[2:N_age]*(1 - F_rate[i]);
    if(use_B[i])
    {
      real S_B = sum(S_W_a[2:N_age]);
      B_take_adj[i] = posdiff(S_B, B_take_all[i], 0.01);
      S_W_a[2:N_age] = S_W_a[2:N_age]*(1 - B_take_adj[i]/S_B);
    }
    S_W[i] = sum(S_W_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];
    q[i,] = S_W_a / S_W[i];

    // Smolt production from brood year i
    M_hat[i] = SR(SR_fun, alpha_Xbeta[i], Mmax_Xbeta[i], S[i], A[i]);
    M0[i] = M_hat[i] * exp(dot_product(X_M[i], beta_M[pop[i],]) + epsilon_M[i]);
  }
}

model {
  // Priors
  
  // smolt recruitment
  alpha ~ lognormal(2.0,2.0);
  Mmax ~ lognormal(mu_Mmax, sigma_Mmax);
  to_vector(beta_M) ~ normal(0,5);
  rho_M ~ pexp(0,0.85,20);  // mildly regularize rho to ensure stationarity
  sigma_M ~ normal(0,3);
  zeta_M ~ std_normal();    // total smolts: log(M) ~ normal(log(M_hat), sigma_M)

  // SAR
  to_vector(beta_MS) ~ normal(0,3);
  rho_MS ~ pexp(0,0.85,20); // mildly regularize rho to ensure stationarity
  sigma_MS ~ normal(0,3);
  zeta_MS ~ std_normal();   // SAR: logit(s_MS) ~ normal(logit(s_MS_hat), sigma_MS)

  // spawner age structure
  to_vector(sigma_p) ~ normal(0,5);
  for(j in 1:N_pop)
    L_p[j] ~ lkj_corr_cholesky(3);
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();
  
  // removals
  // penalty to force obs and adj broodstock take to match
  B_take_obs ~ lognormal(log(B_take_adj[which_B]), 0.01); 

  // initial states
  // (accounting for amalgamation of q_init to q_orphan)
  M_init ~ lognormal(mu_M_init, sigma_M_init);
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_ocean_age*N_pop] q_init_mat;
    
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // observation error
  tau_M ~ pexp(1,0.85,30);   // rule out tau < 0.1 to avoid divergences 
  tau_S ~ pexp(1,0.85,30);   // rule out tau < 0.1 to avoid divergences 

  // Observation model
  M_obs[which_M_obs] ~ lognormal(log(M[which_M_obs]), tau_M[pop[which_M_obs]]);  // observed smolts
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S[pop[which_S_obs]]);  // observed spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_p[N_pop]; // correlation matrices of within-pop cohort log-ratio age distns
  vector[N] LL_M_obs;              // pointwise log-likelihood of smolts
  vector[N] LL_S_obs;              // pointwise log-likelihood of spawners
  vector[N_H] LL_n_H_obs;          // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;          // pointwise log-likelihood of wild age frequencies
  vector[N] LL;                    // total pointwise log-likelihood                              
  
  for(j in 1:N_pop)
    R_p[j] = multiply_lower_tri_self_transpose(L_p[j]);
  
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M[pop[which_M_obs[i]]]); 
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S[pop[which_M_obs[i]]]); 
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_M_obs + LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
