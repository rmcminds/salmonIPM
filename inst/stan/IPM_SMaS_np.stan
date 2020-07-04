functions {
  // spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      // discrete exponential
      R = alpha*S/A;
    else if(SR_fun == 2) // Beverton-Holt
      R = alpha*S/(A + alpha*S/Rmax);
    else if(SR_fun == 3) // Ricker
      R = alpha*(S/A)*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  // Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(vector y, real mu, real sigma, real shape) {
    vector[num_elements(y)] LL;
    
    for(i in 1:num_elements(LL))
      LL[i] = -pow(fabs(y[i] - mu)/sigma, shape);
      
    return(sum(LL));
  }
  
  // Column sums of matrix
  row_vector col_sums(matrix X) {
    row_vector[cols(X)] s;
    s = rep_row_vector(1, rows(X)) * X;
    return s;
  }
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // calendar year identifier
  // smolt production
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  int<lower=0> N_X_M;                  // number of spawner-smolt productivity covariates
  matrix[max(year),N_X_M] X_M;         // spawner-smolt covariates (if none, use vector of zeros)
  // smolt abundance
  int<lower=1,upper=N> N_M_obs;        // number of cases with non-missing smolt abundance obs 
  int<lower=1,upper=N> which_M_obs[N_M_obs]; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            // observed annual smolt abundance (not density)
  // smolt age structure
  int<lower=2> N_Mage;                 // number of smolt age classes
  int<lower=2> max_Mage;               // maximum smolt age
  matrix<lower=0>[N,N_Mage] n_Mage_obs; // observed smolt age frequencies (all zero row = NA)  
  // SAR (sMolt-Spawner survival)
  int<lower=0> N_X_MS;                 // number of SAR productivity covariates
  matrix[max(year),N_X_MS] X_MS;       // SAR covariates (if none, use vector of zeros)
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  // spawner ocean age and Gilbert-Rich age structure
  int<lower=2> N_MSage;                // number of ocean age classes
  int<lower=1> max_MSage;              // maximum ocean age
  matrix<lower=0>[N,N_MSage] n_MSage_obs; // observed ocean age frequencies (all zero row = NA)  
  matrix<lower=0>[N,N_Mage*N_MSage] n_GRage_obs;  // obs W spawner Gilbert-Rich age freqs (all zero row = NA)  
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples (assumes no NAs)
}

transformed data {
  int<lower=1,upper=N> N_pop;        // number of populations
  int<lower=1,upper=N> N_year;       // number of years
  int<lower=0> smolt_ages[N_Mage];   // smolt ages
  int<lower=0> ocean_ages[N_MSage];  // ocean ages
  int<lower=2> max_age;              // maximum adult age
  int<lower=2> N_GRage;              // number of Gilbert-Rich age classes
  int<lower=1> pop_year_indx[N];     // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];        // total sample sizes for H/W frequencies
  
  N_pop = max(pop);
  N_year = max(year);
  for(a in 1:N_Mage)
    smolt_ages[a] = max_Mage - N_Mage + a;
  for(a in 1:N_MSage)
    ocean_ages[a] = max_MSage - N_MSage + a;
  max_age = max_Mage + max_MSage;
  N_GRage = N_Mage*N_MSage;
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
}

parameters {
  //?// indicates params that could be arrays instead of matrices
  // smolt recruitment
  vector<lower=0>[N_pop] alpha;               // intrinsic spawner-smolt productivity
  vector<lower=0>[N_pop] Rmax;                // asymptotic smolt recruitment
  matrix[N_pop,N_X_M] beta_M;                 //?// regression coefs for spawner-smolt productivity
  vector<lower=-1,upper=1>[N_pop] rho_M;      // AR(1) coefs for spawner-smolt productivity
  vector<lower=0>[N_pop] sigma_M;             // spawner-smolt process error SDs
  vector[N] zeta_M;                           // smolt recruitment process errors (z-scored)
  // smolt age structure
  simplex[N_Mage] mu_p_M[N_pop];              // population mean smolt age distributions
  matrix<lower=0>[N_pop,N_Mage-1] sigma_p_M;  //?// log-ratio cohort smolt age distribution SDs
  cholesky_factor_corr[N_Mage-1] L_p_M[N_pop]; // Cholesky-factored corr matrices of log-ratio smolt age distns
  matrix[N,N_Mage-1] zeta_p_M;                //?// log-ratio cohort smolt age distn errors (Z-scored)
  // SAR
  matrix<lower=0,upper=1>[N_pop,N_Mage] mu_MS; //?// mean SAR for each smolt age
  matrix[N_pop,N_X_MS] beta_MS;               //?// regression coefs for SAR
  matrix<lower=-1,upper=1>[N_pop,N_Mage] rho_MS; //?// AR(1) coefs of SAR for each smolt age
  matrix<lower=0>[N_pop,N_Mage] sigma_MS;     //?// SAR process error SDs for each smolt age
  cholesky_factor_corr[N_Mage] L_MS[N_pop];   // Cholesky-factored corr matrices of SAR across smolt ages
  matrix[N,N_Mage] zeta_MS;                   //?// SAR process errors for each smolt age (z-scored)
  // ocean age structure
  simplex[N_MSage] mu_p_MS[N_pop,N_Mage];     // pop mean ocean age distributions for each smolt age
  vector<lower=0>[N_MSage-1] sigma_p_MS[N_pop,N_Mage]; // log-ratio ocean age SDs for each smolt age
  cholesky_factor_corr[N_Mage*(N_MSage-1)] L_p_MS[N_pop]; // Cholesky-factored corr matrices of log-ratio ocean age
  matrix[N,N_Mage*(N_MSage-1)] zeta_p_MS;     //?// log-ratio ocean age errors (Z-scored)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;         // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;        // true broodstock take rate when B_take > 0
  // initial states, observation error
  vector<lower=0>[max_Mage*N_pop] M_init;     // true smolt abundance in years 1:max_Mage
  simplex[N_Mage] q_M_init[max_Mage*N_pop];   // true smolt age distns in years 1:max_Mage
  vector<lower=0>[max_age*N_pop] S_init;      // true total spawner abundance in years 1:max_age
  simplex[N_GRage] q_GR_init[max_age*N_pop];  // true wild spawner age distns in years 1:max_age
  vector<lower=0>[N_pop] tau_M;               // smolt observation error SDs
  vector<lower=0>[N_pop] tau_S;               // spawner observation error SDs
}

transformed parameters {
  //?// indicates transformed params that could be arrays instead of matrices
  // smolt recruitment
  vector<lower=0>[N] M_hat;              // expected smolt abundance (not density) by brood year
  vector[N] epsilon_M;                   // process error in smolt abundance by brood year
  vector<lower=0>[N] M0;                 // true smolt abundance (not density) by brood year
  vector<lower=0>[N] M;                  // true smolt abundance (not density) by outmigration year
  // smolt age structure
  matrix[N_pop,N_Mage-1] gamma_M;        //?// population mean log ratio smolt age distributions
  matrix<lower=0,upper=1>[N,N_Mage] p_M; //?// true smolt age distributions by brood year
  matrix<lower=0,upper=1>[N,N_Mage] q_M; // true smolt age distributions by calendar year
  // SAR
  matrix[N,N_Mage] epsilon_MS;           //?// SAR process errors by smolt age and outmigration year
  matrix<lower=0,upper=1>[N,N_Mage] s_MS; //?// true SAR by smolt age and outmigration year
  vector[N_MSage-1] gamma_MS[N_pop,N_Mage]; // population mean log ratio age distributions
  simplex[N_MSage] p_MS[N,N_Mage];       // true ocean age distns by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  matrix<lower=0,upper=1>[N,N_GRage] q_GR; // true Gilbert-Rich age distns of spawners
  matrix<lower=0,upper=1>[N,N_MSage] q_MS; // true ocean age distns of spawners
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
  {
    // Smolt age 
    gamma_M[j,] = to_row_vector(log(mu_p_M[j,1:(N_Mage-1)]) - log(mu_p_M[j,N_Mage]));
    
    // Ocean age
    for(a in 1:N_Mage)
    {
      gamma_MS[j,a] = log(mu_p_MS[j,a,1:(N_MSage-1)]) - log(mu_p_MS[j,a,N_MSage]);
    }
  }
  
  // Calculate true smolts and total wild and hatchery spawners by age,
  // and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_Mage] alr_p_M;              // temp: alr(p_M[i,])
    row_vector[N_Mage] M_a;                  // temp: true smolts by age
    vector[N_Mage*(N_MSage-1)] gamma_MS_i;   // temp: gamma_MS[i,,] 
    vector[N_Mage*(N_MSage-1)] sigma_p_MS_i; // temp: sigma_p_MS[i,,]
    vector[N_Mage*(N_MSage-1)] alr_p_MS;     // temp: alr(p_MS[i,])    
    matrix[N_Mage,N_MSage] S_W_a;            // temp: true W spawners by smolt and ocean age
    int ii;                                  // temp variable: index into M_init, S_init and q_init
    
    // Time-varying IID age vectors (multivariate Matt trick)
    // Smolt age
    alr_p_M = rep_row_vector(0,N_Mage);
    alr_p_M[1:(N_Mage-1)] = gamma_M[pop[i],] + sigma_p_M[pop[i],] .* (L_p_M[pop[i]] * zeta_p_M[i,]')';
    alr_p_M = exp(alr_p_M);
    p_M[i,] = alr_p_M/sum(alr_p_M);
    
    // Ocean age
    // assemble mean and SD vectors by flattening across smolt age
    for(a in 1:N_Mage)
    {
      gamma_MS_i[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))] = gamma_MS[pop[i],a];
      sigma_p_MS_i[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))] = sigma_p_MS[pop[i],a];
    }
    alr_p_MS = gamma_MS_i + sigma_p_MS_i .* (L_p_MS[pop[i]] * zeta_p_MS[i,]'); //'
    // inverse log-ratio transform and assign back to array
    for(a in 1:N_Mage)
    {
      p_MS[i,a] = exp(append_row(alr_p_MS[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))], 0));
      p_MS[i,a] = p_MS[i,a]/sum(p_MS[i,a]);
    }
    
    // AR(1) smolt recruitment process errors 
    // MAR(1) SAR process errors  
    if(pop_year_indx[i] == 1) // initial process error
    {
      epsilon_M[i] = zeta_M[i]*sigma_M[pop[i]]/sqrt(1 - rho_M[pop[i]]^2);
      // cheat: doesn't use MAR(1) stationary covariance
      epsilon_MS[i,] = zeta_MS[i,] .* sigma_MS[pop[i],] ./ sqrt(1 - square(rho_MS[pop[i],]));
    }
    else
    {
      epsilon_M[i] = rho_M[pop[i]]*epsilon_M[i-1] + zeta_M[i]*sigma_M[pop[i]];
      epsilon_MS[i,] = rho_MS[pop[i],] .* epsilon_MS[i-1,] + sigma_MS[pop[i],] .* (L_MS[pop[i]] * zeta_MS[i,]')';
    }
    // SAR for outmigration year i
    s_MS[i,] = inv_logit(logit(mu_MS[pop[i],]) + dot_product(X_MS[year[i],], beta_MS[pop[i],]) + epsilon_MS[i,]);
    
    // Smolt recruitment
    if(pop_year_indx[i] <= max_Mage)
    {
      // use initial values
      ii = (pop[i] - 1)*max_Mage + pop_year_indx[i];
      M[i] = M_init[ii];  
      q_M[i,] = to_row_vector(q_M_init[ii,]);
    }
    else
    {
      for(a in 1:N_Mage)
      {
        // age-a smolts from appropriate brood year
        M_a[a] = M0[i-smolt_ages[a]]*p_M[i-smolt_ages[a],a]; 
      }
      M[i] = sum(M_a);
      q_M[i,] = M_a/M[i];
    }
    
    // Spawners and age structure
    if(pop_year_indx[i] <= max_age)
    {
      // use initial values
      ii = (pop[i] - 1)*max_age + pop_year_indx[i];
      S_W[i] = S_init[ii]*(1 - p_HOS_all[i]);        
      S_H[i] = S_init[ii]*p_HOS_all[i];
      q_GR[i,] = to_row_vector(q_GR_init[ii,]);
      S_W_a = to_matrix(S_W[i]*q_GR[i,], N_Mage, N_MSage, 0);
      q_MS[i,] = col_sums(to_matrix(q_GR[i,], N_Mage, N_MSage, 0));
    }
    else
    {
      for(sa in 1:N_Mage)
      {
        for(oa in 1:N_MSage)
          S_W_a[sa,oa] = M[i-ocean_ages[oa]]*q_M[i-ocean_ages[oa],sa]*s_MS[i-ocean_ages[oa],sa]*p_MS[i-ocean_ages[oa],sa][oa];
      }
      // catch and broodstock removal (assumes no take of ocean age 1)
      S_W_a[,2:N_MSage] = S_W_a[,2:N_MSage]*(1 - F_rate[i])*(1 - B_rate_all[i]);
      S_W[i] = sum(S_W_a);
      S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
      q_GR[i,] = to_row_vector(S_W_a')/S_W[i]; //'
      q_MS[i,] = col_sums(to_matrix(q_GR[i,], N_Mage, N_MSage, 0));                                
    }
    
    S[i] = S_W[i] + S_H[i];
    
    // Smolt production from brood year i
    M_hat[i] = A[i] * SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    M0[i] = M_hat[i]*exp(dot_product(X_M[year[i],], beta_M[pop[i],]) + epsilon_M[i]);
  }
}

model {
  vector[N_B] B_take; // true broodstock take when B_take_obs > 0
  
  // Priors
  
  // smolt recruitment
  alpha ~ lognormal(2.0,2.0);
  Rmax ~ lognormal(8.0,1.0);
  to_vector(beta_M) ~ normal(0,5);
  rho_M ~ pexp(0,0.85,20); // mildly regularize rho to ensure stationarity
  sigma_M ~ normal(0,5);
  zeta_M ~ std_normal();   // total smolts: log(M) ~ normal(log(M_hat), sigma_M)

  // smolt age structure
  to_vector(sigma_p_M) ~ normal(0,5);
  for(j in 1:N_pop)
    L_p_M[j] ~ lkj_corr_cholesky(1);
  // smolt age probs logistic MVN: 
  // alr(p_M[i,]) ~ MVN(gamma_M[pop[i],], D*R_p_M*D), where D = diag_matrix(sigma_p_M[pop[i],])
  to_vector(zeta_p_M) ~ std_normal();

  // SAR
  to_vector(beta_MS) ~ normal(0,5);
  to_vector(sigma_MS) ~ normal(0,5);
  to_vector(rho_MS) ~ pexp(0,0.85,20);  // mildly regularize rho to ensure stationarity
  for(j in 1:N_pop)
    L_MS[j] ~ lkj_corr_cholesky(1);
  to_vector(zeta_MS) ~ std_normal(); // SAR: logit(s_MS) ~ normal(logit(s_MS_hat), sigma_MS)

  // ocean age structure
  for(j in 1:N_pop)
  {
    for(a in 1:N_Mage)
      sigma_p_MS[j,a] ~ normal(0,5);
    L_p_MS[j] ~ lkj_corr_cholesky(1);
  }
  // ocean age probs logistic MVN: 
  // alr(p_MS[i,]) ~ MVN(gamma_MS[pop[i],,], D*R_p_MS*D), where D = diag_matrix(sigma_p_MS[pop[i],,])
  to_vector(zeta_p_MS) ~ std_normal();

  // removals
  B_take = B_rate .* S_W[which_B] .* (1 - q_MS[which_B,1]) ./ (1 - B_rate);
  B_take_obs ~ lognormal(log(B_take), 0.1); // penalty to force pred and obs broodstock take to match 

  // initial states, observation error
  M_init ~ lognormal(0.0,5.0);
  S_init ~ lognormal(0.0,5.0);
  tau_M ~ lognormal(-3.0,0.2);
  // tau_S ~ pexp(1,0.85,30); // rule out tau < 0.1 to avoid divergences 
  tau_S ~ lognormal(-3.0,0.2);

  // Observation model
  M_obs[which_M_obs] ~ lognormal(log(M[which_M_obs]), tau_M[pop[which_M_obs]]);  // observed smolts
  target += sum(n_Mage_obs .* log(q_M));  // obs wild age freq: n_Mage_obs[i] ~ multinomial(q_M[i])
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S[pop[which_S_obs]]);  // observed spawners
  target += sum(n_MSage_obs .* log(q_MS));  // obs wild age freq: n_MSage_obs[i] ~ multinomial(q_MS[i])
  target += sum(n_GRage_obs .* log(q_GR));  // obs wild age freq: n_age_obs[i] ~ multinomial(q_GR[i])
  n_H_obs ~ binomial(n_HW_obs, p_HOS);  // observed counts of hatchery vs. wild spawners
}

generated quantities {
  corr_matrix[N_Mage-1] R_p_M[N_pop]; // correlation matrices of log-ratio smolt age distns
  corr_matrix[N_Mage] R_MS[N_pop]; // correlation matrices of logit SAR by smolt age
  corr_matrix[N_Mage*(N_MSage-1)] R_p_MS[N_pop]; // correlation matrices of log-ratio ocean age distns
  vector[N] LL_M_obs;           // pointwise log-likelihood of smolts
  vector[N] LL_n_smolt_age_obs; // pointwise log-likelihood of smolt age frequencies
  vector[N] LL_S_obs;           // pointwise log-likelihood of spawners
  vector[N] LL_n_ocean_age_obs; // pointwise log-likelihood of ocean age frequencies
  vector[N] LL_n_GR_age_obs;    // pointwise log-likelihood of Gilbert-Rich age frequencies
  vector[N_H] LL_n_H_obs;       // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL;                 // total pointwise log-likelihood                              
  
  for(j in 1:N_pop)
  {
    R_p_M[j] = multiply_lower_tri_self_transpose(L_p_M[j]);
    R_MS[j] = multiply_lower_tri_self_transpose(L_MS[j]);
    R_p_MS[j] = multiply_lower_tri_self_transpose(L_p_MS[j]);
  }
  
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M); 
  LL_n_smolt_age_obs = (n_Mage_obs .* log(q_M)) * rep_vector(1,N_Mage);
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S); 
  LL_n_ocean_age_obs = (n_MSage_obs .* log(q_MS)) * rep_vector(1,N_MSage);
  LL_n_GR_age_obs = (n_GRage_obs .* log(q_GR)) * rep_vector(1,N_GRage);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_M_obs + LL_n_smolt_age_obs + LL_S_obs + LL_n_ocean_age_obs + LL_n_GR_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
