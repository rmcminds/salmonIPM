functions {
  #include /include/SR.stan
  #include /include/pexp_lpdf_vec.stan
  #include /include/quantile.stan
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // brood year identifier
  // recruitment
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  vector<lower=0>[N] A;                // habitat area associated with each spawner abundance obs
  int<lower=0> K_alpha;                // number of intrinsic productivity covariates
  matrix[N,K_alpha] X_alpha;           // intrinsic productivity covariates
  int<lower=0> K_Rmax;                 // number of maximum recruitment covariates
  matrix[N,K_Rmax] X_Rmax;             // maximum recruitment covariates
  int<lower=0> K_R;                    // number of recruitment covariates
  row_vector[K_R] X_R[N];              // brood-year productivity covariates
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed total spawner abundance (not density)
  // spawner age structure
  int<lower=2> N_age;                  // number of maiden adult age classes
  int<lower=2> max_age;                // maximum maiden adult age
  matrix<lower=0>[N,N_age*2] n_MRage_obs; // observed wild spawner age frequencies [maiden | repeat]
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=2> ages[N_age];                // adult ages
  int<lower=1> min_age;                    // minimum maiden adult age
  int<lower=4> N_MRage = N_age*2;          // number of maiden + repeat adult age classes
  int<lower=1> pop_year_indx[N];           // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  real mu_Rmax = quantile(log(S_obs[which_S_obs]), 0.9); // prior log-mean of Rmax
  real sigma_Rmax = sd(log(S_obs[which_S_obs])); // prior log-SD of Rmax
  vector[max_age*N_pop] mu_S_init;         // prior mean of total spawner abundance in years 1:max_age
  real sigma_S_init = 2*sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_age
  matrix[N_MRage,max_age*N_pop] mu_q_init; // prior counts of age distributions in years 1:max_age
  
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
  min_age = min(ages);  
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  
  for(i in 1:max_age)
  {
    int N_orphan_age = N_age - max(i - min_age, 0); // number of orphan maiden age classes
    int N_amalg_age = N_age - N_orphan_age + 1;     // number of amalgamated maiden age classes
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_age + i; // index into S_init, q_MR_init
      // S_init prior mean scales observed log-mean by fraction of orphan maiden age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_age) - log(N_age);
      // prior on q_MR_init that implies q_orphan ~ Dir(1)
      // (second stacked copy is for orphan repeat age classes in year 1,
      //  otherwise the paired maiden and repeat classes are amalgamated)
      mu_q_init[1:N_age,ii] = append_row(rep_vector(0.5/N_amalg_age, N_amalg_age), rep_vector(0.5, N_orphan_age - 1));
      mu_q_init[(N_age+1):N_MRage,ii] = mu_q_init[1:N_age,ii];
    }
  }
}

parameters {
  // recruitment
  vector<lower=0>[N_pop] alpha;           // intrinsic productivity
  matrix[N_pop,K_alpha] beta_alpha;       // regression coefs for log alpha
  vector<lower=0>[N_pop] Rmax;            // maximum recruitment
  matrix[N_pop,K_Rmax] beta_Rmax;         // regression coefs for log Rmax
  matrix[N_pop,K_R] beta_R;               // regression coefs for log recruitment
  vector<lower=-1,upper=1>[N_pop] rho_R;  // AR(1) coefs for log productivity anomalies
  vector<lower=0>[N_pop] sigma_R;         // process error SDs
  vector[N] zeta_R;                       // recruitment process errors (Z-scored)
  // kelt survival
  vector<lower=0,upper=1>[N_pop] mu_SS;   // mean kelt (spawner-to-spawner) survival
  vector<lower=-1,upper=1>[N_pop] rho_SS; // AR(1) coefs for kelt survival
  vector<lower=0>[N_pop] sigma_SS;        // kelt survival process error SDs
  vector[N] zeta_SS;                      // kelt survival process errors (Z-scored)
  // maiden spawner age structure
  simplex[N_age] mu_p[N_pop];             // population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p; // log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop]; // Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;               // log-ratio cohort age distribution errors (Z-scores)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;     // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;    // true broodstock take rate when B_take > 0
  // initial spawners, observation error
  vector<lower=0>[max_age*N_pop] S_init;  // initial total spawner abundance in years 1:max_age
  simplex[N_MRage] q_MR_init[max_age*N_pop]; // initial age distribution in years 1:max_age [maiden | repeat]
  vector<lower=0>[N_pop] tau;             // observation error SDs of total spawners
}

transformed parameters {
  // recruitment
  vector<lower=0>[N] alpha_Xbeta;      // intrinsic productivity including covariate effects
  vector<lower=0>[N] Rmax_Xbeta;       // maximum recruitment including covariate effects
  vector<lower=0>[N] R_hat;            // expected recruit abundance (not density) by brood year
  vector[N] epsilon_R;                 // process error in recruit abundance by brood year
  vector<lower=0>[N] R;                // true recruit abundance (not density) by brood year
  // kelt survival
  vector[N] epsilon_SS;                // process error in kelt survival by outmigration year
  vector<lower=0,upper=1>[N] s_SS;     // true kelt survival by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                 // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;              // true total wild spawner abundance
  vector[N] S_H;                       // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  vector[N_age-1] mu_alr_p[N_pop];     // population mean log ratio age distributions
  simplex[N_age] p[N];                 // true cohort maiden age distributions
  matrix<lower=0,upper=1>[N,N_MRage] q_MR; // true spawner age distributions [maiden | repeat]

  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;

  // S-R parameters including covariate effects
  alpha_Xbeta = alpha[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
  Rmax_Xbeta = Rmax[pop] .* exp(rows_dot_product(X_Rmax, beta_Rmax[pop,]));

  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    mu_alr_p[j] = log(head(mu_p[j], N_age-1)) - log(mu_p[j,N_age]);

  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    int ii; // index into S_init and q_MR_init
    // number of orphan maiden age classes <lower=0,upper=N_age>
    vector[N_age] q_init; // initial age dist after year 1 (maiden and repeat cells amalgamated)
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_age, 0), N_age);
    vector[N_orphan_age] q_orphan; // orphan maiden age distribution (amalgamated simplex)
    vector[N_age] alr_p; // alr(p[i,])
    row_vector[N_age] S_M_a; // true wild maiden spawners by age
    row_vector[N_age] S_R_a; // true wild repeat spawners by age

    // Multivariate Matt trick for within-pop, time-varying maiden age vectors
    alr_p = rep_vector(0,N_age);
    alr_p[1:(N_age-1)] = mu_alr_p[pop[i]] + sigma_p[pop[i],]' .* (L_p[pop[i]] * zeta_p[i,]');
    alr_p = exp(alr_p);
    p[i] = alr_p / sum(alr_p);

    // AR(1) recruitment and kelt survival process errors
    if(pop_year_indx[i] == 1) // initial process error
    {
      epsilon_R[i] = zeta_R[i] * sigma_R[pop[i]] / sqrt(1 - rho_R[pop[i]]^2);
      epsilon_SS[i] = zeta_SS[i] * sigma_SS[pop[i]] / sqrt(1 - rho_SS[pop[i]]^2);
    }
    else
    {
      epsilon_R[i] = rho_R[pop[i]] * epsilon_R[i-1] + sigma_R[pop[i]] * zeta_R[i];
      epsilon_SS[i] = rho_SS[pop[i]] * epsilon_SS[i-1] + sigma_SS[pop[i]] * zeta_SS[i];
    }
    // kelt survival
    s_SS[i] = inv_logit(logit(mu_SS[pop[i]]) + epsilon_SS[i]);

    // Spawners and age structure
    // use initial values for orphan age classes, otherwise use process model
    if(pop_year_indx[i] <= max_age)
    {
      ii = (pop[i] - 1)*max_age + pop_year_indx[i];
      
      // in year 1 all maiden and repeat classes are orphans, so use full initial age dist
      // otherwise amalgamate to give maiden-only dist q_init before further amalgamation
      if(pop_year_indx[i] > 1)
      {
        q_init = q_MR_init[ii][1:N_age] + q_MR_init[ii][(N_age+1):N_MRage];
        q_orphan = append_row(sum(head(q_init, N_age - N_orphan_age + 1)), tail(q_init, N_orphan_age - 1));
      }
    }

    for(a in 1:N_age)
    {
      // Maiden spawners
      if(pop_year_indx[i] <= ages[a]) // use initial values
      {
        if(pop_year_indx[i] == 1) // use full initial age dist
          S_M_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_MR_init[ii][a];
        else // use amalgamated maiden-only age dist
          S_M_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
      }
      else // use recruitment process model
        S_M_a[a] = R[i-ages[a]]*p[i-ages[a],a];

      // Repeat spawners
      if(pop_year_indx[i] == 1) // use initial values
        S_R_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_MR_init[ii][N_age + a];
      else // use recruitment process model
        S_R_a[a] = S_W[i-1]*q_MR[i-1,a]*s_SS[i-1];
    }

    // Total spawners
    // Catch and broodstock removal (assumes no take of age 1)
    S_M_a[2:N_age] = S_M_a[2:N_age]*(1 - F_rate[i])*(1 - B_rate_all[i]);
    S_R_a[2:N_age] = S_R_a[2:N_age]*(1 - F_rate[i])*(1 - B_rate_all[i]);
    S_W[i] = sum(S_M_a + S_R_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];
    q_MR[i,] = append_col(S_M_a, S_R_a)/S_W[i];

    // Recruitment
    R_hat[i] = SR(SR_fun, alpha_Xbeta[i], Rmax_Xbeta[i], S[i], A[i]);
    R[i] = R_hat[i] * exp(dot_product(X_R[i], beta_R[pop[i],]) + epsilon_R[i]);
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0

  // Priors

  // recruitment
  alpha ~ lognormal(1.5,0.5);
  Rmax ~ lognormal(mu_Rmax, sigma_Rmax);
  to_vector(beta_R) ~ normal(0,5);
  rho_R ~ pexp(0,0.85,20); // mildly regularize rho_R to ensure stationarity
  sigma_R ~ normal(0,5);
  zeta_R ~ std_normal();   // total recruits: R ~ lognormal(log(R_hat), sigma_R)

  // kelt survival
  sigma_SS ~ normal(0,3);
  rho_SS ~ pexp(0,0.85,20); // mildly regularize rho_R to ensure stationarity
  zeta_SS ~ normal(0,1);

  // recruit age structure
  to_vector(sigma_p) ~ normal(0,5);
  for(j in 1:N_pop)
    L_p[j] ~ lkj_corr_cholesky(3);
  // age probs logistic MVN: alr_p[i,] ~ MVN(mu_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();

  // removals
  log_B_take = log(S_W[which_B]) + log1m(q_MR[which_B,1]) + logit(B_rate); // B_take = S_W*(1 - q_MR[,1])*B_rate/(1 - B_rate)
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match

  // initial spawners and maiden spawner age distribution
  // (accounting for amalgamation of q_MR_init to q_orphan)
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_MRage,max_age*N_pop] q_init_mat;
    for(j in 1:size(q_MR_init)) q_init_mat[,j] = q_MR_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_MR_init[i] ~ Dir(mu_q_init[,i])
  }

  // spawner observation error
  tau ~ pexp(1,0.85,30);  // rule out tau < 0.1 to avoid divergences

  // Observation model
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau[pop[which_S_obs]]);  // observed total spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // obs counts of hatchery vs. wild spawners
  target += sum(n_MRage_obs .* log(q_MR)); // obs wild spawner age freq: n_MRage_obs[i] ~ multinomial(q_MR[i])
}

generated quantities {
  corr_matrix[N_age-1] R_p[N_pop]; // correlation matrices of within-pop cohort log-ratio age distns
  vector[N] LL_S_obs;              // pointwise log-likelihood of total spawners
  vector[N_H] LL_n_H_obs;          // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_MRage_obs;        // pointwise log-likelihood of wild spawner age age frequencies [maiden | repeat]
  vector[N] LL;                    // total pointwise log-likelihood

  for(j in 1:N_pop)
    R_p[j] = multiply_lower_tri_self_transpose(L_p[j]);

  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau[pop[which_S_obs]]);
  LL_n_MRage_obs = (n_MRage_obs .* log(q_MR)) * rep_vector(1,N_MRage);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_S_obs + LL_n_MRage_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
