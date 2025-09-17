data {
  
  // ---- Dimensions ----
  int<lower=1> n;                     // number of sites
  int<lower=1> p;                     // number of covariates 
  int<lower=1> n_z;                   // number of presence/absence sites
  int<lower=1> n_N;                   // number of abundance sites

  // ---- Matrices ----
  matrix[n, n] D;                     // distance matrix 
  matrix[n, p] X;                     // covariate matrix 

  // ---- Data ----
  array[n_z]int<lower=0, upper=1> z;  // presence/absence observations 
  array[n_N]int<lower=0> N;           // abundance observations 
  
  // ---- Fixed Parameters ----
  real<lower=0> phi;                  // spatial range (treating as known to facilitate comparison)

}

transformed data {
  
  // geostatistical covariance 
  matrix[n, n] R = exp(-D / phi);
  matrix[n, n] L_R = cholesky_decompose(R);

}


parameters {
  
  // regression coefficients 
  vector[p] beta;
  // spatial parameters 
  real<lower=0> sigma2;
  // linear predictor 
  vector[n] log_lambda;
  
}

transformed parameters {
  
  // intensity 
  vector[n_N] lambda;
  
  // probability of occupancy 
  vector[n_z] pi;
  
  // presence/absence
  for (i in 1:n_z) {
    pi[i] = 1-exp(-exp(log_lambda[i]));
  }
  
  // abundance 
  for(i in 1:n_N){
    lambda[i] = exp(log_lambda[i+n_z]);
  }

}

model {
  
  //
  // Priors 
  //
  
  // variance parameters 
  sigma2 ~ inv_gamma(0.001, 0.001);

  // regression coefficients 
  beta ~ normal(0, 2.25);


  //
  // Process
  //
  
  log_lambda ~ multi_normal_cholesky(X * beta, sqrt(sigma2) * L_R);
  
  //
  // Likelihood
  // 
  
  // abundance likelihood
  for (i in 1:n_N) {
    target += poisson_lpmf(N[i] | lambda[i]);
  }
  
  // presence/absence likelihood
  for (i in 1:n_z) {
    target += bernoulli_lpmf(z[i] | pi[i]);
  }
  

}
