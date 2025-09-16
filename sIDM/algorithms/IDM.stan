functions {
  
  // AR(1) covariance matrix 
  matrix ar1_cov_matrix(int T, real rho, real sigma2_xi) {
    matrix[T, T] Sigma;
    for (i in 1:T) {
      for (j in 1:T) {
        Sigma[i, j] = sigma2_xi * pow(rho, abs(i - j));
      }
    }
    return Sigma;
  }
  
  // ICAR covariance matrix 
  matrix build_icar_precision(int n_zeta, real rho, matrix A) {
    matrix[n_zeta, n_zeta] Q;
    for (i in 1:n_zeta) {
      for (j in 1:n_zeta) {
        if (i == j)
          Q[i, j] = sum(A[i, ]);  // Diagonal = row sum
        else
          Q[i, j] = -rho * A[i, j];  // Off-diagonal = -rho * A
      }
    }
    
    return Q;
  }
  
real depletion_lpmf_partial_sum(
  array[] int slice_idx,
  int start,
  int end,
  array[] int y_sum_flat,
  array[] int N_max_flat,
  array[] int J_y_flat,
  array[] int y_idx,
  array[] int y_flat,
  array[] int N_deplete_flat,
  vector log_lambda,
  vector logit_p
) {
  real lp = 0;
  int n_slices = end - start + 1;  // number of elements in slice_idx
  for (i in 1:n_slices) {
    int site = slice_idx[i];  // get the actual site index

    real lambda = exp(log_lambda[site]);
    real p = logit_p[site];
    vector[N_max_flat[site]] lps;

    for (l in 1:N_max_flat[site]) {
      lps[l] = poisson_lpmf(y_sum_flat[site] + l - 1 | lambda);

      for (j in 1:J_y_flat[site]) {
        lps[l] += binomial_logit_lpmf(
          y_flat[y_idx[site] + j] | 
          y_sum_flat[site] - N_deplete_flat[y_idx[site] + j] + l - 1,
          p);
      }
    }
    lp += log_sum_exp(lps);
  }
  return lp;
}
  
}

data {
  // ---- Dimensions ----
  int<lower=1> n;                     // Number of sites
  int<lower=1> T;                     // Number of years
  int<lower=1> Q;                     // Number of stream networks
  int<lower=1> p_X;                   // Number of covariates (abundance)
  int<lower=1> p_H;                   // Number of covariates (detection)
  int<lower=1> n_zeta;                // Number of watersheds

  // ---- Watershed variables and indices ----
  matrix[n, n] D;                     // Tail-down covariance matric 
  matrix[n_zeta, n_zeta] A;           // Watershed adjacency matrix
  matrix[n_zeta, n_zeta] zeta_COV;
  array[n] int<lower=1> q_idx;        // Stream network index for each site
  array[Q] int<lower=1> n_q;          // Number of sites in each network 
  array[n] int<lower=1> zeta_idx;     // Watershed index for each site

  // ---- Depletion model ----
  int<lower=1> n_y;                       // Number of valid sites and years in depletion surveys 
  int<lower=1> n_y_all;                   // Number of all passes across all sites and years 
  array[n_y_all] int<lower=0> y_flat;     // Flattened depletion counts
  array[n_y_all] int<lower=0> N_deplete_flat; // Flattened depletion history
  array[n_y] int<lower=0> y_sum_flat;     // Flattened summed depletion counts
  array[n_y] int<lower=1> N_max_flat;     // Flattened N_max
  array[n_y] real<lower=0> A_y_flat;      // Flattened sample area
  array[n_y] int<lower=1> y_site_idx;     // Site index for each depletion obs
  array[n_y] int<lower=1> y_year_idx;     // Year index for each depletion obs
  array[n_y] int<lower=1> J_y_flat;       // Number of depletion occasions per site year
  array[n_y] int<lower=0> y_idx;          // Start index location of each site year
  matrix[n_y, p_X] X_y_flat;              // Flattened abundance covariates
  matrix[n_y, p_H] H_y_flat;              // Flattened detection covariates

  // ---- Presence model ----
  int<lower=1> n_z;                   // Number of valid presence observations
  array[n_z] int<lower=0> z_flat;     // Flattened detections
  array[n_z] real<lower=0> A_z_flat;  // Flattened sample area
  array[n_z] int<lower=1> z_site_idx; // Site index for each presence obs
  array[n_z] int<lower=1> z_year_idx; // Year index for each presence obs
  matrix[n_z, p_X] X_z_flat;          // Flattened abundance covariates
}

transformed data {
  
  // Identity matrix 
  matrix[n, n] I = diag_matrix(rep_vector(1.0, n));
  
}

parameters {
  // Regression Coefficients 
  vector[p_X] beta;
  vector[p_H] gamma;
  // Temporal random effects
  vector[T] xi_raw;
  real<lower=0, upper=1> rho_xi;
  real<lower=0> sigma2_xi;
  // Spatial random effects
  real<lower=0> sigma2_zeta;
  vector[n_zeta] zeta_raw;
  // Tail-down random effect
  vector[n] eta_raw;
  vector<lower=0>[Q] sigma2_eta;
  vector<lower=0>[Q] phi;
}

transformed parameters {
  
  // Linear Predictors 
  vector[n_y] logit_p;
  vector[n_y] log_lambda;
  vector[n_z] log_prob;
  
  // Center random effects (to be identifiable)
  vector[T] xi=xi_raw-mean(xi_raw);
  vector[n_zeta] zeta=zeta_raw-mean(zeta_raw);
  vector[n] eta=eta_raw-mean(eta_raw);
  
  //
  // Calculate linear predictors 
  //
  
  // Detection
  logit_p = H_y_flat * gamma; 
  
  // presence/absence
  for (i in 1:n_z) {
    log_prob[i] = log(A_z_flat[i]) + X_z_flat[i,] * beta + eta_raw[z_site_idx[i]] + zeta_raw[zeta_idx[z_site_idx[i]]] + xi_raw[z_year_idx[i]];
  }
  
  // depletion
  for(i in 1:n_y){
    log_lambda[i] = log(A_y_flat[i]) + X_y_flat[i,] * beta + eta_raw[y_site_idx[i]] + zeta_raw[zeta_idx[y_site_idx[i]]] + xi_raw[y_year_idx[i]];
  }

}

model {
  
  //
  // Priors 
  //
  
  // variance parameters 
  sigma2_zeta ~ inv_gamma(1, 1);
  sigma2_xi ~ inv_gamma(1, 1);
  sigma2_eta ~ inv_gamma(1, 1);
  
  // spatial range 
  phi ~ gamma(1, 1);

  // regression coefficients 
  beta ~ normal(0, 1.5);
  gamma ~ normal(0, 1.5);
  
  //
  // Process
  //
  
  profile("Spatial") {
  
  // Spatial Dependence (across watersheds)
  zeta_raw ~ multi_normal(rep_vector(0.0, n_zeta), sigma2_zeta * zeta_COV);
  
  // Spatial Dependence (within stream network)
{
  int pos;                 // starting position in flat vector
  pos = 1;
  for (q in 1:Q) {
    int nq = n_q[q];       // size of this network
    vector[nq] mu = rep_vector(0.0, nq);
    vector[nq] eta_sub = segment(eta_raw, pos, nq);
    matrix[nq, nq] D_sub = block(D, pos, pos, nq, nq);
    matrix[nq, nq] Sigma;

    // build covariance
    for (i in 1:nq) {
      for (j in 1:nq) {
        Sigma[i, j] = sigma2_eta[q] * exp(-D_sub[i, j] / phi[q]);
      }
    }

    eta_sub ~ multi_normal(mu, Sigma);

    pos += nq;             // update position in flat vector
  }
}

  }
  
  profile("Temporal") {
  
  // Temporal Dependence 
  matrix[T,T] L = cholesky_decompose(ar1_cov_matrix(T, rho_xi, sigma2_xi));
  xi_raw ~ multi_normal_cholesky(rep_vector(0, T), L);
  
  }

  //
  // Likelihood
  // 
  
  profile("Depletion") {

    target += reduce_sum(depletion_lpmf_partial_sum, 
                         linspaced_int_array(n_y, 1, n_y), // slice indices
                         1,
                         y_sum_flat,
                         N_max_flat,
                         J_y_flat,
                         y_idx,
                         y_flat,
                         N_deplete_flat,
                         log_lambda,
                         logit_p);
  
  }
  
  profile("Presence") {
  
  // Presence Likelihood (cloglog link)
  for (i in 1:n_z) {
    real prob = -expm1(-exp(log_prob[i]));  // numerically stable
    target += bernoulli_logit_lpmf(z_flat[i] | prob);
  }
  
  }

}

