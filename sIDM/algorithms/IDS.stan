functions {
  real weighted_poisson_loglik(row_vector pi, real lambda) {
    int M = num_elements(pi) - 1;
    real log_lambda = log(lambda);
    real log_lik = negative_infinity();

    for (k in 0:M) {
      real log_poisson = -lambda + k * log_lambda - lgamma(k + 1);
      if (pi[k + 1] > 0) {
        real log_weighted = log(pi[k + 1]) + log_poisson;
        log_lik = log_sum_exp(log_lik, log_weighted);
      }
    }
    return log_lik;
  }
}

data {
  int<lower=1> n;              // Total number of sites
  int<lower=1> n_z;            // Sites with binary detection
  int<lower=1> n_N;            // Sites with count data
  int<lower=1> T;              // Number of years
  int<lower=1> p_X;            // Number of covariates
  int<lower=1> M;              // Max count (superpopulation size)

  array[n_z] int z;            // Binary presence-absence
  matrix[n_N, M + 1] pi;       // Weighted probabilities for count data
  matrix[n, p_X] X;            // Site-level covariates
  matrix[n, T] year_mat;       // Year design matrix
  vector[n] S;                 // Survey areas (not used in model yet)
}

parameters {
  vector[p_X] beta;
  vector[T] xi_raw;
  real<lower=0> tau;
}

transformed parameters {
  vector[n] log_lambda;
  vector[n] lambda;
  vector[T] xi = xi_raw - mean(xi_raw);  // Centered temporal random effects

  for (i in 1:n) {
    log_lambda[i] = X[i] * beta + year_mat[i] * xi_raw;
    lambda[i] = S[i] * exp(log_lambda[i]);
  }
}

model {
  // Priors
  tau ~ inv_gamma(0.001, 0.001);
  beta ~ normal(0, 1.5);
  xi_raw ~ normal(0, sqrt(tau));

  // Binary detection (distance sampling-like component)
  for (i in 1:n_z) {
    real prob = -expm1(-lambda)[i];  // numerically stable
    z[i] ~ bernoulli(prob);
  }

  // Count-based (presence) likelihood using weighted Poisson
  for (i in 1:n_N) {
    target += weighted_poisson_loglik(pi[i], lambda[i + n_z]);
  }
}
