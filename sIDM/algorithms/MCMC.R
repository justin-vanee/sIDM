###
### Functions 
###

logit_inv <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  log(p / (1 - p))
}

pad_square_matrices <- function(mat_list) {
  # Determine the maximum dimension
  max_dim <- max(sapply(mat_list, nrow))  # same as ncol for square matrices
  
  # Pad each matrix
  padded_list <- lapply(mat_list, function(mat) {
    d <- nrow(mat)
    padded <- matrix(0, nrow = max_dim, ncol = max_dim)
    padded[1:d, 1:d] <- mat
    return(padded)
  })
  
  return(padded_list)
}

blockdiag_to_list <- function(D, block_sizes) {
  out <- vector("list", length(block_sizes))
  start <- 1
  for (i in seq_along(block_sizes)) {
    end <- start + block_sizes[i] - 1
    out[[i]] <- D[start:end, start:end, drop = FALSE]
    start <- end + 1
  }
  out
}

replace_low_off_diagonals <- function(mat_list, value) {
  lapply(mat_list, function(mat) {
    if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
      stop("All elements must be square matrices.")
    }
    off_diag_idx <- row(mat) != col(mat)
    mat[off_diag_idx & mat < value] <- value
    return(mat)
  })
}

plot_traceplots <- function(mcmc_mat, param_names = NULL, burnin = 0) {
  if (is.null(param_names)) {
    param_names <- paste0("beta[", seq_len(nrow(mcmc_mat)), "]")
  }
  
  # Apply burn-in trimming
  if (burnin > 0) {
    mcmc_mat <- mcmc_mat[, -(1:burnin), drop = FALSE]
  }
  
  # Convert to long format for ggplot
  library(tidyverse)
  trace_df <- as.data.frame(mcmc_mat) %>%
    mutate(Parameter = param_names) %>%
    pivot_longer(
      cols = -Parameter,
      names_to = "Iteration",
      values_to = "Value"
    ) %>%
    mutate(
      Iteration = as.integer(gsub("V", "", Iteration)) + burnin,
      Parameter = factor(Parameter, levels = param_names)
    )
  
  # Plot
  ggplot(trace_df, aes(x = Iteration, y = Value)) +
    geom_line(alpha = 0.6, color = "steelblue") +
    facet_wrap(~Parameter, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(
      title = "Traceplots for MCMC Samples",
      x = "Iteration", y = "Sample Value"
    )
}

make_ar1_corr <- function(T, rho) {
  mat <- matrix(0, nrow = T, ncol = T)
  for (i in 1:T) {
    for (j in i:T) {
      mat[i, j] <- rho^abs(i - j)
      mat[j, i] <- mat[i, j]  # enforce symmetry
    }
  }
  return(mat)
}

fast_rmvn <- function(mean, cov) {
  L <- chol(cov)                      # Cholesky decomposition: cov = L' * L
  z <- rnorm(length(mean))           # Standard normal vector
  mean + backsolve(L, z, transpose = TRUE)  # Equivalent to: mean + L'^{-1} z
}

NB_polson <- function(N, nu, log_lambda){
  psi = logit_inv(log_lambda)
  dense = nu * log(1-psi) + lgamma(N+nu) - lgamma(nu)
  return(dense)
}

NB_polson_N <- function(N, nu, log_lambda){
  psi = logit_inv(log_lambda)
  dense = N * log(psi) + lgamma(N+nu) - lfactorial(N)
  return(dense)
}

NB_polson_II <- function(N, nu, log_lambda){
  psi = logit_inv(log_lambda)
  dense = N * log(psi) + nu * log(1-psi)
  return(dense)
}

make_taildown_spherical_R <- function(D_list, A_list, B_list, range) {
  R_list <- vector("list", length(D_list))
  
  for (i in seq_along(D_list)) {
    D <- D_list[[i]]
    A <- A_list[[i]]
    B <- B_list[[i]]
    
    # Normalize distances
    r  <- D / range
    r1 <- A / range
    r2 <- B / range
    
    # Initialize R matrix
    R <- matrix(0, nrow = nrow(D), ncol = ncol(D))
    
    # Identity matrix for diagonals
    diag(R) <- 1
    
    # Flow-connected: A == 0, D > 0
    conn_mask <- (A == 0) & (D > 0) & (r <= 1)
    R[conn_mask] <- (1 - 1.5 * r[conn_mask] + 0.5 * r[conn_mask]^3)
    
    # Flow-unconnected: A > 0, B > 0, r2 <= 1
    unconn_mask <- (A > 0) & (B > 0) & (r2 <= 1)
    R[unconn_mask] <- (1 - 1.5 * r1[unconn_mask] + 0.5 * r2[unconn_mask]) *
      (1 - r2[unconn_mask])^2
    
    # Save result
    R_list[[i]] <- R
  }
  
  return(R_list)
}

compute_crps_matrices <- function(S_test, X_test, year_mat_test, watershed_mat_test, site_mat_test,
                                  v_test, PG_out, burnin, n_mcmc) {
  
  # Helper function for linear predictor
  predict_log_lambda <- function(S, X, year, watershed, site, beta, xi, zeta, eta) {
    log(S) + X %*% beta + year %*% xi + watershed %*% zeta + site %*% eta
  }
  
  n_pred <- length(v_test)
  n_iter <- n_mcmc - burnin
  
  crps       <- matrix(NA, nrow = n_pred, ncol = n_iter)
  crps_logit <- matrix(NA, nrow = n_pred, ncol = n_iter)
  
  for (l in (burnin + 1):n_mcmc) {
    log_lambda <- predict_log_lambda(
      S_test, X_test, year_mat_test, watershed_mat_test, site_mat_test,
      PG_out$beta[, l], PG_out$xi[, l], PG_out$zeta[, l], PG_out$eta[, l]
    )
    
    idx <- l - burnin
    crps[, idx]       <- (v_test - (1 - exp(-exp(log_lambda))))^2
    crps_logit[, idx] <- (v_test - logit_inv(log_lambda))^2
  }
  
  return(list(crps = crps, crps_logit = crps_logit))
}

extract_trace_df <- function(fit_list, param, rows = 1:4) {
  map2_dfr(fit_list, names(fit_list), ~ {
    mat <- .x[[param]][rows, sample_idx, drop = FALSE]
    mat <- if (is.null(dim(mat))) matrix(mat, nrow = 1) else mat
    colnames(mat) <- paste0("iter_", seq_len(ncol(mat)))
    
    as_tibble(mat) %>%
      mutate(param_row = row_number()) %>%
      pivot_longer(-param_row, names_to = "iter", values_to = "value") %>%
      mutate(
        iter = as.integer(gsub("iter_", "", iter)),
        model = .y,
        param = paste0(param, "_", param_row)
      )
  })
}

extract_scalar_trace_df <- function(fit_list, param) {
  map2_dfr(fit_list, names(fit_list), ~ {
    vec <- .x[[param]][sample_idx]
    tibble(
      iter = seq_along(vec),
      value = vec,
      model = .y,
      param = param
    )
  })
}

plot_trace_df <- function(df, burnin = 4000) {
  df %>%
    filter(iter > burnin) %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line(alpha = 0.5) +
    facet_grid(param ~ model, scales = "free_y") +
    theme_minimal() +
    labs(x = "Iteration", y = "Value", title = "Traceplots by Parameter and Model")
}

extract_log_lambda_means <- function(fit_list) {
  map2_dfr(fit_list, names(fit_list), ~ {
    log_lambda_sub <- .x$log_lambda[, sample_idx, drop = FALSE]
    means <- rowMeans(log_lambda_sub)
    tibble(
      value = means,
      model = .y
    )
  })
}

## Simulate count and presence/absence data from shared process
simulate_data <- 
  function(n, z_prop, p, nu=NULL){
    
    # Calculate number of absence/presence sites
    n_z <- round(z_prop * n)
    n_N <- n - n_z
    z_idx <- 1:n_z
    
    # Simulate n random points in (0,1) x (0,1)
    coords <- matrix(runif(2 * n), ncol = 2)
    
    # Compute the Euclidean distance matrix
    D <- as.matrix(dist(coords))
    
    # Matrices 
    X <- 
      cbind(
        rep(1, n),
        matrix(rnorm(p*n), n, p)
      )
    
    # Parameters values 
    sigma2 = rgamma(1, 1, 1)
    phi = runif(1)
    beta = rnorm(p+1)
    
    # Simulate log_lambda
    R = exp( -D / phi )
    log_lambda = rmvn(1, X %*% beta, sigma2 * R)
    
    # Simulate observations 
    if(is.null(nu)){
      # cloglog link
      prob = 1 - exp(-exp(log_lambda[z_idx]))
      v = rbinom(n_z, 1, prob)
      N = rpois(n_N, exp(log_lambda[-z_idx]))
    } else {
      # Polson 2013 (NB)
      #psi = logit_inv(log_lambda)
      # link function implied by negative binomial model
      #prob = 1 - (1 - psi[z_idx])^nu
      #z = rbinom(n_z, 1, prob)
      # parameterized in terms of probability of failure Polson 2013
      #N = rnbinom(n_N, size=nu, prob=(1-psi[-z_idx])) 
      # Typical mean parameterization 
      lambda = exp(log_lambda)
      # link function implied by negative binomial model
      prob = 1 - (nu / ( nu + lambda ))^nu
      z = rbinom(n_z, 1, prob[z_idx])
      N = rnbinom(n_N, size=nu, mu=lambda[-z_idx]) 
    }
    
    # ---- Assemble list ----
    return(
      list(
        n = n,
        p = p+1, # +1 for intercept
        n_z = n_z,
        n_N = n_N,
        coords = coords,
        D = D,
        X = X,
        z = z,
        N = N,
        # Parameter values
        phi = phi,
        sigma2 = sigma2,
        beta = beta,
        log_lambda = log_lambda,
        nu = nu
      )
    )
    
  }

## Fit simulation models 
fit_models <- 
  function(data,
           iter_sampling, iter_warmup, chains=1,
           n_mcmc, burnin,
           VI_iter, VI_draws,
           Models = c("HMC", "gIDM", "gIDM(proposal)", "VI")) {
    
    # Create tibble for outputs
    metrics_tib <- tibble::tibble(
      Model = Models,
      dataset = data$dataset,
      time = NA,
      cov_beta = NA,
      cov_lambda = NA,
      bias_beta = NA,
      bias_lambda = NA,
      sd_beta = NA,
      sd_lambda = NA,
      ess_beta = NA,
      ess_lambda = NA,
      rejection = NA,
      nu = data$nu
    ) %>% tibble::column_to_rownames("Model")
    
    
    if ("HMC" %in% Models) {
  
    model <- cmdstan_model(
      '/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/algorithms/IDM_sim.stan'
    )
  
    fit <- 
      model$sample(
        data = data,
        chains = chains,
        iter_sampling = iter_sampling,
        iter_warmup = iter_warmup
      )
  
    # --- check success here ---  
    if (all(fit$return_codes() == 0)) {
    
      metrics_tib["HMC", "time"] <- pluck(fit$time(), "total")
    
      # Extract posterior samples (HMC) 
      HMC_posterior_beta <-
        fit$draws(variables = c("beta")) %>%
        aperm(c(3, 1, 2)) %>%
        matrix(nrow = data$p, chains * iter_sampling)
      
      HMC_posterior_log_lambda <-
        fit$draws(variables = c("log_lambda")) %>%
        aperm(c(3, 1, 2)) %>%
        matrix(nrow = data$n, chains * iter_sampling)
    
    } else {
    
      message("⚠️ VI failed (non-zero return code). Leaving metrics as NA.")
      HMC_posterior_beta <- NULL
      MC_posterior_log_lambda <- NULL
    
    }
}

    
    # Fit Geometric Regression 
    if ("gIDM" %in% Models) {
      source("/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/algorithms/MCMC.R")
      mu_beta <- numeric(data$p)
      Sigma_beta <- diag(1.5^2, data$p)
      
      Start <- Sys.time()
      gIDM_out <- PG_joint_sim(
        z = data$v, N = data$N, X = data$X, D = data$D,
        mu_beta = mu_beta, Sigma_beta = Sigma_beta,
        Poisson = FALSE, phi = data$phi,
        n_mcmc = n_mcmc
      )
      End <- Sys.time()
      metrics_tib["gIDM", "time"] <- difftime(End, Start, units = "secs")
    }
    
    # Fit Geometric Regression (proposal)
    if ("gIDM(proposal)" %in% Models) {
      source("/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/algorithms/MCMC.R")
      mu_beta <- numeric(data$p)
      Sigma_beta <- diag(1.5^2, data$p)
      
      Start <- Sys.time()
      gIDM_proposal_out <- PG_joint_sim(
        z = data$v, N = data$N, X = data$X, D = data$D,
        mu_beta = mu_beta, Sigma_beta = Sigma_beta,
        Poisson = TRUE, phi = data$phi,
        n_mcmc = n_mcmc
      )
      End <- Sys.time()
      metrics_tib["gIDM(proposal)", "time"] <- difftime(End, Start, units = "secs")
      metrics_tib["gIDM(proposal)", "rejection"] <- mean(gIDM_proposal_out$rejections)
    }
    
    # Fit VI model (mean field approximation)
    if ("VI" %in% Models) {
      if (!exists("model")) {
        model <- cmdstan_model(
          '/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/algorithms/IDM_sim.stan'
        )
      }
      
      # Run Pathfinder (other VI will often fail)
      pathfinder_fit <- model$pathfinder(
        data = data,
        num_paths = 1
      )
      
      # Extract good initial values
      init_list <- pathfinder_fit$draws(format = "list")
      
      fit_VI <-
        model$variational(
        data = data,
        algorithm = "meanfield",
        iter = VI_iter,
        init = init_list,
        draws = VI_draws
      )
      
      # --- check success here ---
      if (all(fit_VI$return_codes() == 0)) {
        metrics_tib["VI", "time"] <- pluck(fit_VI$time(), "total")
        
        VI_posterior_beta <-
          fit_VI$draws(variables = "beta") %>% t()
        
        VI_posterior_log_lambda <-
          fit_VI$draws(variables = "log_lambda") %>% t()
      } else {
        message("⚠️ VI failed (non-zero return code). Leaving metrics as NA.")
        VI_posterior_beta <- NULL
        VI_posterior_log_lambda <- NULL
      }
      
    }
    
    # Define thinning indices
    thin <- floor(seq(from = burnin, to = n_mcmc, length.out = iter_sampling))
    data$beta <- data$beta[-1] # Remove intercept
    
    # Loop over models and compute metrics
    for (model_name in Models) {
      if (model_name == "HMC" && !is.null(HMC_posterior_beta)) {
        beta_post <- HMC_posterior_beta[-1,] # Remove intercept
        lambda_post <- HMC_posterior_log_lambda
      } else if (model_name == "gIDM(proposal)") {
        beta_post <- gIDM_proposal_out$beta[-1, thin] # Remove intercept
        lambda_post <- gIDM_proposal_out$log_lambda[, thin]
      } else if (model_name == "gIDM") {
        beta_post <- gIDM_out$beta[-1, thin] # Remove intercept
        lambda_post <- gIDM_out$log_lambda[, thin]
      } else if (model_name == "VI" && !is.null(VI_posterior_beta)) {
        beta_post <- VI_posterior_beta[-1,] # Remove intercept
        lambda_post <- VI_posterior_log_lambda
      } else {
        next
      }
      
      # Coverage
      lwr <- apply(lambda_post, 1, quantile, 0.025)
      upr <- apply(lambda_post, 1, quantile, 0.975)
      metrics_tib[model_name, "cov_lambda"] <- mean(lwr <= data$log_lambda & data$log_lambda <= upr)
      
      lwr <- apply(beta_post, 1, quantile, 0.025)
      upr <- apply(beta_post, 1, quantile, 0.975)
      metrics_tib[model_name, "cov_beta"] <- mean(lwr <= data$beta & data$beta <= upr)
      
      # Bias
      metrics_tib[model_name, "bias_lambda"] <- mean(abs(data$log_lambda - rowMeans(lambda_post)))
      metrics_tib[model_name, "bias_beta"]   <- mean(abs(data$beta - rowMeans(beta_post)))
      
      # Posterior SD
      metrics_tib[model_name, "sd_lambda"] <- mean(apply(lambda_post, 1, sd))
      metrics_tib[model_name, "sd_beta"]   <- mean(apply(beta_post, 1, sd))
      
      # ESS
      metrics_tib[model_name, "ess_lambda"] <- mean(apply(lambda_post, 1, mcmcse::ess))
      metrics_tib[model_name, "ess_beta"]   <- mean(apply(beta_post, 1, mcmcse::ess))
    }
    
    # Clear up memory 
    fit_VI=fit=gIDM_out=gIDM_proposal_out=NULL
    beta_post=lambda_post=VI_posterior_beta=VI_posterior_log_lambda=HMC_posterior_beta=HMC_posterior_log_lambda=NULL
    gc()
    
    # Return final tibble with "Model" column restored
    metrics_tib <- metrics_tib %>% tibble::rownames_to_column("Model")
    return(metrics_tib)
  }

## Fit 3 models to real data 
fit_models_real_parallel <- function(n_mcmc = 2e4, 
                                     c_1 = 0.7, 
                                     n_cores = 3) {
  
  ### Packages
  library(tidyverse)
  library(scoringRules)
  library(Matrix)
  library(foreach)
  library(doParallel)
  
  ### Load data
  load("/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/Outputs/data.RData")
  
  ### Priors & Hyperparameters
  mu_beta <- numeric(data$p_X)
  Sigma_beta <- diag(1.5^2, data$p_X)
  mu_gamma <- numeric(data$p_H)
  Sigma_gamma <- diag(1.5^2, data$p_H)
  
  ### Indices
  z_zeta_idx <- data$zeta_idx[data$z_site_idx]
  y_zeta_idx <- data$zeta_idx[data$y_site_idx]
  zeta_idx <- c(z_zeta_idx, y_zeta_idx)
  year_idx <- c(data$z_year_idx, data$y_year_idx)
  site_idx <- c(data$z_site_idx, data$y_site_idx)
  y_idx <- data$y_idx
  
  ### Dimensions
  year_levels <- sort(unique(year_idx))
  watershed_levels <- sort(unique(zeta_idx))
  site_levels <- sort(unique(site_idx))
  T <- length(year_levels)
  n <- length(site_levels)
  J_y <- data$J_y_flat
  n_q <- table(data$q_idx[unique(site_idx)])
  
  ### Matrices
  A <- data$A[watershed_levels, watershed_levels]
  eta_COV <- bdiag(data$eta_COV)
  eta_COV <- as.matrix(eta_COV[site_levels, site_levels])
  D <- data$D %>% 
    bdiag() %>% 
    as.matrix() %>% 
    .[site_levels, site_levels] %>% 
    blockdiag_to_list(n_q)
  X <- rbind(data$X_z_flat, data$X_y_flat)
  S <- c(data$A_z_flat, data$A_y_flat)
  H <- data$H_y_flat
  
  ### Model matrices
  year_mat <- model.matrix(~ factor(year_idx, levels = year_levels) - 1)
  watershed_mat <- model.matrix(~ factor(zeta_idx, levels = watershed_levels) - 1)
  site_mat <- model.matrix(~ factor(site_idx, levels = site_levels) - 1)
  
  ### Data
  z <- data$z_flat
  y <- data$y_flat
  N_deplete <- data$N_deplete_flat
  y_sum <- data$y_sum_flat
  
  ### Register parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  ### Parallel model fitting
  model_fits <- foreach(i = 1:3, .packages = c("Matrix")) %dopar% {
    
    source("/Users/justinvanee/Library/Mobile Documents/com~apple~CloudDocs/Documents/brook_trout_postdoc/algorithms/MCMC.R")
    
    start <- Sys.time()
    
    fit <- switch(i,
                  ### Load algorithms 
                  
                  # Model 1: Joint
                  PG_joint(
                    z = z, y = y, J_y = J_y, y_idx = y_idx, y_sum = y_sum, N_deplete = N_deplete,
                    H = H, X = X, D = D, A = A, S = S,
                    year_mat = year_mat, watershed_mat = watershed_mat, site_mat = site_mat,
                    mu_gamma = mu_gamma, Sigma_gamma = Sigma_gamma,
                    mu_beta = mu_beta, Sigma_beta = Sigma_beta,
                    var_shape = 2, var_rate = 1.5^2,
                    c_1 = c_1, n_mcmc = n_mcmc
                  ),
                  
                  # Model 2: Binary
                  PG_binary(
                    z = z, 
                    X = X, D = D, A = A, S = S,
                    year_mat = year_mat, watershed_mat = watershed_mat, site_mat = site_mat,
                    mu_beta = mu_beta, Sigma_beta = Sigma_beta,
                    var_shape = 2, var_rate = 1.5^2,
                    c_1 = c_1, n_mcmc = n_mcmc
                  ),
                  
                  # Model 3: Count
                  PG_count(
                    y = y, J_y = J_y, y_idx = y_idx, y_sum = y_sum, N_deplete = N_deplete,
                    H = H, X = X, D = D, A = A, S = S,
                    year_mat = year_mat, watershed_mat = watershed_mat, site_mat = site_mat,
                    mu_gamma = mu_gamma, Sigma_gamma = Sigma_gamma,
                    mu_beta = mu_beta, Sigma_beta = Sigma_beta,
                    var_shape = 2, var_rate = 1.5^2,
                    c_1 = c_1, n_mcmc = n_mcmc
                  )
    )
    
    end <- Sys.time()
    list(fit = fit, time = difftime(end, start, units = "mins"))
  }
  
  stopCluster(cl)
  
  ### Return results
  PG_joint_out  <- model_fits[[1]]$fit
  PG_binary_out <- model_fits[[2]]$fit
  PG_count_out  <- model_fits[[3]]$fit
  time          <- sapply(model_fits, function(x) x$time)
  
  return(list(PG_count_out=PG_count_out, PG_joint_out=PG_joint_out, PG_binary_out=PG_binary_out, time=time))
}

###
### MCMC Algorithms
###

###
### Simulation Study
###

## Pólya-gamma (sim)
PG_joint_sim <- 
  function(# Data
    z, N, 
    X, D,
    # Priors 
    mu_beta, Sigma_beta,
    var_shape=0.001, var_rate=0.001,
    Poisson=FALSE, phi, 
    # Hyperparameters 
    n_mcmc){
    
    ###
    ### Packages
    ###
    
    library(BayesLogit)
    library(Boom)
    library(Matrix)
    
    ###
    ### Loop Variables
    ###
    
    ### Dimensions 
    n_z = length(z)
    n_N = length(N)
    p = ncol(X)
    n = nrow(X)
    
    ### Priors 
    Sigma_beta_inv = solve(Sigma_beta)
    Sigma_beta_inv_times_mu = Sigma_beta_inv %*% mu_beta
    
    ### Objects 
    kappa = z - rep(1, n_z) / 2
    I_n = diag(1, n)
    q_d = n / 2 + var_shape
    X_N = X[-(1:n_z),]
    M = cbind(X, I_n)
    M_z = M[1:n_z,]
    M_N = M[-(1:n_z),]
    mean_vec = c(Sigma_beta_inv_times_mu, numeric(n))
    
    ###
    ### Starting Values
    ###
    
    beta = mu_beta
    eta = numeric(n)
    log_lambda = numeric(n)
    sigma2_eta = 1
    nu = 1
    R = exp( -D / phi)
    R_inv = solve(R)
    
    ###
    ### Acceptance
    ###
    
    nu_acc = 0 
    nu_tune = 2.4^2
    ### Hyperparameters for tuning adaptive sampler 
    r_opt = 0.234
    rejections = numeric(n_mcmc)
    rejections = rejections - 1
    
    ###
    ### Storage
    ###
    
    log_lambda_save = matrix(NA, n, n_mcmc)
    beta_save = matrix(NA, p, n_mcmc)
    eta_save = matrix(NA, n, n_mcmc)
    sigma2_eta_save = numeric(n_mcmc)
    
    ###
    ### MCMC loop
    ###
    
    for (q in 1:n_mcmc) {
      
      ### Sample omegas
      omega1 = rpg(n_z, 1, log_lambda[1:n_z])
      omega2 = rpg(n_N, N+nu, log_lambda[-(1:n_z)])
      
      ### Sample beta & eta
      var_inv = as.matrix(bdiag(Sigma_beta_inv, (1 / sigma2_eta) * R_inv))
      tmp_chol = chol(t(M_z) %*% (omega1 * M_z) + t(M_N) %*% (omega2 * M_N) +  var_inv)
      
      if(Poisson){ 
        
        # Calculate repeat quantities for efficiency 
        mh2 <-
          # Add log(Likelihood) Geometric 
          sum(NB_polson_II(N, 1, log_lambda[-(1:n_z)])) +
          sum(dbinom(z, 1, logit_inv(log_lambda[(1:n_z)]), TRUE)) -
          # Subtract log(Likelihood) Poisson 
          sum(dpois(N, exp(log_lambda[-(1:n_z)]), TRUE)) -
          sum(dbinom(z, 1, prob=1-exp(-exp(log_lambda[(1:n_z)])), TRUE)) 

        
        repeat {
          
          # Count number of samples 
          rejections[q] <- rejections[q] + 1
          
          coef = backsolve(tmp_chol, 
                           backsolve(tmp_chol,  
                                     t(M_z) %*% kappa + t(M_N) %*% (N-nu)/2 + mean_vec,
                                     transpose=TRUE) + rnorm(p+n)) 
          
          log_lambda_star = M %*% coef
          
          mh1 <-
            # Add log(Likelihood) Poisson 
            sum(dpois(N, exp(log_lambda_star[-(1:n_z)]), TRUE)) +
            sum(dbinom(z, 1, prob=1-exp(-exp(log_lambda_star[(1:n_z)])), TRUE)) -
            # Subtract log(Likelihood) Geometric 
            sum(NB_polson_II(N, 1, log_lambda_star[-(1:n_z)])) -
            sum(dbinom(z, 1, logit_inv(log_lambda_star[(1:n_z)]), TRUE)) 
          
          # MH ratio
          mh = exp(mh1 + mh2)
          
          # Accept or keep looping
          if (runif(1) < mh | rejections[q] > 1e6) {
            log_lambda = log_lambda_star
            break
          }
        }
        
      } else {
        coef = backsolve(tmp_chol, 
                         backsolve(tmp_chol,  
                                   t(M_z) %*% kappa + t(M_N) %*% (N-nu)/2 + mean_vec,
                                   transpose=TRUE) + rnorm(p+n)) 
        log_lambda = M %*% coef
      }
      beta = coef[1:p]
      eta = coef[-(1:p)]
      
      
      ### Sample sigma2_eta
      r_tilde = t(eta) %*% R_inv %*% eta / 2 + var_rate
      sigma2_eta = 1 / rgamma(1, shape = q_d, scale = 1 / r_tilde)
      
      ### Save Samples
      log_lambda_save[, q] = log_lambda
      beta_save[, q] = beta
      eta_save[, q] = eta
      sigma2_eta_save[q] = sigma2_eta
      
      ### Timer
      if(Poisson){
        if (q %% 10 == 0) cat(q, " ")
      } else {
        if (q %% 1000 == 0) cat(q, " ")
      }
    }
    
    list(
      ## MCMC samples
      log_lambda=log_lambda_save,
      beta=beta_save,
      eta=eta_save,
      sigma2_eta=sigma2_eta_save,
      rejections=rejections
    )
  }

###
### Kery et al. (2024) IDS
###

## Distance-sampling model only (Stage I)
ds_MCMC <- function(
    # Data
  Y, D, H, 
  # Priors 
  u_d, 
  mu_gamma, Sigma_gamma,
  # Hyperparameters 
  c_1=0.7, n_mcmc){   
  
  ####  
  #### Setup Variables  
  ####  
  
  M <- nrow(Y) 
  n_N <- ncol(Y)
  
  # Storage
  gamma_save <- matrix(NA, nrow=ncol(H), ncol=n_mcmc)
  chi_save <- matrix(NA, nrow=n_N, ncol=n_mcmc)
  N_save <- matrix(NA, nrow=n_N, ncol=n_mcmc)
  
  #### Starting Values
  gamma <- glm(colSums(Y) ~ H - 1, family = "poisson")$coefficients
  sigma2 <- t(exp(H %*% gamma))
  V <- Y
  D_miss <- is.na(D)
  n_D_miss <- sum(D_miss)
  D[D_miss] <- runif(n_D_miss, 0, u_d)
  chi <- colMeans(Y)
  Chi <- matrix(rep(chi, each = M), nrow = M)
  P <- exp(-sweep(D^2, 2, sigma2, FUN = "/"))
  
  #### Tuning
  gamma_acc = 0
  gamma_tune = 2.4^2 / 100
  Sigma_tune = diag(1, ncol(H))
  ### Hyperparameters for tuning adaptive sampler 
  r_opt = 0.234
  
  #### Begin MCMC Loop
  for (q in 1:n_mcmc) {
    
    ### For automatic tuning 
    autotune = 1.0 / q^c_1
    
    #### Sample V
    chi_tmp <- Chi * (1 - P) / (Chi * (1 - P) + 1 - Chi)
    V[Y == 0] <- rbinom(n_D_miss, 1, chi_tmp[Y == 0])
    
    #### Sample chi
    tmp_sum <- apply(V, 2, sum)
    chi <- rbeta(n_N, tmp_sum + 1, M - tmp_sum + 1) # Discrete-Uniform(0,M) implied for N
    Chi <- matrix(rep(chi, each = M), nrow = M)
    
    #### Sample gamma
    gamma_star <- Boom::rmvn(1, gamma, gamma_tune * Sigma_tune)
    P_star <- exp(-sweep(D^2, 2, exp(H %*% gamma_star), FUN = "/"))
    mh1 <- sum(dbinom(Y[V == 1], 1, P_star[V == 1], log = TRUE))
    mh2 <- sum(dbinom(Y[V == 1], 1, P[V == 1], log = TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      gamma <- gamma_star
      P <- P_star
      sigma2 <- t(exp(H %*% gamma))
      gamma_acc <- gamma_acc + 1
      gamma_tune <- exp(log(gamma_tune) + autotune*(1.0 - r_opt))
    } else {
      gamma_tune <- exp(log(gamma_tune) + autotune*(0.0 - r_opt))
    }
    
    #### Sample d
    D_star <- D
    D_star[D_miss] <- runif(n_D_miss, 0, u_d)
    p_star <- exp(-sweep(D_star^2, 2, sigma2, FUN = "/")[D_miss])
    mh1 <- V[D_miss] * log(1 - p_star)
    mh2 <- V[D_miss] * log(1 - P[D_miss])
    keep_idx <- exp(mh1 - mh2) > runif(n_D_miss)
    keep_idx[is.na(keep_idx)] <- FALSE # Fixes numerical issues
    D[D_miss][keep_idx] <- D_star[D_miss][keep_idx]
    P[D_miss][keep_idx] <- p_star[keep_idx]
    
    #### Save Samples
    gamma_save[,q] <- gamma
    chi_save[,q] <- chi
    N_save[,q] <- colSums(V)
    
    #### Update covariance tuning matrix
    if(q>=1000){
      Sigma_tune <- Sigma_tune + autotune * (cov(t(gamma_save[,1:q])) - Sigma_tune)
    }
    
    #### Timer
    if (q %% 1000 == 0) cat(q, " ")
  }
  
  #### Write Output
  list(
    gamma = gamma_save,
    chi = chi_save,
    N = N_save,
    gamma_acc = gamma_acc,
    Sigma = gamma_tune * Sigma_tune
  )
}


## Joint model (Stage II)
ds_II_MCMC <- 
  function(# Data
    z, N_I, X, S, year_mat,
    # Priors 
    mu_beta, Sigma_beta,
    var_shape=0.001, var_rate=0.001,
    Poisson=FALSE,
    # Hyperparameters 
    n_mcmc){
    
    ###
    ### Packages
    ###
    
    library(BayesLogit)
    library(Boom)
    library(Matrix)
    
    ###
    ### Loop Variables
    ###
    
    ### Dimensions 
    n_z = length(z)
    n_N = nrow(N_I)
    n_samp = ncol(N_I)
    p_X = ncol(X)
    T = ncol(year_mat)
    
    ### Priors 
    Sigma_beta_inv = solve(Sigma_beta)
    Sigma_beta_inv_times_mu = Sigma_beta_inv %*% mu_beta
    ### How we effectively set the coefficient on area (S) to 1 
    coef_precision = 1e5
    
    ### Objects 
    kappa = z - rep(1, n_z) / 2
    I = diag(1, n_z)
    I_T = diag(1, T)
    M = cbind(log(S), X, year_mat)
    M_z = M[1:n_z,]
    M_N = M[-(1:n_z),]
    
    ###
    ### Starting Values & Loop Variables
    ###
    
    q_xi = T / 2 + var_shape
    beta = mu_beta
    xi = numeric(T)
    log_lambda = numeric(n_z+n_N)
    N = N_I[,n_samp]
    tau = 1
    nu = 1
    mean_vec = c(coef_precision, Sigma_beta_inv_times_mu, rep(0, T))
    coef_inv = diag(coef_precision, 1)
    
    ###
    ### Acceptance
    ###
    
    N_acc = numeric(n_N)
    
    ###
    ### Storage
    ###
    
    beta_save = matrix(NA, p_X, n_mcmc)
    N_save = matrix(NA, n_N, n_mcmc)
    xi_save = matrix(NA, T, n_mcmc)
    log_lambda_save = matrix(NA, n_z+n_N, n_mcmc)
    tau_save = numeric(n_mcmc)
    coef_save = numeric(n_mcmc)
    
    ###
    ### MCMC loop
    ###
    
    for (q in 1:n_mcmc) {
      
      ###
      ### Sample fixed and random effects 
      ###
      
      ### Sample omegas
      omega1 = rpg(n_z, 1, log_lambda[1:n_z])
      omega2 = rpg(n_N, N+nu, log_lambda[-(1:n_z)])
      
      ### Sample beta & eta
      var_inv = as.matrix(bdiag(coef_inv, Sigma_beta_inv, (1 / tau) * I_T))
      tmp_chol = chol(t(M_z) %*% (omega1 * M_z) + t(M_N) %*% (omega2 * M_N) +  var_inv)
      
      if(Poisson){ 
        
        # Calculate repeat quantities for efficiency 
        mh2 <-
          # Add log(Likelihood) Geometric 
          sum(NB_polson_II(N, 1, log_lambda[-(1:n_z)])) +
          sum(dbinom(z, 1, logit_inv(log_lambda[(1:n_z)]), TRUE)) -
          # Subtract log(Likelihood) Poisson 
          sum(dpois(N, exp(log_lambda[-(1:n_z)]), TRUE)) -
          sum(dbinom(z, 1, prob=1-exp(-exp(log_lambda[(1:n_z)])), TRUE)) 
        
          
        coef = backsolve(tmp_chol, 
                         backsolve(tmp_chol,  
                                   t(M_z) %*% kappa + t(M_N) %*% (N-nu)/2 + mean_vec,
                                   transpose=TRUE) + rnorm(p_X + T + 1)) 
          
        log_lambda_star = M %*% coef
          
        mh1 <-
          # Add log(Likelihood) Poisson 
          sum(dpois(N, exp(log_lambda_star[-(1:n_z)]), TRUE)) +
          sum(dbinom(z, 1, prob=1-exp(-exp(log_lambda_star[(1:n_z)])), TRUE)) -
          # Subtract log(Likelihood) Geometric 
          sum(NB_polson_II(N, 1, log_lambda_star[-(1:n_z)])) -
          sum(dbinom(z, 1, logit_inv(log_lambda_star[(1:n_z)]), TRUE)) 
          
        # MH ratio
        mh = exp(mh1 + mh2)
          
        # Accept or keep looping
        if (runif(1) < mh) {
          log_lambda = log_lambda_star
          coef_save[q] = coef[1]
          coef = coef[-1]
          beta = coef[1:p_X]
          xi = coef[(p_X + 1):(p_X + T)]
        }
        
      } else {
        coef = backsolve(tmp_chol, 
                         backsolve(tmp_chol,  
                                   t(M_z) %*% kappa + t(M_N) %*% (N-nu)/2 + mean_vec,
                                   transpose=TRUE) + rnorm(p_X + T + 1)) 
        log_lambda = M %*% coef
        coef_save[q] = coef[1]
        coef = coef[-1]
        beta = coef[1:p_X]
        xi = coef[(p_X + 1):(p_X + T)]
      }

      
      ### Sample tau
      r_tilde = t(xi) %*% xi / 2 + var_rate
      tau = 1 / rgamma(1, shape = q_xi, scale = 1 / r_tilde)
      
      ###
      ### Sample N (refine Stage I samples)
      ###
      
      for(i in 1:n_N){
        
        # Proposal value
        N_star = N_I[i,sample(n_samp, 1, TRUE)]
        
        # Metropolis-Hastings ratio
        if(Poisson){
          mh1 = dpois(N_star, exp(log_lambda[i+n_z]), TRUE)
          mh2 = dpois(N[i], exp(log_lambda[i+n_z]), TRUE)
        } else {
          mh1 = NB_polson_N(N_star, nu, log_lambda[i+n_z])
          mh2 = NB_polson_N(N[i], nu, log_lambda[i+n_z])
        }
        mh = exp(mh1 - mh2)
        
        # Accept or reject
        if (runif(1) < mh) {
          N[i] = N_star
          N_acc[i] = N_acc[i] + 1
        }
        
      }
      
      ### Save Samples
      log_lambda_save[, q] = log_lambda
      N_save[, q] = N
      beta_save[, q] = beta
      xi_save[, q] = xi
      tau_save[q] = tau

      ### Timer
      ### Timer
      if(Poisson){
        if (q %% 100 == 0) cat(q, " ")
      } else {
        if (q %% 1000 == 0) cat(q, " ")
      }
    }
    
    
    list(
      ## MCMC samples
      coef=coef_save,
      log_lambda=log_lambda_save,
      N=N_save,
      beta=beta_save,
      xi=xi_save,
      tau=tau_save
    )
  }

###
### Brook Trout (depletion sampling + presence/absence)
###

## Pólya-gamma (presence/absence only)
PG_binary <- 
  function(# Data
    z,
    X, D, A, S, 
    year_mat, watershed_mat, site_mat,
    # Priors 
    mu_beta, Sigma_beta,
    var_shape=0.001, var_rate=0.001,
    phi_shape=0.1, phi_rate=0.0001,
    # Hyperparameters 
    c_1, n_mcmc){
    
    ###
    ### Packages
    ###
    
    library(BayesLogit)
    library(Boom)
    library(Matrix)
    library(purrr)
    
    ###
    ### Loop Variables
    ###
    
    ### Dimensions 
    n_z = length(z)
    n_y = nrow(X) - n_z
    p_X = ncol(X)
    n_zeta = ncol(watershed_mat)
    T = ncol(year_mat)
    n = ncol(site_mat)
    n_q = length(D)
    
    ### Priors 
    Sigma_beta_inv = solve(Sigma_beta)
    Sigma_beta_inv_times_mu = Sigma_beta_inv %*% mu_beta
    ### How we effectively set the coefficient on area (S) to 1 
    coef_precision = 1e5
    
    ### Objects 
    kappa = z - rep(1, n_z) / 2
    I = diag(1, n_z)
    zeta_COV_inv = (diag(rowSums(A)) - 0.999 * A)
    zeta_COV = solve(zeta_COV_inv)
    q_zeta = n_zeta / 2 + var_shape
    q_xi = T / 2 + var_shape
    q_eta = map_dbl(D, ncol) / 2 + var_shape
    M = cbind(log(S), X, watershed_mat, year_mat, site_mat)
    M_z = M[1:n_z,]
    
    ###
    ### Starting Values
    ###
    
    beta = mu_beta
    zeta = numeric(n_zeta)
    xi = numeric(T)
    eta = numeric(n)
    log_lambda = numeric(n_z + n_y)
    nu = 1
    sigma2_zeta = 1
    sigma2_xi = 1
    sigma2_eta = 1
    phi = rep(100, n_q)
    eta_COV = map2(D, phi, ~exp(-.x/.y))
    eta_COV_chol = chol(bdiag(eta_COV))
    eta_COV_inv = chol2inv(eta_COV_chol)
    rho = 0.5
    xi_COV = make_ar1_corr(T, rho)
    xi_COV_chol = chol(xi_COV)
    xi_COV_inv = chol2inv(xi_COV_chol)
    mean_vec = c(coef_precision, Sigma_beta_inv_times_mu, rep(0, n_zeta + T + n))
    coef_inv = diag(coef_precision, 1)
    
    ### Get eta index
    sizes <- vapply(eta_COV, ncol, integer(1))  # block sizes
    starts <- c(1, cumsum(sizes[-length(sizes)]) + 1)
    ends <- cumsum(sizes)
    # Return list of index vectors
    eta_idx <- Map(seq, starts, ends)

    ###
    ### Acceptance
    ###
    
    rho_acc = 0
    rho_tune = 2.4^2 
    phi_acc = numeric(n_q)
    phi_tune = rep(2.4^2, n_q)
    ### Hyperparameters for tuning adaptive sampler 
    r_opt = 0.234
    
    ###
    ### Storage
    ###
    
    beta_save = matrix(NA, p_X, n_mcmc)
    zeta_save = matrix(NA, n_zeta, n_mcmc)
    xi_save = matrix(NA, T, n_mcmc)
    eta_save = matrix(NA, n, n_mcmc)
    log_lambda_save = matrix(NA, n_z+n_y, n_mcmc)
    sigma2_zeta_save = numeric(n_mcmc)
    sigma2_eta_save = matrix(NA, n_q, n_mcmc)
    sigma2_xi_save = numeric(n_mcmc)
    rho_save = numeric(n_mcmc)
    phi_save = matrix(NA, n_q, n_mcmc)
    coef_save = numeric(n_mcmc)
    
    ###
    ### MCMC loop
    ###
    
    for (q in 1:n_mcmc) {
      
      ### For automatic tuning 
      gamma = 1.0 / q^c_1
      
      ### Sample omega
      omega = rpg(n_z, 1, log_lambda[1:n_z])
      
      ### Sample beta, zeta, xi
      var_inv = as.matrix(bdiag(coef_inv, Sigma_beta_inv, (1 / sigma2_zeta) * zeta_COV_inv, (1 / sigma2_xi) * xi_COV_inv, eta_COV_inv))
      # eta_COV_inv already has sigma2_eta included
      tmp_chol = chol(t(M_z) %*% (omega * M_z) +  var_inv)
      coef = backsolve(tmp_chol, 
                       backsolve(tmp_chol,  
                                 t(M_z) %*% kappa + mean_vec,
                                 transpose=TRUE) + rnorm(p_X + n_zeta + T + n + 1)) 
      log_lambda = M %*% coef
      coef_save[q] = coef[1]
      coef = coef[-1]
      beta = coef[1:p_X]
      zeta = coef[(p_X + 1):(p_X + n_zeta)]
      xi = coef[(p_X + n_zeta + 1):(p_X + n_zeta + T)]
      eta = coef[(p_X + n_zeta + T + 1):(p_X + n_zeta + T + n)]
      
      ### Sample sigma2_w
      r_tilde = t(zeta) %*% zeta_COV_inv %*% zeta / 2 + var_rate
      sigma2_zeta = 1 / rgamma(1, shape = q_zeta, scale = 1 / r_tilde)
      
      ### Sample sigma2_t
      r_tilde = t(xi) %*% xi_COV_inv %*% xi / 2 + var_rate
      sigma2_xi = 1 / rgamma(1, shape = q_xi, scale = 1 / r_tilde)
      
      ### Sample sigma2_d
      r_tilde = map2_dbl(eta_idx, sigma2_eta, ~sum(t(eta[.x]) %*% (.y * eta_COV_inv[.x,.x]) %*% eta[.x] / 2 + var_rate)) 
      sigma2_eta = 1 / rgamma(n_q, shape = q_eta, scale = 1 / r_tilde)
      
      ###
      ### Sample rho
      ###
      
      # random walk proposals
      rho_star = logit_inv(rnorm(1, logit(rho), sd=sqrt(rho_tune)))
      
      if(0.01 < rho_star & rho_star < 0.99){
        
        # Calculate AR(1) correlation matrix 
        xi_COV_star = sigma2_xi * make_ar1_corr(T, rho_star)
        xi_COV_chol_star = chol(xi_COV_star)
        xi_COV_chol = chol(sigma2_xi * make_ar1_corr(T, rho))
        
        # Under proposal 
        mh1 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol_star, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol_star))) 
        
        # Under current values
        mh2 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol))) 
        
        # MH ratio
        mh = exp(mh1-mh2)
        
        # Accept or reject
        if (runif(1) < mh) {
          rho = rho_star
          xi_COV = xi_COV_star
          xi_COV_inv = chol2inv(xi_COV_chol_star)
          rho_acc = rho_acc + 1
          rho_tune = exp(log(rho_tune) + gamma*(1.0 - r_opt));
        } else { 
          rho_tune = exp(log(rho_tune) + gamma*(0.0 - r_opt));
        } 
        
      } else {
        rho_tune = exp(log(rho_tune) + gamma*(0.0 - r_opt));
      }
      
      ###
      ### Sample phi
      ###
      
      for (i in 1:n_q) {
        
        # Random walk proposal for network i
        phi_star <- exp(rnorm(1, log(phi[i]), sd = sqrt(phi_tune[i])))
        
        # Calculate correlation matrix for this network
        eta_COV_star_i <- exp(-D[[i]] / phi_star)
        
        # Attempt Cholesky
        chol_success <- TRUE
        eta_COV_chol_star_i <- tryCatch({
          suppressWarnings(chol(sigma2_eta[i] * eta_COV_star_i))
        }, error = function(e) {
          chol_success <<- FALSE
          NULL
        })
        
        # Only compute MH ratio if Cholesky succeeded
        if (chol_success & phi_star > 10) {
          
          # Current correlation matrix Cholesky
          eta_COV_chol_i <- chol(sigma2_eta[i] * eta_COV[[i]])
          eta_tmp <- eta[eta_idx[[i]]]
          
          # Likelihood under proposal
          mh1 <- (-0.5) * sum(backsolve(eta_COV_chol_star_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_star_i))) +
            # Prior on phi 
            dgamma(phi_star, shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # Likelihood under current
          mh2 <- (-0.5) * sum(backsolve(eta_COV_chol_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_i))) +
            # Prior on phi 
            dgamma(phi[i], shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # MH ratio
          mh <- exp(mh1 - mh2)
          
          # Accept or reject
          if (runif(1) < mh) {
            phi[i] <- phi_star
            eta_COV[[i]] <- eta_COV_star_i
            phi_acc[i] <- phi_acc[i] + 1
            phi_tune[i] <- exp(log(phi_tune[i]) + gamma*(1.0 - r_opt))
          } else {
            phi_tune[i] <- exp(log(phi_tune[i]) + gamma*(0.0 - r_opt))
          }
          
        } else {
          # Skip update if Cholesky failed
          phi_tune[i] <- exp(log(phi_tune[i]) + gamma*(0.0 - r_opt))
        }
      }
      
      ### Update eta_COV_inv 
      eta_COV_inv <- bdiag(map2(eta_COV, sigma2_eta, ~solve(.y * .x)))
      
      
      ### Save Samples
      log_lambda_save[, q] = log_lambda
      beta_save[, q] = beta
      zeta_save[, q] = zeta
      xi_save[, q] = xi
      eta_save[, q] = eta
      sigma2_zeta_save[q] = sigma2_zeta
      sigma2_eta_save[, q] = sigma2_eta
      sigma2_xi_save[q] = sigma2_xi
      rho_save[q] = rho
      phi_save[, q] = phi
      
      ### Timer
      if (q %% 1000 == 0) cat(q, " ")
    }
    
    
    list(
      ## MCMC samples
      coef=coef_save,
      log_lambda=log_lambda_save,
      beta=beta_save,
      zeta=zeta_save,
      xi=xi_save,
      eta=eta_save,
      sigma2_zeta=sigma2_zeta_save,
      sigma2_eta=sigma2_eta_save,
      sigma2_xi=sigma2_xi_save,
      rho=rho_save,
      phi=phi_save,
      ## Acceptance 
      rho_tune=rho_tune,
      rho_acc=rho_acc,
      phi_tune=phi_tune,
      phi_acc=phi_acc
    )
  }

## Pólya-gamma (count only)
PG_count <- 
  function(# Data
    y, J_y, y_idx, y_sum, N_deplete,
    H, X, D, A, S, 
    year_mat, watershed_mat, site_mat,
    # Priors 
    mu_gamma, Sigma_gamma,
    mu_beta, Sigma_beta,
    var_shape=0.001, var_rate=0.001,
    phi_shape=0.1, phi_rate=0.0001,
    # Hyperparameters 
    c_1, n_mcmc){
    
    ###
    ### Packages
    ###
    
    library(BayesLogit)
    library(Boom)
    library(Matrix)
    library(purrr)
    
    ###
    ### Loop Variables
    ###
    
    ### Dimensions 
    n_N = length(y_idx)
    n_v = nrow(X) - n_N
    n_y = length(y)
    p_X = ncol(X)
    p_H = ncol(H)
    n_zeta = ncol(watershed_mat)
    T = ncol(year_mat)
    n = ncol(site_mat)
    n_q = length(D)
    
    ### Priors 
    Sigma_beta_inv = solve(Sigma_beta)
    Sigma_beta_inv_times_mu = Sigma_beta_inv %*% mu_beta
    Sigma_gamma_inv = solve(Sigma_gamma)
    Sigma_gamma_inv_times_mu = Sigma_gamma_inv %*% mu_gamma
    ### How we effectively set the coefficient on area (S) to 1 
    coef_precision = 1e5
    
    ### Objects 
    zeta_COV_inv = (diag(rowSums(A)) - 0.999 * A)
    zeta_COV = solve(zeta_COV_inv)
    q_zeta = n_zeta / 2 + var_shape
    q_xi = T / 2 + var_shape
    q_eta = map_dbl(D, ncol) / 2 + var_shape
    M = cbind(log(S), X, watershed_mat, year_mat, site_mat)
    M_v = M[1:n_v,]
    M_N = M[-(1:n_v),]
    H_big <- H[rep(seq_len(nrow(H)), times = J_y), ]
    
    ###
    ### Starting Values
    ###
    
    beta = mu_beta
    gamma = mu_gamma
    pi_logit = H %*% gamma
    pi = logit_inv(pi_logit)
    zeta = numeric(n_zeta)
    xi = numeric(T)
    eta = numeric(n)
    log_lambda = numeric(n_v+n_N)
    N = y_sum
    nu = 1
    sigma2_zeta = 1
    sigma2_xi = 1
    sigma2_eta = 1
    phi = rep(100, n_q)
    eta_COV = map2(D, phi, ~exp(-.x/.y))
    eta_COV_chol = chol(bdiag(eta_COV))
    eta_COV_inv = chol2inv(eta_COV_chol)
    rho = 0.5
    xi_COV = make_ar1_corr(T, rho)
    xi_COV_chol = chol(xi_COV)
    xi_COV_inv = chol2inv(xi_COV_chol)
    mean_vec = c(coef_precision, Sigma_beta_inv_times_mu, rep(0, n_zeta + T + n))
    coef_inv = diag(coef_precision, 1)
    
    ### Get eta index
    sizes <- vapply(eta_COV, ncol, integer(1))  # block sizes
    starts <- c(1, cumsum(sizes[-length(sizes)]) + 1)
    ends <- cumsum(sizes)
    # Return list of index vectors
    eta_idx <- Map(seq, starts, ends)
    
    ###
    ### Acceptance
    ###
    
    rho_acc = 0
    phi_acc = numeric(n_q)
    N_acc = numeric(n_N)
    rho_tune = 2.4^2 
    phi_tune = rep(2.4^2 , n_q)
    ### Hyperparameters for tuning adaptive sampler 
    r_opt = 0.234
    
    ###
    ### Storage
    ###
    
    beta_save = matrix(NA, p_X, n_mcmc)
    gamma_save = matrix(NA, p_H, n_mcmc)
    N_save = matrix(NA, n_N, n_mcmc)
    zeta_save = matrix(NA, n_zeta, n_mcmc)
    xi_save = matrix(NA, T, n_mcmc)
    eta_save = matrix(NA, n, n_mcmc)
    log_lambda_save = matrix(NA, n_N+n_v, n_mcmc)
    sigma2_zeta_save = numeric(n_mcmc)
    sigma2_eta_save = matrix(NA, n_q, n_mcmc)
    sigma2_xi_save = numeric(n_mcmc)
    rho_save = numeric(n_mcmc)
    phi_save = matrix(NA, n_q, n_mcmc)
    coef_save = numeric(n_mcmc)
    
    ###
    ### MCMC loop
    ###
    
    for (q in 1:n_mcmc) {
      
      ### For automatic tuning 
      tune = 1.0 / q^c_1
      
      ###
      ### Sample fixed and random effects 
      ###
      
      ### Sample omega
      omega = rpg(n_N, N+nu, log_lambda[-(1:n_v)])
      
      ### Sample beta, zeta, xi, and eta
      var_inv = as.matrix(bdiag(coef_inv, Sigma_beta_inv, (1 / sigma2_zeta) * zeta_COV_inv, (1 / sigma2_xi) * xi_COV_inv, eta_COV_inv))
      tmp_chol = chol(t(M_N) %*% (omega * M_N) +  var_inv)
      coef = backsolve(tmp_chol, 
                       backsolve(tmp_chol,  
                                 t(M_N) %*% (N-nu)/2 + mean_vec,
                                 transpose=TRUE) + rnorm(p_X + n_zeta + T + n + 1)) 
      log_lambda = M %*% coef
      coef_save[q] = coef[1]
      coef = coef[-1]
      beta = coef[1:p_X]
      zeta = coef[(p_X + 1):(p_X + n_zeta)]
      xi = coef[(p_X + n_zeta + 1):(p_X + n_zeta + T)]
      eta = coef[(p_X + n_zeta + T + 1):(p_X + n_zeta + T + n)]
      
      ### Sample sigma2_zeta
      r_tilde = t(zeta) %*% zeta_COV_inv %*% zeta / 2 + var_rate
      sigma2_zeta = 1 / rgamma(1, shape = q_zeta, scale = 1 / r_tilde)
      
      ### Sample sigma2_xi
      r_tilde = t(xi) %*% xi_COV_inv %*% xi / 2 + var_rate
      sigma2_xi = 1 / rgamma(1, shape = q_xi, scale = 1 / r_tilde)
      
      ### Sample sigma2_eta
      r_tilde = map2_dbl(eta_idx, sigma2_eta, ~sum(t(eta[.x]) %*% (.y * eta_COV_inv[.x,.x]) %*% eta[.x] / 2 + var_rate)) 
      sigma2_eta = 1 / rgamma(n_q, shape = q_eta, scale = 1 / r_tilde)
      
      ###
      ### Sample rho
      ###
      
      # random walk proposals
      rho_star = logit_inv(rnorm(1, logit(rho), sd=sqrt(rho_tune)))
      
      if(0.01 < rho_star & rho_star < 0.99){
        
        # Calculate AR(1) correlation matrix 
        xi_COV_star = sigma2_xi * make_ar1_corr(T, rho_star)
        xi_COV_chol_star = chol(xi_COV_star)
        xi_COV_chol = chol(sigma2_xi * make_ar1_corr(T, rho))
        
        # Under proposal 
        mh1 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol_star, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol_star))) 
        
        # Under current values
        mh2 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol))) 
        
        # MH ratio
        mh = exp(mh1-mh2)
        
        # Accept or reject
        if (runif(1) < mh) {
          rho = rho_star
          xi_COV = xi_COV_star
          xi_COV_inv = chol2inv(xi_COV_chol_star)
          rho_acc = rho_acc + 1
          rho_tune = exp(log(rho_tune) + tune*(1.0 - r_opt))
        } else { 
          rho_tune = exp(log(rho_tune) + tune*(0.0 - r_opt))
        } 
        
      } else {
        rho_tune = exp(log(rho_tune) + tune*(0.0 - r_opt))
      }
      
      ###
      ### Sample N
      ###
      
      for(i in 1:n_N){
        
        # Proposal value
        N_star = rpois(1, N[i])
        
        if(N_star >= y_sum[i]){
          
          # Metropolis-Hastings ratio
          mh1 = NB_polson_N(N_star, nu, log_lambda[i+n_v]) +
            dpois(N[i], N_star, TRUE)
          mh2 = NB_polson_N(N[i], nu, log_lambda[i+n_v]) +
            dpois(N_star, N[i], TRUE)
          
          for(j in 1:J_y[i]){
            mh1 = dbinom(y[y_idx[i]+j], N_star-N_deplete[y_idx[i]+j], pi[i], TRUE) + mh1
            mh2 = dbinom(y[y_idx[i]+j], N[i]-N_deplete[y_idx[i]+j], pi[i], TRUE) + mh2
          }
          
          mh = exp(mh1 - mh2)
          
          # Accept or reject
          if (runif(1) < mh) {
            N[i] = N_star
            N_acc[i] = N_acc[i] + 1
          }
          
        }
        
      }
      N_big = rep(N, J_y) - N_deplete
      non_zero_idx = which(N_big!=0)
      
      ###
      ### Sample gamma
      ###
      
      ### Sample omega
      omega = rpg(length(non_zero_idx), N_big[non_zero_idx], rep(pi_logit, J_y)[non_zero_idx])
      
      ### Sample gamma
      tmp_chol = chol(t(H_big[non_zero_idx,]) %*% (omega * H_big[non_zero_idx,]) + Sigma_gamma_inv)
      gamma = backsolve(tmp_chol, 
                        backsolve(tmp_chol,  
                                  t(H_big[non_zero_idx,]) %*% (y[non_zero_idx]-N_big[non_zero_idx]/2),
                                  transpose=TRUE) + rnorm(p_H)) 
      pi_logit = H %*% gamma
      pi = logit_inv(pi_logit)
      
      ###
      ### Sample phi
      ###
      
      for (i in 1:n_q) {
        
        # Random walk proposal for network i
        phi_star <- exp(rnorm(1, log(phi[i]), sd = sqrt(phi_tune[i])))
        
        # Calculate correlation matrix for this network
        eta_COV_star_i <- exp(-D[[i]] / phi_star)
        
        # Attempt Cholesky
        chol_success <- TRUE
        eta_COV_chol_star_i <- tryCatch({
          suppressWarnings(chol(sigma2_eta[i] * eta_COV_star_i))
        }, error = function(e) {
          chol_success <<- FALSE
          NULL
        })
        
        # Only compute MH ratio if Cholesky succeeded
        if (chol_success & phi_star > 10) {
          
          # Current correlation matrix Cholesky
          eta_COV_chol_i <- chol(sigma2_eta[i] * eta_COV[[i]])
          eta_tmp <- eta[eta_idx[[i]]]
          
          # Likelihood under proposal
          mh1 <- (-0.5) * sum(backsolve(eta_COV_chol_star_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_star_i))) +
            # Prior on phi 
            dgamma(phi_star, shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # Likelihood under current
          mh2 <- (-0.5) * sum(backsolve(eta_COV_chol_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_i))) +
            # Prior on phi 
            dgamma(phi[i], shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # MH ratio
          mh <- exp(mh1 - mh2)
          
          # Accept or reject
          if (runif(1) < mh) {
            phi[i] <- phi_star
            eta_COV[[i]] <- eta_COV_star_i
            phi_acc[i] <- phi_acc[i] + 1
            phi_tune[i] <- exp(log(phi_tune[i]) + tune*(1.0 - r_opt))
          } else {
            phi_tune[i] <- exp(log(phi_tune[i]) + tune*(0.0 - r_opt))
          }
          
        } else {
          # Skip update if Cholesky failed
          phi_tune[i] <- exp(log(phi_tune[i]) + tune*(0.0 - r_opt))
        }
      }
      
      ### Update eta_COV_inv 
      eta_COV_inv <- bdiag(map2(eta_COV, sigma2_eta, ~solve(.y * .x)))

      ### Save Samples
      log_lambda_save[, q] = log_lambda
      N_save[, q] = N
      gamma_save[, q] = gamma
      beta_save[, q] = beta
      zeta_save[, q] = zeta
      xi_save[, q] = xi
      eta_save[, q] = eta
      sigma2_zeta_save[q] = sigma2_zeta
      sigma2_eta_save[, q] = sigma2_eta
      sigma2_xi_save[q] = sigma2_xi
      rho_save[q] = rho
      phi_save[, q] = phi
      
      ### Timer
      if (q %% 1000 == 0) cat(q, " ")
    }
    
    
    list(
      ## MCMC samples
      coef=coef_save,
      log_lambda=log_lambda_save,
      N=N_save,
      gamma=gamma_save,
      beta=beta_save,
      zeta=zeta_save,
      xi=xi_save,
      eta=eta_save,
      sigma2_zeta=sigma2_zeta_save,
      sigma2_eta=sigma2_eta_save,
      sigma2_xi=sigma2_xi_save,
      rho=rho_save,
      phi=phi_save,
      ## Acceptance 
      rho_tune=rho_tune,
      rho_acc=rho_acc,
      phi_tune=phi_tune,
      phi_acc=phi_acc,
      N_acc=N_acc
    )
  }

## Pólya-gamma (count + presence/absence)
PG_joint <- 
  function(# Data
    z, y, J_y, y_idx, y_sum, N_deplete,
    H, X, D, A, S, 
    year_mat, watershed_mat, site_mat,
    # Priors 
    mu_gamma, Sigma_gamma,
    mu_beta, Sigma_beta,
    var_shape=0.001, var_rate=0.001,
    phi_shape=0.1, phi_rate=0.0001,
    # Hyperparameters 
    c_1, n_mcmc){
    
    ###
    ### Packages
    ###
    
    library(BayesLogit)
    library(Boom)
    library(Matrix)
    library(purrr)
    
    ###
    ### Loop Variables
    ###
    
    ### Dimensions 
    n_z = length(z)
    n_N = length(y_idx)
    n_y = length(y)
    p_X = ncol(X)
    p_H = ncol(H)
    n_zeta = ncol(watershed_mat)
    T = ncol(year_mat)
    n = ncol(site_mat)
    n_q = length(D)
    
    ### Priors 
    Sigma_beta_inv = solve(Sigma_beta)
    Sigma_beta_inv_times_mu = Sigma_beta_inv %*% mu_beta
    Sigma_gamma_inv = solve(Sigma_gamma)
    Sigma_gamma_inv_times_mu = Sigma_gamma_inv %*% mu_gamma
    ### How we effectively set the coefficient on area (S) to 1 
    coef_precision = 1e5
    
    ### Objects 
    kappa = z - rep(1, n_z) / 2
    I = diag(1, n_z)
    zeta_COV_inv = (diag(rowSums(A)) - 0.999 * A)
    zeta_COV = solve(zeta_COV_inv)
    q_zeta = n_zeta / 2 + var_shape
    q_xi = T / 2 + var_shape
    q_eta = map_dbl(D, ncol) / 2 + var_shape
    M = cbind(log(S), X, watershed_mat, year_mat, site_mat)
    M_z = M[1:n_z,]
    M_N = M[-(1:n_z),]
    H_big <- H[rep(seq_len(nrow(H)), times = J_y), ]
    
    ###
    ### Starting Values
    ###
    
    beta = mu_beta
    gamma = mu_gamma
    pi_logit = H %*% gamma
    pi = logit_inv(pi_logit)
    zeta = numeric(n_zeta)
    xi = numeric(T)
    eta = numeric(n)
    log_lambda = numeric(n_z+n_y)
    N = y_sum
    nu = 1
    sigma2_zeta = 1
    sigma2_xi = 1
    sigma2_eta = 1
    phi = rep(100, n_q)
    eta_COV = map2(D, phi, ~exp(-.x/.y))
    eta_COV_chol = chol(bdiag(eta_COV))
    eta_COV_inv = chol2inv(eta_COV_chol)
    rho = 0.5
    xi_COV = make_ar1_corr(T, rho)
    xi_COV_chol = chol(xi_COV)
    xi_COV_inv = chol2inv(xi_COV_chol)
    mean_vec = c(coef_precision, Sigma_beta_inv_times_mu, rep(0, n_zeta + T + n))
    coef_inv = diag(coef_precision, 1)
    
    ### Get eta index
    sizes <- vapply(eta_COV, ncol, integer(1))  # block sizes
    starts <- c(1, cumsum(sizes[-length(sizes)]) + 1)
    ends <- cumsum(sizes)
    # Return list of index vectors
    eta_idx <- Map(seq, starts, ends)
    
    ###
    ### Acceptance
    ###
    
    rho_acc = 0
    phi_acc = numeric(n_q)
    N_acc = numeric(n_N)
    rho_tune = 2.4^2 
    phi_tune = rep(2.4^2, n_q)
    ### Hyperparameters for tuning adaptive sampler 
    r_opt = 0.234
    
    ###
    ### Storage
    ###
    
    beta_save = matrix(NA, p_X, n_mcmc)
    gamma_save = matrix(NA, p_H, n_mcmc)
    N_save = matrix(NA, n_N, n_mcmc)
    zeta_save = matrix(NA, n_zeta, n_mcmc)
    xi_save = matrix(NA, T, n_mcmc)
    eta_save = matrix(NA, n, n_mcmc)
    log_lambda_save = matrix(NA, n_z+n_N, n_mcmc)
    sigma2_zeta_save = numeric(n_mcmc)
    sigma2_eta_save = matrix(NA, n_q, n_mcmc)
    sigma2_xi_save = numeric(n_mcmc)
    rho_save = numeric(n_mcmc)
    phi_save = matrix(NA, n_q, n_mcmc)
    coef_save = numeric(n_mcmc)
    
    ###
    ### MCMC loop
    ###
    
    for (q in 1:n_mcmc) {
      
      ### For automatic tuning 
      tune = 1.0 / q^c_1
      
      ###
      ### Sample fixed and random effects 
      ###
      
      ### Sample omegas
      omega1 = rpg(n_z, 1, log_lambda[1:n_z])
      omega2 = rpg(n_N, N+nu, log_lambda[-(1:n_z)])
      
      ### Sample beta, zeta, xi, and eta
      var_inv = as.matrix(bdiag(coef_inv, Sigma_beta_inv, (1 / sigma2_zeta) * zeta_COV_inv, (1 / sigma2_xi) * xi_COV_inv, eta_COV_inv))
      tmp_chol = chol(t(M_z) %*% (omega1 * M_z) + t(M_N) %*% (omega2 * M_N) +  var_inv)
      coef = backsolve(tmp_chol, 
                       backsolve(tmp_chol,  
                                 t(M_z) %*% kappa + t(M_N) %*% (N-nu)/2 + mean_vec,
                                 transpose=TRUE) + rnorm(p_X + n_zeta + T + n + 1)) 
      log_lambda = M %*% coef
      coef_save[q] = coef[1]
      coef = coef[-1]
      beta = coef[1:p_X]
      zeta = coef[(p_X + 1):(p_X + n_zeta)]
      xi = coef[(p_X + n_zeta + 1):(p_X + n_zeta + T)]
      eta = coef[(p_X + n_zeta + T + 1):(p_X + n_zeta + T + n)]
      
      ### Sample sigma2_w
      r_tilde = t(zeta) %*% zeta_COV_inv %*% zeta / 2 + var_rate
      sigma2_zeta = 1 / rgamma(1, shape = q_zeta, scale = 1 / r_tilde)
      
      ### Sample sigma2_t
      r_tilde = t(xi) %*% xi_COV_inv %*% xi / 2 + var_rate
      sigma2_xi = 1 / rgamma(1, shape = q_xi, scale = 1 / r_tilde)
      
      ### Sample sigma2_eta
      r_tilde = map2_dbl(eta_idx, sigma2_eta, ~sum(t(eta[.x]) %*% (.y * eta_COV_inv[.x,.x]) %*% eta[.x] / 2 + var_rate)) 
      sigma2_eta = 1 / rgamma(n_q, shape = q_eta, scale = 1 / r_tilde)
      
      ###
      ### Sample rho
      ###
      
      # random walk proposals
      rho_star = logit_inv(rnorm(1, logit(rho), sd=sqrt(rho_tune)))
      
      if(0.01 < rho_star & rho_star < 0.99){
        
        # Calculate AR(1) correlation matrix 
        xi_COV_star = sigma2_xi * make_ar1_corr(T, rho_star)
        xi_COV_chol_star = chol(xi_COV_star)
        xi_COV_chol = chol(sigma2_xi * make_ar1_corr(T, rho))
        
        # Under proposal 
        mh1 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol_star, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol_star))) 
        
        # Under current values
        mh2 <-
          # Likelihood
          (-0.5) * sum(backsolve(xi_COV_chol, xi, transpose = TRUE)^2) + 
          (-0.5) * 2 * sum(log(diag(xi_COV_chol))) 
        
        # MH ratio
        mh = exp(mh1-mh2)
        
        # Accept or reject
        if (runif(1) < mh) {
          rho = rho_star
          xi_COV = xi_COV_star
          xi_COV_inv = chol2inv(xi_COV_chol_star)
          rho_acc = rho_acc + 1
          rho_tune = exp(log(rho_tune) + tune*(1.0 - r_opt));
        } else { 
          rho_tune = exp(log(rho_tune) + tune*(0.0 - r_opt));
        } 
        
      } else {
        rho_tune = exp(log(rho_tune) + tune*(0.0 - r_opt));
      }
      
      ###
      ### Sample N
      ###
      
      for(i in 1:n_N){
        
        # Proposal value
        N_star = rpois(1, N[i])
        
        if(N_star >= y_sum[i]){
          
          # Metropolis-Hastings ratio
          mh1 = NB_polson_N(N_star, nu, log_lambda[i+n_z]) +
            dpois(N[i], N_star, TRUE)
          mh2 = NB_polson_N(N[i], nu, log_lambda[i+n_z]) +
            dpois(N_star, N[i], TRUE)
          
          for(j in 1:J_y[i]){
            mh1 = dbinom(y[y_idx[i]+j], N_star-N_deplete[y_idx[i]+j], pi[i], TRUE) + mh1
            mh2 = dbinom(y[y_idx[i]+j], N[i]-N_deplete[y_idx[i]+j], pi[i], TRUE) + mh2
          }
          
          mh = exp(mh1 - mh2)
          
          # Accept or reject
          if (runif(1) < mh) {
            N[i] = N_star
            N_acc[i] = N_acc[i] + 1
          }
          
        }
        
      }
      N_big = rep(N, J_y) - N_deplete
      non_zero_idx = which(N_big!=0)
      
      ###
      ### Sample gamma
      ###
      
      ### Sample omega
      omega = rpg(length(non_zero_idx), N_big[non_zero_idx], rep(pi_logit, J_y)[non_zero_idx])
      
      ### Sample gamma
      tmp_chol = chol(t(H_big[non_zero_idx,]) %*% (omega * H_big[non_zero_idx,]) + Sigma_gamma_inv)
      gamma = backsolve(tmp_chol, 
                        backsolve(tmp_chol,  
                                  t(H_big[non_zero_idx,]) %*% (y[non_zero_idx]-N_big[non_zero_idx]/2),
                                  transpose=TRUE) + rnorm(p_H)) 
      pi_logit = H %*% gamma
      pi = logit_inv(pi_logit)
      
      
      ###
      ### Sample phi
      ###
      
      for (i in 1:n_q) {
        
        # Random walk proposal for network i
        phi_star <- exp(rnorm(1, log(phi[i]), sd = sqrt(phi_tune[i])))
        
        # Calculate correlation matrix for this network
        eta_COV_star_i <- exp(-D[[i]] / phi_star)
        
        # Attempt Cholesky
        chol_success <- TRUE
        eta_COV_chol_star_i <- tryCatch({
          suppressWarnings(chol(sigma2_eta[i] * eta_COV_star_i))
        }, error = function(e) {
          chol_success <<- FALSE
          NULL
        })
        
        # Only compute MH ratio if Cholesky succeeded
        if (chol_success & phi_star > 10) {
          
          # Current correlation matrix Cholesky
          eta_COV_chol_i <- chol(sigma2_eta[i] * eta_COV[[i]])
          eta_tmp <- eta[eta_idx[[i]]]
          
          # Likelihood under proposal
          mh1 <- (-0.5) * sum(backsolve(eta_COV_chol_star_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_star_i))) +
            # Prior on phi 
            dgamma(phi_star, shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # Likelihood under current
          mh2 <- (-0.5) * sum(backsolve(eta_COV_chol_i, eta_tmp, transpose = TRUE)^2) +
                 (-0.5) * 2 * sum(log(diag(eta_COV_chol_i))) +
            # Prior on phi 
            dgamma(phi[i], shape=phi_shape, rate=phi_rate, log=TRUE)
          
          # MH ratio
          mh <- exp(mh1 - mh2)
          
          # Accept or reject
          if (runif(1) < mh) {
            phi[i] <- phi_star
            eta_COV[[i]] <- eta_COV_star_i
            phi_acc[i] <- phi_acc[i] + 1
            phi_tune[i] <- exp(log(phi_tune[i]) + tune*(1.0 - r_opt))
          } else {
            phi_tune[i] <- exp(log(phi_tune[i]) + tune*(0.0 - r_opt))
          }
          
        } else {
          # Skip update if Cholesky failed
          phi_tune[i] <- exp(log(phi_tune[i]) + tune*(0.0 - r_opt))
        }
      }
      
      ### Update eta_COV_inv 
      eta_COV_inv <- bdiag(map2(eta_COV, sigma2_eta, ~solve(.y * .x)))

      ### Save Samples
      log_lambda_save[, q] = log_lambda
      N_save[, q] = N
      gamma_save[, q] = gamma
      beta_save[, q] = beta
      zeta_save[, q] = zeta
      xi_save[, q] = xi
      eta_save[, q] = eta
      sigma2_zeta_save[q] = sigma2_zeta
      sigma2_eta_save[, q] = sigma2_eta
      sigma2_xi_save[q] = sigma2_xi
      rho_save[q] = rho
      phi_save[, q] = phi
      
      ### Timer
      if (q %% 1000 == 0) cat(q, " ")
    }
    
    
    list(
      ## MCMC samples
      coef=coef_save,
      log_lambda=log_lambda_save,
      N=N_save,
      gamma=gamma_save,
      beta=beta_save,
      zeta=zeta_save,
      xi=xi_save,
      eta=eta_save,
      sigma2_zeta=sigma2_zeta_save,
      sigma2_eta=sigma2_eta_save,
      sigma2_xi=sigma2_xi_save,
      rho=rho_save,
      phi=phi_save,
      ## Acceptance 
      rho_tune=rho_tune,
      rho_acc=rho_acc,
      phi_tune=phi_tune,
      phi_acc=phi_acc,
      N_acc=N_acc
    )
  }
