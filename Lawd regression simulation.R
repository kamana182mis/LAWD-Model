rm(list=ls())
library(doParallel)
library(foreach)



####################### Uncensored data #############################

short_mle_analysis <- function(i_range, true_params, lower_bound, upper_bound, batch_size = 250) {
  k_range <- 1000
  num_batches <- ceiling(k_range / batch_size)
  results <- matrix(0, ncol = length(true_params), nrow = k_range, dimnames = list(NULL, names(true_params)))
  
  
  gen_data <- function() {
    ds <- c()  # Uncensored samples
    q_u = 0.9
    for (i in 1:i_range) {
      sol <- function(x) {
        zp <- rnorm(1)
        meu <- exp(true_params['phi1'] + true_params['phi2'] * zp)
        alph <- ((q_u / (1 - q_u))^(1 / true_params['eta']) - 
                   true_params['beta'] * meu^true_params['gamma']) / meu^true_params['theta']
        alph * x^true_params['theta'] + 
          true_params['beta'] * x^true_params['gamma'] - 
          ((u) / (1 - u))^(1 / true_params['eta'])
      }
      u <- runif(1)
      sample <- uniroot(sol, c(0.001, 100), extendInt = "yes")$root
      ds <- c(ds, sample)
    }
    return(list(ds = ds))
  }
  
  log_likelihood <- function(x, data) {
    the <- x[1]; bet <- x[2]; gam <- x[3]; et <- x[4]; ph1 <- x[5]; ph2 <- x[6]
    q_u = 0.9
    ds <- data[[1]]  # Uncensored samples
    sum_ds <- 0
    sum1_ds <- 0
    sum2_ds <- 0
    for (j in 1:length(ds)) {
      z1 <- rnorm(1)
      meu <- exp(ph1 + ph2 * z1)
      alp <- ((q_u / (1 - q_u))^(1 / et) - bet * meu^gam) / meu^the
      if (alp > 0) {
        sum_ds <- sum_ds + log(alp * the * ds[j]^(the - 1) + bet * gam * ds[j]^(gam - 1))
        sum1_ds <- sum1_ds + log(alp * ds[j]^the + bet * ds[j]^gam)
        sum2_ds <- sum2_ds + log(1 + (alp * ds[j]^the + bet * ds[j]^gam)^et)
      } else {
        return(-1e10)  # Penalty
      }
    }
    ff <- length(ds) * log(et) + sum_ds + (et - 1) * sum1_ds - 2 * sum2_ds
    return(ff)
  }
  
  # Set up parallel backend
  num_cores <- detectCores() - 1  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 24) 
  
  # Run the analysis in batches
  for (batch in 1:num_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, k_range)
    batch_results <- foreach(k = batch_start:batch_end, .combine = rbind, .packages = "stats") %dopar% {
      successful <- FALSE
      result <- NULL
      while (!successful) {
        data <- gen_data()
        optim_result <- optim(true_params, log_likelihood, data = data,
                              lower = lower_bound, upper = upper_bound,
                              method = "L-BFGS-B", control = list(fnscale = -1), hessian = TRUE)
        if (optim_result$value > -1e10) {  
          successful <- TRUE
          result <- optim_result$par
        }
      }
      return(result)
    }
    
    results[batch_start:batch_end, ] <- batch_results
    print(paste("Batch", batch, "completed"))
  }
  
  stopCluster(cl)
  
  est_list <- list()
  summary_stats <- sapply(seq_along(true_params), function(i) {
    est <- results[, i]
    est_list[[names(true_params)[i]]] <<- est  
    rmse <- sqrt(mean((est - true_params[i])^2))
    sd <- sd(est)
    rel_bias <- mean(abs(est - true_params[i]) / true_params[i])
    c(mean = mean(est), rmse = rmse, sd = sd, rel_bias = rel_bias)
  })
  
  list(summary_stats = as.data.frame(t(summary_stats), row.names = names(true_params)),
       estimates = est_list)  
}



tp <- c(theta = 1.8, beta = 1.3, gamma = 0.5, eta = 1.1, phi1 = 0.8, phi2 = 0.7)
lb = tp - 0.49
ub = tp + 0.51



results6 <- short_mle_analysis(i_range = 20, tp, lb, ub)
results6$summary_stats

results7 <- short_mle_analysis(i_range = 30, tp, lb, ub)
results7$summary_stats


results0 <- short_mle_analysis(i_range = 50, tp, lb, ub)
results0$summary_stats

results1 <- short_mle_analysis(i_range = 100, tp, lb, ub)
results1$summary_stats


results3 <- short_mle_analysis(i_range = 300, tp, lb, ub)
results3$summary_stats


####################### Censored data #############################

short_mle_analysis <- function(i_range, true_params, lower_bound, upper_bound, batch_size = 250) {
  k_range <- 1000
  num_batches <- ceiling(k_range / batch_size)
  results <- matrix(0, ncol = length(true_params), nrow = k_range, dimnames = list(NULL, names(true_params)))
  
  
  gen_data <- function() {
    ds <- c()  # Uncensored samples
    dc <- c()  # Censored samples
    q_u = 0.9
    for (i in 1:i_range) {
      sol <- function(x) {
        zp <- rnorm(1)
        meu <- exp(true_params['phi1'] + true_params['phi2'] * zp)
        alph <- ((q_u / (1 - q_u))^(1 / true_params['eta']) - 
                   true_params['beta'] * meu^true_params['gamma']) / meu^true_params['theta']
        alph * x^true_params['theta'] + 
          true_params['beta'] * x^true_params['gamma'] - 
          ((u) / (1 - u))^(1 / true_params['eta'])
      }
      u <- runif(1)
      sample <- uniroot(sol, c(0.001, 100), extendInt = "yes")$root
      
      u <- 0.7                                                           #### change censoring cp = 1-u
      t_yc <- uniroot(sol, c(0.001, 100), extendInt = "yes")$root
      censor_time <- runif(1, 0, t_yc)
      
      if (sample < censor_time) {
        ds <- c(ds, sample)
      } else {
        dc <- c(dc, censor_time)
      }
    }
    
    return(list(ds = ds, dc = dc))
  }
  
  log_likelihood <- function(x, data) {
    the <- x[1]; bet <- x[2]; gam <- x[3]; et <- x[4]; ph1 <- x[5]; ph2 <- x[6]
    q_u = 0.9
    ds <- data[[1]]  # Uncensored samples
    dc <- data[[2]]  # Censored samples
    sum_ds <- 0
    sum1_ds <- 0
    sum2_ds <- 0
    for (j in 1:length(ds)) {
      z1 <- rnorm(1)
      meu <- exp(ph1 + ph2 * z1)
      alp <- ((q_u / (1 - q_u))^(1 / et) - bet * meu^gam) / meu^the
      if (alp > 0) {
        sum_ds <- sum_ds + log(alp * the * ds[j]^(the - 1) + bet * gam * ds[j]^(gam - 1))
        sum1_ds <- sum1_ds + log(alp * ds[j]^the + bet * ds[j]^gam)
        sum2_ds <- sum2_ds + log(1 + (alp * ds[j]^the + bet * ds[j]^gam)^et)
      } else {
        return(-1e10)  # Penalty
      }
    }
    sum3_dc <- 0
    for (k in 1:length(dc)) {
      z1 <- rnorm(1)
      meu <- exp(ph1 + ph2 * z1)
      alp <- ((q_u / (1 - q_u))^(1 / et) - bet * meu^gam) / meu^the
      if (alp > 0) {
        sum3_dc <- sum3_dc + log(1 + (alp * dc[k]^the + bet * dc[k]^gam)^et)
      } else {
        return(-1e10)  # Penalty
      }
    }
    ff <- length(ds) * log(et) + sum_ds + (et - 1) * sum1_ds - 2 * sum2_ds
    ffc <- -sum3_dc
    return(ff + ffc)
  }
  
  # Set up parallel backend
  num_cores <- detectCores() - 1  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 24) 
  
  # Run the analysis in batches
  for (batch in 1:num_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, k_range)
    batch_results <- foreach(k = batch_start:batch_end, .combine = rbind, .packages = "stats") %dopar% {
      successful <- FALSE
      result <- NULL
      while (!successful) {
        data <- gen_data()
        optim_result <- optim(true_params, log_likelihood, data = data,
                              lower = lower_bound, upper = upper_bound,
                              method = "L-BFGS-B", control = list(fnscale = -1), hessian = TRUE)
        if (optim_result$value > -1e10) {  
          successful <- TRUE
          result <- optim_result$par
        }
      }
      return(result)
    }
    
    results[batch_start:batch_end, ] <- batch_results
    print(paste("Batch", batch, "completed"))
  }
  
  stopCluster(cl)
  
  est_list <- list()
  summary_stats <- sapply(seq_along(true_params), function(i) {
    est <- results[, i]
    est_list[[names(true_params)[i]]] <<- est  
    rmse <- sqrt(mean((est - true_params[i])^2))
    sd <- sd(est)
    rel_bias <- mean(abs(est - true_params[i]) / true_params[i])
    c(mean = mean(est), rmse = rmse, sd = sd, rel_bias = rel_bias)
  })
  
  list(summary_stats = as.data.frame(t(summary_stats), row.names = names(true_params)),
       estimates = est_list)  
}



tp <- c(theta = 1.8, beta = 1.3, gamma = 0.5, eta = 1.1, phi1 = 0.8, phi2 = 0.7)
lb = tp - 0.49
ub = tp + 0.51



results6 <- short_mle_analysis(i_range = 20, tp, lb, ub)
results6$summary_stats

results7 <- short_mle_analysis(i_range = 30, tp, lb, ub)
results7$summary_stats


results0 <- short_mle_analysis(i_range = 50, tp, lb, ub)
results0$summary_stats

results1 <- short_mle_analysis(i_range = 100, tp, lb, ub)
results1$summary_stats


results3 <- short_mle_analysis(i_range = 300, tp, lb, ub)
results3$summary_stats
