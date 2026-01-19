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
      
      u <- 0.7 #### change censoring cp = 1-u
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

##################### check alp > 0
# z_values <- numeric(1000)
# alp_values <- numeric(1000)
# 
# # Loop and calculate z and alp for each iteration
# for (i in 1:1000) {
#   z <- rnorm(1, 0, 1)
#   meu <- exp(phi1 + z * phi2)
#   alp <- ((u / (1 - u))^(1 / eta) - beta * (meu^gamma)) / (meu^eta)
#   
#   # Store z and alp values
#   z_values[i] <- z
#   alp_values[i] <- alp
# }
# alp1 = alp_values[alp_values>0]
# alp2 = (alp1[alp1<3])
# 
# # Plot z vs. alp
# #plot(z_values, alp_values, main="Plot of z vs. alp",
# xlab="z", ylab="alp", pch=19, col="blue")


##################### cp = 0
# theta = 1.8; beta = 1.3;  gamma = 0.6; eta = 1.1; phi1 = 0.5; phi2 = 0.7
# q_u = 0.8
# 
# samples <- vector("list", 1000)
# for (sample_idx in 1:1000) {
#   data <- numeric(50)
#   meu <- numeric(50)
#   for (i in 1:50) {
#     u <- runif(1, 0, 1)
#     sol <- function(x) {
#       zp <- rnorm(1, 0, 1)
#       meu[i] <- exp(phi1 + phi2 * zp)
#       alph <- ((q_u / (1 - q_u))^(1 / eta) - beta * (meu[i]^gamma)) / (meu[i]^theta)
#       alph * (x^theta) + beta * (x^gamma) - (u / (1 - u))^(1 / eta)
#     }
#     data[i] <- uniroot(sol, c(0.001, 100), extendInt = "yes")$root
#   }
#   samples[[sample_idx]] <- data
# }
# 
# 
# 
# results <- matrix(NA, nrow = 1000, ncol = 6)
# 
# for(i in 1:1000) {
#   data <- samples[[i]]
#   LogL <- function(x) {
#     the <- x[1]
#     bet <- x[2]
#     gam <- x[3]
#     et <- x[4]
#     ph1 <- x[5]
#     ph2 <- x[6]
# 
#     sum <- 0
#     sum1 <- 0
#     sum2 <- 0
#     mu <- numeric(length(data))
# 
#     for(j in 1:length(data)) {
#       z1 <- rnorm(1, 0, 1)
#       mu[j] <- exp(ph1 + ph2 * z1)
# 
#       # Calculate alp and check if it is positive
#       alp <- ((q_u / (1 - q_u))^(1 / et) - bet * (mu[j]^gam)) / (mu[j]^the)
# 
#       # Only proceed if alp is positive
#       if (alp > 0) {
#         sum <- sum + log(alp * the * (data[j]^(the - 1)) + bet * gam * (data[j]^(gam - 1)))
#         sum1 <- sum1 + log(alp * (data[j]^the) + bet * (data[j]^gam))
#         sum2 <- sum2 + log(1 + (alp * (data[j]^the) + bet * (data[j]^gam))^et)
#       } else {
#         # Penalize if alp is not positive
#         return(-1e10)  # Large negative penalty to discourage these parameter values
#       }
#     }
# 
#     # Compute final log-likelihood
#     ff <- length(data) * log(et) + sum + (et - 1) * sum1 - 2 * sum2
#     return(ff)
#   }
# 
#   res <- optim(
#     par = c(1.8, 1.3, 0.6, 1.1, 0.5, 0.7),
#     fn = LogL,
#     lower = rep(0.01,6),
#       #c(0.8, 0.3, 0.1, 0.01, 0.1, 0.2),  # Define lower bounds for all parameters
#     upper = rep(5,6),
#       #c(2.8, 2.3, 1.6, 2, 1, 1.2),     # Define upper bounds for all parameters
#     method = "L-BFGS-B",
#     control = list(fnscale = -1),
#     hessian = TRUE
#   )
#   results[i, ] <- res$par
# }
# 
# 
# col_means <- colMeans(results, na.rm = TRUE)
# print(col_means)


# results <- matrix(NA, nrow = 1000, ncol = 6)
# 
# 
# for(i in 1:1000) {
#   data <- y_new_list[[i]]
#   dc <- y_newc_list[[i]]
#   LogL <- function(x) {
#     the <- x[1]
#     bet <- x[2]
#     gam <- x[3]
#     et <- x[4]
#     ph1 <- x[5]
#     ph2 <- x[6]
#     
#     sum <- 0
#     sum1 <- 0
#     sum2 <- 0
#     mu <- numeric(length(data)) 
#     
#     for(j in 1:length(data)) {
#       z1 <- rnorm(1, 0, 1)
#       mu[j] <- exp(ph1 + ph2 * z1)
#       
#       # Calculate alp and check if it is positive
#       alp <- ((q_u / (1 - q_u))^(1 / et) - bet * (mu[j]^gam)) / (mu[j]^the)
#       
#       # Only proceed if alp is positive
#       if (alp > 0) {
#         sum <- sum + log(alp * the * (data[j]^(the - 1)) + bet * gam * (data[j]^(gam - 1)))
#         sum1 <- sum1 + log(alp * (data[j]^the) + bet * (data[j]^gam))
#         sum2 <- sum2 + log(1 + (alp * (data[j]^the) + bet * (data[j]^gam))^et)
#       } else {
#         # Penalize if alp is not positive
#         return(-1e10)  # Large negative penalty to discourage these parameter values
#       }
#     }
#     
#     sum3 <- 0 
#     mu <- numeric(length(dc)) 
#     for(k in 1:length(dc)){
#       z1 <- rnorm(1, 0, 1)
#       mu[k] <- exp(ph1 + ph2 * z1)
#       alp <- ((q_u / (1 - q_u))^(1 / et) - bet * (mu[k]^gam)) / (mu[k]^the)
#       if (alp > 0) {
#         sum3 <- sum3 + log(1 + (alp * (dc[k]^the) + bet * (dc[k]^gam))^et)
#       } else {
#         return(-1e10)  
#       }
#     }
#     
#     ff <- length(data) * log(et) + sum + (et - 1) * sum1 - 2 * sum2
#     ffc <- - sum3
#     
#     return(ff + ffc)
#   }
#   
#   res <- optim(
#     par = c(1.8, 1.3, 0.6, 1.0, 0.5, 0.7),
#     fn = LogL,
#     lower = rep(0.01,6),
#     #c(0.8, 0.3, 0.1, 0.01, 0.1, 0.2),  
#     upper = rep(5,6),
#     #c(2.8, 2.3, 1.6, 2, 1, 1.2),     
#     method = "L-BFGS-B",
#     control = list(fnscale = -1),
#     hessian = TRUE
#   )
#   results[i, ] <- res$par
# }
# 
# 
# col_means <- colMeans(results, na.rm = TRUE)
# print(col_means)






# 
# set.seed(24)
# 
# theta = 1.8; beta = 1.3;  gamma = 0.6; eta = 1.1; phi1 = 0.5; phi2 = 0.7
# q_u = 0.85
# 
# generate_simulation_samples <- function(n, num_simulations) {
#   samples <- vector("list", num_simulations)
#   
#   for (sample_idx in seq_len(num_simulations)) {
#     samples[[sample_idx]] <- replicate(n, {
#       u <- runif(1)
#       uniroot(function(x) {
#         zp <- rnorm(1)
#         meu <- exp(phi1 + phi2 * zp)
#         alph <- ((q_u / (1 - q_u))^(1 / eta) - beta * meu^gamma) / meu^theta
#         alph * x^theta + beta * x^gamma - (u / (1 - u))^(1 / eta)
#       }, c(0.001, 100), extendInt = "yes")$root
#     })
#   }
#   
#   return(samples)
# }
# 
# 
# generate_censoring_times <- function(n, num_simulations, u ) {
#   samples_yc <- vector("list", num_simulations)
#   for (i in 1:num_simulations) {
#     meu <- numeric(n)
#     t_values <- numeric(n)
#     yc <- numeric(n)
#     inverse_cdf <- function(t) {
#       zp <- rnorm(1, 0, 1)
#       meu_t <- exp(phi1 + phi2 * zp)
#       alph_t <- ((q_u / (1 - q_u))^(1 / eta) - beta * (meu_t^gamma)) / (meu_t^eta)
#       alph_t * (t^theta) + beta * (t^gamma) - (u / (1 - u))^(1 / eta)  
#     }
#     for (j in 1:n) {
#       t_values[j] <- uniroot(inverse_cdf, c(0.001, 100), extendInt = "yes")$root
#       yc[j] <- runif(1, 0, t_values[j])
#     }
#     samples_yc[[i]] <- yc
#   }
#   
#   return(samples_yc)
# }
# 
# n = 100
# num_simulations = 1000
# 
# samples <- generate_simulation_samples(n, num_simulations)
# samples_yc <- generate_censoring_times(n, num_simulations, u = 0.9)
# 
# 
# 
# 
# y_newc_list = vector("list", num_simulations)  
# y_new_list = vector("list", num_simulations) 
# for(i in 1:num_simulations){
#   y_newc = numeric() 
#   y_new = numeric()
#   for(j in 1:n){
#     if(samples[[i]][j] > samples_yc[[i]][j]){
#       y_newc = c(y_newc, samples_yc[[i]][j])
#     }
#     else{
#       y_new = c(y_new, samples[[i]][j])
#     }
#   }
#   y_newc_list[[i]] = y_newc 
#   y_new_list[[i]] = y_new
# }
# 
# 
# 
# penalty_count <- 0
# success_count <- 0
# results <- matrix(NA, nrow = 1000, ncol = 6)
# 
# i <- 1 
# while (success_count < 1000) {
#   data <- y_new_list[[i]]
#   dc <- y_newc_list[[i]]
#   
#   LogL <- function(x) {
#     the <- x[1]
#     bet <- x[2]
#     gam <- x[3]
#     et <- x[4]
#     ph1 <- x[5]
#     ph2 <- x[6]
#     
#     sum <- 0
#     sum1 <- 0
#     sum2 <- 0
#     mu <- numeric(length(data)) 
#     
#     for(j in 1:length(data)) {
#       z1 <- rnorm(1, 0, 1)
#       mu[j] <- exp(ph1 + ph2 * z1)
#       
#       alp <- ((q_u / (1 - q_u))^(1 / et) - bet * (mu[j]^gam)) / (mu[j]^the)
#       
#       if (alp > 0) {
#         sum <- sum + log(alp * the * (data[j]^(the - 1)) + bet * gam * (data[j]^(gam - 1)))
#         sum1 <- sum1 + log(alp * (data[j]^the) + bet * (data[j]^gam))
#         sum2 <- sum2 + log(1 + (alp * (data[j]^the) + bet * (data[j]^gam))^et)
#       } else {
#         return(-1e10)  
#       }
#     }
#     
#     sum3 <- 0 
#     for(k in 1:length(dc)) {
#       z1 <- rnorm(1, 0, 1)
#       mu[k] <- exp(ph1 + ph2 * z1)
#       
#       alp <- ((q_u / (1 - q_u))^(1 / et) - bet * (mu[k]^gam)) / (mu[k]^the)
#       if (alp > 0) {
#         sum3 <- sum3 + log(1 + (alp * (dc[k]^the) + bet * (dc[k]^gam))^et)
#       } else {
#         return(-1e10)  
#       }
#     }
#     
#     ff <- length(data) * log(et) + sum + (et - 1) * sum1 - 2 * sum2
#     ffc <- -sum3
#     
#     return(ff + ffc)
#   }
#   
#   
#   res <- optim(
#     par = c(1.8, 1.3, 0.6, 1.0, 0.5, 0.7),
#     fn = LogL,
#     lower = c(1.8, 1.3, 0.6, 1.0, 0.5, 0.7)-0.5,
#     upper = c(1.8, 1.3, 0.6, 1.0, 0.5, 0.7)+0.5,
#     method = "L-BFGS-B",
#     control = list(fnscale = -1),
#     hessian = TRUE
#   )
#   
#   if (res$value == -1e10) {
#     penalty_count <- penalty_count + 1
#   } else {
#     results[success_count + 1, ] <- res$par
#     success_count <- success_count + 1
#     i <- i + 1
#   }
# }
# 
# true_vals <- c(theta = 1.8, beta = 1.3,  gamma = 0.6, eta = 1.1, phi1 = 0.5, phi2 = 0.7)
# 
# estimates <- list(theta = results[,1], beta = results[,2], gamma = results[,3],
#                   eta = results[,4], phi1 = results[,5], phi2 = results[,6])
# 
# 
# compute_metrics <- function(est, true_val) {
#   c(mean = mean(est), RMSE = sqrt(mean((est - true_val)^2)), sd = sd(est),
#     rel_bias = mean(abs(est - true_val)/true_val))
# }
# 
# results_df <- do.call(rbind, lapply(seq_along(estimates), 
#                                     function(i) data.frame(Parameter = names(estimates)[i], 
#                                                            t(compute_metrics(estimates[[i]], true_vals[i])))))
# results_df

