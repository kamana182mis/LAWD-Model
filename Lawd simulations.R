# Load necessary libraries
rm(list=ls())

library(doParallel)
library(foreach)

short_mle_analysis <- function(i_range, true_params, lower_bound, upper_bound, batch_size = 100) {
  k_range <- 1000
  num_batches <- ceiling(k_range / batch_size)
  results <- matrix(0, ncol = length(true_params), nrow = k_range, dimnames = list(NULL, names(true_params)))
  gen_data <- function() {
    sapply(1:i_range, function(i) {
      u <- runif(1)
      sol <- function(x) true_params['alpha'] * x^true_params['theta'] + 
        true_params['beta'] * x^true_params['gamma'] - 
        ((u) / (1 - u))^(1 / true_params['eta'])
      uniroot(sol, c(0.001, 100), extendInt = "yes")$root
    })
  }
  
  log_likelihood <- function(x, data) {
    alp <- x[1]; the <- x[2]; bet <- x[3]; gam <- x[4]; et <- x[5]
    sum_term1 <- sum(log(alp * the * data^(the - 1) + bet * gam * data^(gam - 1)))
    sum_term2 <- sum(log(alp * data^the + bet * data^gam))
    sum_term3 <- sum(log(1 + (alp * data^the + bet * data^gam)^et))
    i_range * log(et) + sum_term1 + (et - 1) * sum_term2 - 2 * sum_term3
  }
  
  # Set up parallel backend
  num_cores <- detectCores() - 1  # Use all but one core
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Run the analysis in batches
  for (batch in 1:num_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, k_range)
    batch_results <- foreach(k = batch_start:batch_end, .combine = rbind, .packages = "stats") %dopar% {
      data <- gen_data()
      optim(true_params, log_likelihood, data = data,
            lower = lower_bound, upper = upper_bound,
            method = "L-BFGS-B", control = list(fnscale = -1), hessian = TRUE)$par
    }
    
    # Store results for the current batch
    results[batch_start:batch_end, ] <- batch_results
    
    # Print progress message
    print(paste("Batch", batch, "completed"))
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Store `est` values and summary statistics
  est_list <- list()
  summary_stats <- sapply(seq_along(true_params), function(i) {
    est <- results[, i]
    est_list[[names(true_params)[i]]] <<- est  # Save `est` values in `est_list`
    rmse <- sqrt(mean((est - true_params[i])^2))
    sd <- sd(est)
    rel_bias <- mean(abs(est - true_params[i]) / true_params[i])
    c(mean = mean(est), rmse = rmse, sd = sd, rel_bias = rel_bias)
  })
  
  list(summary_stats = as.data.frame(t(summary_stats), row.names = names(true_params)),
       estimates = est_list)  # Return both summary stats and est values
}

# Example of running the function with lower and upper bounds

set.seed(22)
tp = c(alpha = 2.4, theta = 2.8, beta = 2.9, gamma = 2.4, eta = 0.8)
# tp = c(alpha = 1.2, theta = 1.6, beta = 1.4, gamma = 1.3, eta = 1.5)
# tp = c(alpha = 2.5, theta = 0.8, beta = 2.6, gamma = 0.7, eta = 1.6)
# tp = c(alpha = 0.8, theta = 2.6, beta = 0.7, gamma = 2.5, eta = 2.6)

lb = tp - 0.5
ub = tp + 0.5



results5 <- short_mle_analysis(i_range = 15, tp, lb, ub)
results5$summary_stats

results6 <- short_mle_analysis(i_range = 30, tp, lb, ub)
results6$summary_stats

results0 <- short_mle_analysis(i_range = 50, tp, lb, ub)
results0$summary_stats


results1 <- short_mle_analysis(i_range = 100, tp, lb, ub)
results1$summary_stats

results2 <- short_mle_analysis(i_range = 200, tp, lb, ub)
results2$summary_stats

results3 <- short_mle_analysis(i_range = 300, tp, lb, ub)
results3$summary_stats

results4 <- short_mle_analysis(i_range = 500, tp, lb, ub)
results4$summary_stats


###################TTT PLOTS###################################

hrf_function = function(al, th, be, ga, et){
  srf = 1- 1/(1 + (al * x^th + be * x^ga)^(-et))
  pdf = (et * (al * th * x^(th-1) + be * ga * x^(ga-1)) * ((al * x^th + be * x^ga)^(et-1))) /
    ((1 + (al * x^th + be * x^ga)^(et))^2)
  hrf = pdf/srf
  return(hrf)
}

png(file.path(folder_path, "hrf_1.png"), width = 6, height = 8, units = "in", res = 300)
par(mar = c(5, 6, 4, 2) + 0.1)
x = seq(0.1,7,.01)



# Compute the y-range for all curves
y_max = max(
  hrf_function(2.4, 2.8, 2.9, 2.4, 0.8),
  hrf_function(1.2, 1.6, 1.4, 1.3, 1.5),
  hrf_function(2.5, 0.8, 2.6, 0.7, 1.6),
  hrf_function(0.8, 2.6, 0.7, 2.5, 2.6)
)

# Now plot with the adjusted y-limits
plot(x, hrf_function(2.4, 2.8, 2.9, 2.4, 0.8), type = "l", lwd = 2, lty = 1, col = "red",
     xlab = "X", ylab = "HF", ylim = c(0, y_max), cex.lab = 2, cex.axis = 2)
lines(x, hrf_function(1.2, 1.6, 1.4, 1.3, 1.5), lwd = 2, lty = 2, col = "blue")
lines(x, hrf_function(2.5, 0.8, 2.6, 0.7, 1.6), lwd = 2, lty = 3, col = "purple")
lines(x, hrf_function(0.8, 2.6, 0.7, 2.5, 2.6), lwd = 2, lty = 4, col = "green")


legend(4, 6, legend = c(
  expression(paste(theta == 4)),
  expression(paste(theta == 2)),
  expression(paste(theta == 1.5)),
  expression(paste(theta == 0.5)),
  expression(paste(theta == 0.8))
), lwd = 2, lty = 1:5, col = c("red", "blue", "purple", "green", "orange"), cex = 1.5, bty = "n")

set.seed(123)

gen_data_once <- function(n, true_params){
  sapply(1:n, function(i) {
    u <- runif(1)
    sol <- function(x) true_params['alpha'] * x^true_params['theta'] + 
      true_params['beta'] * x^true_params['gamma'] - 
      ((u)/(1-u))^(1/true_params['eta'])
    uniroot(sol, c(0.001, 100), extendInt = "yes")$root
  })
}




set.seed(22)
tp = c(alpha = 2.4, theta = 2.8, beta = 2.9, gamma = 2.4, eta = 0.8)
lb = tp - 0.5
ub = tp + 0.5

x_sim <- gen_data_once(500, tp)

TTT <- function(x){
  x <- sort(x)
  n <- length(x)
  x <- c(0, x)
  num <- cumsum((n:1) * diff(x))
  den <- sum((n:1) * diff(x))
  c(0, num/den)
}

ttt <- TTT(x_sim)
u <- seq(0, 1, length.out = length(ttt))

plot(u, ttt, type = "l", lwd = 2,
     xlab = "u", ylab = "TTT(u)",
     main = "TTT plot for LAWD (bathtub case)")
abline(0, 1, lty = 2)



set.seed(22)
tp = c(alpha = 1.2, theta = 1.6, beta = 1.4, gamma = 1.3, eta = 1.5)
lb = tp - 0.5
ub = tp + 0.5

x_sim <- gen_data_once(500, tp)

TTT <- function(x){
  x <- sort(x)
  n <- length(x)
  x <- c(0, x)
  num <- cumsum((n:1) * diff(x))
  den <- sum((n:1) * diff(x))
  c(0, num/den)
}

ttt <- TTT(x_sim)
u <- seq(0, 1, length.out = length(ttt))

plot(u, ttt, type = "l", lwd = 2,
     xlab = "u", ylab = "TTT(u)",
     main = "TTT plot for LAWD (bathtub case)")
abline(0, 1, lty = 2)



set.seed(22)
tp = c(alpha = 2.5, theta = 0.8, beta = 2.6, gamma = 0.7, eta = 1.6)
lb = tp - 0.5
ub = tp + 0.5

x_sim <- gen_data_once(500, tp)

TTT <- function(x){
  x <- sort(x)
  n <- length(x)
  x <- c(0, x)
  num <- cumsum((n:1) * diff(x))
  den <- sum((n:1) * diff(x))
  c(0, num/den)
}

ttt <- TTT(x_sim)
u <- seq(0, 1, length.out = length(ttt))

plot(u, ttt, type = "l", lwd = 2,
     xlab = "u", ylab = "TTT(u)",
     main = "TTT plot for LAWD (Decreasing)")
abline(0, 1, lty = 2)



set.seed(22)
tp = c(alpha = 0.8, theta = 2.6, beta = 0.7, gamma = 2.5, eta = 2.6)
lb = tp - 0.5
ub = tp + 0.5

x_sim <- gen_data_once(500, tp)

TTT <- function(x){
  x <- sort(x)
  n <- length(x)
  x <- c(0, x)
  num <- cumsum((n:1) * diff(x))
  den <- sum((n:1) * diff(x))
  c(0, num/den)
}

ttt <- TTT(x_sim)
u <- seq(0, 1, length.out = length(ttt))

plot(u, ttt, type = "l", lwd = 2,
     xlab = "u", ylab = "TTT(u)",
     main = "TTT plot for LAWD (Increasing)")
abline(0, 1, lty = 2)

