rm(list = ls())
data = c(4690, 740, 1010, 1190, 2450, 1390, 350,6095, 3000, 1458,  
         550, 1690, 745, 1225, 1480, 245, 600, 246, 1805,
         258, 114, 312, 772, 498, 162, 444, 1464, 132, 1266, 300,
         520, 1240,  222, 144, 745, 396,
         510, 252, 408, 528, 690, 
         714, 348, 546, 174, 696, 
         294, 234, 288,444, 390, 168, 558, 288 )/100

dc = c(6200, 1740, 2440, 2600,1000,900)/100

zd = append(append(rep(52.5,19), rep(55,17)), rep(57.5,18))

zdc = c(52.5, 55, 55, 55, 57.5, 57.5)

####################################################
#         LAWD Regression 
####################################################

LogL <- function(x) {
  the <- x[1]
  bet <- x[2]
  gam <- x[3]
  et <- x[4]
  ph1 <- x[5]
  ph2 <- x[6]
  sum <- 0
  sum1 <- 0
  sum2 <- 0
  mu <- numeric(length(data))
  
  for(j in 1:length(data)) {
    z1 <- zd[j]
    mu[j] <- exp(ph1 + ph2 * z1)
    
    alp <- (1 - bet * (mu[j]^gam)) / (mu[j]^the)
    
    sum <- sum + log(alp * the * (data[j]^(the - 1)) + bet * gam * (data[j]^(gam - 1)))
    sum1 <- sum1 + log(alp * (data[j]^the) + bet * (data[j]^gam))
    sum2 <- sum2 + log(1 + (alp * (data[j]^the) + bet * (data[j]^gam))^et)
  }
  
  sum3 <- 0
  for(k in 1:length(dc)) {
    z1 <- zdc[k]
    mu[k] <- exp(ph1 + ph2 * z1)
    
    alp <- (1 - bet * (mu[k]^gam)) / (mu[k]^the)
    
    sum3 <- sum3 + log(1 + (alp * (dc[k]^the) + bet * (dc[k]^gam))^et)
  }
  
  ff <- length(data) * log(et) + sum + (et - 1) * sum1 - 2 * sum2
  ffc <- -sum3
  
  return(-(ff + ffc))
}
res = optim(par = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
            fn = LogL,
            lower = c(0.0001, 0.001, 0.0001, 0.0001, -Inf, -Inf),
            upper = rep(Inf, 6),
            method = "L-BFGS-B", hessian = T
)


#1.0481676  0.0010000  0.2358246  1.8037058 13.2850121 -0.2073515 # again estimate beta with these

res = optim(par = c(1.0481676,  0.0010000*2,  0.2358246,  1.8037058, 13.2850121, -0.2073515),
            fn = LogL,
            lower = c(0.001, 0.01, 0.001, 0.001, -Inf, -Inf)/10,
            upper = rep(Inf, 6),
            method = "L-BFGS-B", hessian = T
)



var_cov_matrix <- solve(res$hessian)  # Inverse of Hessian
std_errors <- sqrt(diag(var_cov_matrix))  # Standard errors of parameters


z_scores <- res$par / std_errors
p_values <- 2 * (1 - pnorm(abs(z_scores)))

results1 <- data.frame(
  Parameter = c("theta", "beta", "gamma", "eta", "phi1", "phi2"),
  Estimate = res$par,
  Std_Error = std_errors,
  z_value = z_scores,
  p_value = p_values
)

print(results1)

the <- res$par[1]
bet <- res$par[2]
gam <- res$par[3]
et <- res$par[4]
ph1 <- res$par[5]
ph2 <- res$par[6]

loglik_vals1 <- numeric(length(data))  # Initialize result vector
mu1 <- numeric(length(data))           # Initialize mu vector

for (k in 1:length(data)) {
  z1 <- zd[k]
  mu1[k] <- exp(ph1 + ph2 * z1)
  
  alp <- (1 - bet * (mu1[k]^gam)) / (mu1[k]^the)
  
  loglik_vals1[k] <- 1 -log(1 + (alp * (data[k]^the) + bet * (data[k]^gam))^et)
}

loglik_vals <- numeric(length(dc))  # Initialize result vector
mu <- numeric(length(dc))           # Initialize mu vector

for (k in 1:length(dc)) {
  z1 <- zdc[k]
  mu[k] <- exp(ph1 + ph2 * z1)
  
  alp <- (1 - bet * (mu[k]^gam)) / (mu[k]^the)
  
  loglik_vals[k] <- -log(1 + (alp * (dc[k]^the) + bet * (dc[k]^gam))^et)
}

delta = append(rep(1,54), rep(0,6))
rMi <- numeric(length(delta))  # intermediate residual (r_Mi)
rDi <- numeric(length(delta))  # final deviance residual (r_Di)


rMi = append(loglik_vals1, loglik_vals)

for (k in 1:length(delta)) {
  
  if ((delta[k] - rMi[k]) > 0) {
    rDi[k] <- sign(rMi[k]) * sqrt(-2 * (rMi[k] + delta[k] * log(delta[k] - rMi[k])))
  } else {
    rDi[k] <- NA  
  }
}


mu2 = append(mu1,mu)
dac = append(data,dc)

u = numeric(60)
for (i in 1:60) {
  sol <- function(u) {
    ((1-bet*mu2[i]^gam)/(mu2[i]^the))*(dac[i]^the) + bet*(dac[i]^gam) - (u/(1-u))^(1/et)
  }
  u[i] <- uniroot(sol, c(0.001, 1), extendInt = "yes")$root
}


# Deviance-residual Plot
plot(u, sort(rDi), xlab = "Index", ylab = "deviance residuals", xaxt = "n",cex.lab = 1.5, cex=1.5,
     cex.axis = 1.3 )
x_ticks <- pretty(u)
axis(1, at = x_ticks, labels = c(0, 10, 20, 30, 40, 60),cex.axis = 1.3)  
abline(h = -2, lty = 3,lwd=2, col="red")
abline(h = 2, lty = 3,lwd=2,col="red")
abline(h = 0, lty = 2,lwd=2)

####################################################
#         Lexp Regression 
####################################################


LogL5 <- function(x) {
  et <- x[1]
  ph1 <- x[2]
  ph2 <- x[3]
  sum <- 0
  sum1 <- 0
  sum2 <- 0
  mu <- numeric(length(data))
  
  for(j in 1:length(data)) {
    z1 <- zd[j]
    mu[j] <- exp(ph1 + ph2 * z1)
    
    alp <- log(2)/mu[j]
    
    sum <- sum + log(alp) + alp*data[j]
    sum1 <- sum1 + log(exp(alp*data[j])-1)
    sum2 <- sum2 + log(1 + (exp(alp*data[j])-1)^et)
  }
  
  sum3 <- 0
  for(k in 1:length(dc)) {
    z1 <- zdc[k]
    mu[k] <- exp(ph1 + ph2 * z1)
    
    alp <-  log(2)/mu[k]
    
    sum3 <- sum3 + log(1 + (exp(alp*dc[k])-1)^et)
  }
  
  ff <- length(data) * log(et) + sum + (et - 1) * sum1 - 2 * sum2
  ffc <- -sum3
  
  return(-(ff + ffc))
}


res1 = optim(
  par = c(0.01, 0.01, 0.01)/100,
  fn = LogL5,
  lower = c(0.001, -Inf, -Inf)/100,
  upper = rep(Inf, 3),
  method = "L-BFGS-B", hessian = T
)


var_cov_matrix <- solve(res1$hessian)  # Inverse of Hessian
std_errors <- sqrt(diag(var_cov_matrix))  # Standard errors of parameters

z_scores <- res1$par / std_errors
p_values <- 2 * (1 - pnorm(abs(z_scores)))

# Combine results into a table
results2 <- data.frame(
  Parameter = c( "eta", "phi1", "phi2"),
  Estimate = res1$par,
  Std_Error = std_errors,
  z_value = z_scores,
  p_value = p_values
)

print(results2)

et <- res1$par[1]
ph1 <- res1$par[2]
ph2 <- res1$par[3]



loglik_vals1 <- numeric(length(data))  # Initialize result vector
mu1 <- numeric(length(data))           # Initialize mu vector

for (k in 1:length(data)) {
  z1 <- zd[k]
  mu1[k] <- exp(ph1 + ph2 * z1)
  
  alp <-  log(2)/mu1[k]
  
  loglik_vals1[k] <- 1 -log(1 + (exp(alp*data[k])-1)^et)
}

loglik_vals <- numeric(length(dc))  # Initialize result vector
mu <- numeric(length(dc))           # Initialize mu vector
for (k in 1:length(dc)) {
  z1 <- zdc[k]
  mu[k] <- exp(ph1 + ph2 * z1)
  
  alp <-  log(2)/mu[k]
  
  loglik_vals[k] <-  -log(1 + (exp(alp*dc[k])-1)^et)
}

rMi <- numeric(length(delta))  # intermediate residual (r_Mi)
rDi <- numeric(length(delta))  # final deviance residual (r_Di)

delta = append(rep(1,54), rep(0,6))

rMi = append(loglik_vals1, loglik_vals)

for (k in 1:length(delta)) {
  
  if ((delta[k] - rMi[k]) > 0) {
    rDi[k] <- sign(rMi[k]) * sqrt(-2 * (rMi[k] + delta[k] * log(delta[k] - rMi[k])))
  } else {
    rDi[k] <- NA  
  }
}

mu2 = append(mu1,mu)

u = numeric(60)

for (i in 1:60) {
  sol <- function(u) {
    (mu2[i]/log(2))*log(1+(u/(1-u))^(1/et)) - dac[i]
  }
  u[i] <- uniroot(sol, c(0.001, 1), extendInt = "yes")$root
}


# Deviance-residual Plot
plot(u, sort(rDi), xlab = "Index", ylab = "deviance residuals", xaxt = "n",cex.lab = 1.5, cex=1.5,
     cex.axis = 1.3 )
x_ticks <- pretty(u)
axis(1, at = x_ticks, labels = c(0, 10, 20, 30, 40, 60),cex.axis = 1.3)  
abline(h = -2, lty = 3,lwd=2, col="red")
abline(h = 2, lty = 3,lwd=2,col="red")
abline(h = 0, lty = 2,lwd=2)

