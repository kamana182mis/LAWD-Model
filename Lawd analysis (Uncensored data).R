rm(list = ls())
library(moments)

# ======================= 1. LAWD =======================
fit_distribution1 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood1 <- function(x) {
    alp <- x[1]; the <- x[2]; bet <- x[3]; gam <- x[4]; et <- x[5]
    sum_term1 <- sum(log(alp * the * data^(the - 1) + bet * gam * data^(gam - 1)))
    sum_term2 <- sum(log(alp * data^the + bet * data^gam))
    sum_term3 <- sum(log(1 + (alp * data^the + bet * data^gam)^et))
    ll <- n * log(et) + sum_term1 + (et - 1) * sum_term2 - 2 * sum_term3
    return(-ll)
  }
  
  res <- optim(log_likelihood1, par = par, lower = lower, upper = upper,
               method = "L-BFGS-B", hessian = TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent = TRUE)
  if (inherits(se, "try-error")) se <- rep(NA, length(res$par))
  
  nll <- -res$value
  k <- length(res$par)
  
  ks_result <- ks.test(data, function(x)
    1 / (1 + (res$par[1]*x^(res$par[2]) + res$par[3]*x^(res$par[4]))^(-res$par[5])))
  
  data.frame(Model="LAWD", Conv=res$convergence, LL=nll,
             AIC=2*k-2*nll, BIC=k*log(n)-2*nll,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 2. LMWD =======================
fit_distribution2 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood2 <- function(x) {
    alp <- x[1]; the <- x[2]; bet <- x[3]; et <- x[4]
    ll <- n*log(et) +
      sum(log(alp*the*data^(the-1) + bet)) +
      (et-1)*sum(log(alp*data^the + bet*data)) -
      2*sum(log(1+(alp*data^the + bet*data)^et))
    -ll
  }
  
  res <- optim(log_likelihood2, par = par, lower = lower, upper = upper,
               method="L-BFGS-B", hessian=TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent=TRUE)
  if(inherits(se,"try-error")) se <- rep(NA,length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x^(res$par[2])+res$par[3]*x)^(-res$par[4])))
  
  data.frame(Model="LMWD", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 3. LLFRD =======================
fit_distribution3 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood3 <- function(x) {
    alp <- x[1]; bet <- x[2]; et <- x[3]
    ll <- n*log(et) +
      sum(log(2*alp*data + bet)) +
      (et-1)*sum(log(alp*data^2 + bet*data)) -
      2*sum(log(1+(alp*data^2 + bet*data)^et))
    -ll
  }
  
  res <- optim(log_likelihood3, par = par, lower = lower, upper = upper,
               method="L-BFGS-B", hessian=TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent=TRUE)
  if(inherits(se,"try-error")) se <- rep(NA,length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x^2+res$par[2]*x)^(-res$par[3])))
  
  data.frame(Model="LLFRD", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 4. LWD =======================
fit_distribution4 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood4 <- function(x) {
    alp <- x[1]; the <- x[2]; et <- x[3]
    ll <- n*log(et) +
      sum(log(alp*the*data^(the-1))) +
      (et-1)*sum(log(alp*data^the)) -
      2*sum(log(1+(alp*data^the)^et))
    -ll
  }
  
  res <- optim(log_likelihood4, par = par, lower = lower, upper = upper,
               method="L-BFGS-B", hessian=TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent=TRUE)
  if(inherits(se,"try-error")) se <- rep(NA,length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x^(res$par[2]))^(-res$par[3])))
  
  data.frame(Model="LWD", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 5.  =======================

fit_distribution5 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood5 <- function(x) {
    alp <- x[1]; the <- 2; bet <- 0; gam <- 1; et <- x[2]
    
    ll <- n*log(et) +
      sum(log(alp*the*data^(the-1))) +
      (et-1)*sum(log(alp*data^the)) -
      2*sum(log(1+(alp*data^the)^et))
    
    -ll
  }
  
  res <- optim(par = par,
               fn = log_likelihood5,
               lower = lower,
               upper = upper,
               method = "L-BFGS-B",
               hessian = TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent = TRUE)
  if(inherits(se,"try-error")) se <- rep(NA, length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x^2)^(-res$par[2])))
  
  data.frame(Model="LRD", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}


# ======================= 6.  =======================

fit_distribution6 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood6 <- function(x) {
    alp <- x[1]; the <- 1; bet <- 0; gam <- 1; et <- x[2]
    
    ll <- n*log(et) +
      sum(log(alp*the*data^(the-1))) +
      (et-1)*sum(log(alp*data^the)) -
      2*sum(log(1+(alp*data^the)^et))
    
    -ll
  }
  
  res <- optim(par = par,
               fn = log_likelihood6,
               lower = lower,
               upper = upper,
               method = "L-BFGS-B",
               hessian = TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent = TRUE)
  if(inherits(se,"try-error")) se <- rep(NA, length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x)^(-res$par[2])))
  
  data.frame(Model="LED", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 7.  =======================

fit_distribution7 <- function(data, par, lower, upper) {
  n <- length(data)
  
  log_likelihood7 <- function(x) {
    alp <- x[1]; the <- 1; lam <- x[2]; gam <- 2; et <- x[3]
    bet <- lam/2
    
    ll <- n*log(et) +
      sum(log(alp*the*data^(the-1) + bet*gam*data^(gam-1))) +
      (et-1)*sum(log(alp*data^the + bet*data^gam)) -
      2*sum(log(1+(alp*data^the + bet*data^gam)^et))
    
    -ll
  }
  
  res <- optim(par = par,
               fn = log_likelihood7,
               lower = lower,
               upper = upper,
               method = "L-BFGS-B",
               hessian = TRUE)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent = TRUE)
  if(inherits(se,"try-error")) se <- rep(NA, length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1/(1+(res$par[1]*x + res$par[2]*(x^2)/2)^(-res$par[3])))
  
  data.frame(Model="LLiED", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(n)+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}



# ======================= 5. Lexp =======================
Lexp_Model <- function(data) {
  ll_lEd <- function(p) {
    al <- p[1]; lam <- p[2]
    n <- length(data)
    ll <- n*log(al*lam) + (al-1)*sum(log(exp(lam*data)-1)) +
      sum(lam*data) - 2*sum(log(1+(exp(lam*data)-1)^al))
    -ll
  }
  
  res <- optim(fn = ll_lEd, par = c(0.5, 0.5), lower = c(0.0001, 0.0001), 
               upper = c(Inf, Inf), method = "L-BFGS-B", hessian = T)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent=TRUE)
  if(inherits(se,"try-error")) se <- rep(NA,length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1-(1+(exp(res$par[2]*x)-1)^res$par[1])^(-1))
  
  data.frame(Model="Lexp", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(length(data))+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}

# ======================= 6. NH =======================
NH_Model <- function(data) {
  ll_lEd <- function(p) {
    al <- p[1]; lam <- p[2]
    n <- length(data)
    ll <- n + n*log(al*lam) + (al-1)*sum(log(lam*data+1)) -
      sum((1+lam*data)^al)
    -ll
  }
  
  
  res <- optim(fn = ll_lEd, par = c(0.5, 0.5)/5, lower = c(0.0001, 0.0001), 
               upper = c(Inf, Inf), method = "L-BFGS-B", hessian = T)
  
  se <- try(sqrt(diag(solve(res$hessian))), silent=TRUE)
  if(inherits(se,"try-error")) se <- rep(NA,length(res$par))
  
  ks_result <- ks.test(data, function(x)
    1-exp(1-(1+res$par[2]*x)^res$par[1]))
  
  data.frame(Model="NH", Conv=res$convergence, LL=-res$value,
             AIC=2*length(res$par)+2*res$value,
             BIC=length(res$par)*log(length(data))+2*res$value,
             KS_p=ks_result$p.value, KS_stat=ks_result$statistic,
             Params=paste(round(res$par,4),collapse=", "),
             SE=paste(round(se,5),collapse=", "))
}












# Example 8
df1 <- c(2, 10, 13, 23, 23, 28, 30, 65, 80, 88, 106, 143, 147, 173, 181, 212, 245, 247, 261, 266, 275, 293, 300, 300, 300, 300, 300, 300, 300, 300)
df1 <- df1 / 100
hist(df1)

par1 <- c(0.1, 0.1, 0.1, 0.08, 0.1) / 10
lower1 <- rep(0.0001, 5)
upper1 <- rep(Inf, 5)
result1 <- fit_distribution1(df1, par1, lower1, upper1)
result1


par2 <- c(0.1, 0.1, 0.1, 0.1) / 100
lower2 <- rep(0.0001, 4)
upper2 <- rep(Inf, 4)
result2 <- fit_distribution2(df1, par2, lower2, upper2)

par3 <- c(0.1, 0.1, 0.1)
lower3 <- rep(0.0001, 3)
upper3 <- rep(Inf, 3)
result3 <- fit_distribution3(df1, par3, lower3, upper3)

result4 <- fit_distribution4(df1, par3, lower3, upper3)

par5 <- c(0.1, 0.1)
lower5 <- rep(0.0001, 2)
upper5 <- rep(Inf, 2)
result5 <- fit_distribution5(df1, par5, lower5, upper5)
result6 <- fit_distribution6(df1, par5, lower5, upper5)

result7 <- fit_distribution7(df1, par3, lower3, upper3)
result8 = Lexp_Model(df1)
result9 = NH_Model(df1)

final_results1 <- rbind(result1, result2, result3, result4, result8, result9)

# Example 7
df2 <- c(0.1, 7, 36, 67, 84, 0.2, 11, 40, 67, 84, 1, 12, 45, 67, 84, 1, 18, 46, 67, 85, 1, 18, 47, 72, 85, 1, 18, 50, 75, 85, 1, 18, 55, 79, 85, 2, 18, 60, 82, 85, 3, 21, 63, 82, 86, 6, 32, 63, 83, 86)
df2 <- df2 / 10
hist(df2)

par1 <- c(0.1, 0.1, 0.1, 0.08, 0.1) / 10
lower1 <- rep(0.0001, 5)
upper1 <- rep(Inf, 5)
result1 <- fit_distribution1(df2, par1, lower1, upper1)

par2 <- c(0.1, 0.1, 0.1, 0.1) / 10
lower2 <- rep(0.0001, 4) * 10
upper2 <- rep(Inf, 4)
result2 <- fit_distribution2(df2, par2, lower2, upper2)

par3 <- c(0.1, 0.1, 0.1)
lower3 <- rep(0.0001, 3)
upper3 <- rep(Inf, 3)
result3 <- fit_distribution3(df2, par3, lower3, upper3)

result4 <- fit_distribution4(df2, par3, lower3, upper3)

par5 <- c(0.1, 0.1)
lower5 <- rep(0.0001, 2)
upper5 <- rep(Inf, 2)
result5 <- fit_distribution5(df2, par5, lower5, upper5)
result6 <- fit_distribution6(df2, par5, lower5, upper5)

result7 <- fit_distribution7(df2, par3, lower3, upper3)

result8 = Lexp_Model(df2)
result9 = NH_Model(df2)

final_results2 <- rbind(result1, result2, result3, result4, result8, result9)





# Example 5
df3 <- c(115, 1277, 181, 1290, 255, 1357, 418, 1369, 441, 1408, 461, 1455, 516, 1478, 739, 1519, 743, 1578, 789, 1578, 807, 1599, 865, 1603, 924, 1605, 983, 1696, 1025, 1735, 1062, 1799, 1063, 1815, 1165, 1852, 1191, 1899, 1222, 1925, 1222, 1965, 1251)
df3 <- df3 / 1000
hist(df3)

par1 <- c(0.1, 0.1, 0.1, 0.08, 0.1)
lower1 <- rep(0.0001, 5)
upper1 <- rep(Inf, 5)
result1 <- fit_distribution1(df3, par1, lower1, upper1)

par2 <- c(0.1, 0.1, 0.1, 0.1)
lower2 <- rep(0.0001, 4)
upper2 <- rep(Inf, 4)
result2 <- fit_distribution2(df3, par2, lower2, upper2)

par3 <- c(0.1, 0.1, 0.1)
lower3 <- rep(0.0001, 3)
upper3 <- rep(Inf, 3)
result3 <- fit_distribution3(df3, par3, lower3, upper3)

result4 <- fit_distribution4(df3, par3, lower3, upper3)

par5 <- c(0.1, 0.1)
lower5 <- rep(0.0001, 2)
upper5 <- rep(Inf, 2)
result5 <- fit_distribution5(df3, par5, lower5, upper5)
result6 <- fit_distribution6(df3, par5, lower5, upper5)

result7 <- fit_distribution7(df3, par3, lower3, upper3)

result8 = Lexp_Model(df3)
result9 = NH_Model(df3)

final_results3 <- rbind(result1, result2, result3, result4, result8, result9)



final_results1
final_results2
final_results3






####################
library(ggplot2)

library(ggpubr)  # For arranging plots

# Data
fdf1 <- data.frame(Value = c(2, 10, 13, 23, 23, 28, 30, 65, 80, 88, 
                             106, 143, 147, 173, 181, 212, 245, 247, 
                             261, 266, 275, 293, 300, 300, 300, 300, 
                             300, 300, 300, 300))

fdf2 <- data.frame(Value = c(0.1, 7, 36, 67, 84, 0.2, 11, 40, 67, 84, 1, 12, 45, 67, 84, 
                             1, 18, 46, 67, 85, 1, 18, 47, 72, 85, 1, 18, 50, 75, 85, 
                             1, 18, 55, 79, 85, 2, 18, 60, 82, 85, 3, 21, 63, 82, 86, 
                             6, 32, 63, 83, 86))

fdf3 <- data.frame(Value = c(115, 1277, 181, 1290, 255, 1357, 418, 1369, 441, 1408, 461, 
                             1455, 516, 1478, 739, 1519, 743, 1578, 789, 1578, 807, 
                             1599, 865, 1603, 924, 1605, 983, 1696, 1025, 1735, 1062, 
                             1799, 1063, 1815, 1165, 1852, 1191, 1899, 1222, 1925, 
                             1222, 1965, 1251))

summary(fdf1)
kurtosis(fdf1); skewness(fdf1)


summary(fdf2)
kurtosis(fdf2); skewness(fdf2)

summary(fdf3)
kurtosis(fdf3); skewness(fdf3)

var(fdf1); var(fdf2); var(fdf3)


# Function to create boxplots
create_boxplot <- function(data, y_label) {
  ggplot(data, aes(y = Value, x = "")) + 
    geom_boxplot(
      color = "blue", fill = "blue", alpha = 0.5, notch = TRUE, 
      notchwidth = 0.6, outlier.colour = "red", outlier.fill = "red", outlier.size = 3, width = 0.5
    ) + 
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15),
      axis.title.y = element_text(vjust = 1.5, hjust = 0.4, size = 15),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    ggtitle(" ") +
    xlab("") +
    ylab(y_label)
}

# Create individual plots
plot1 <- create_boxplot(fdf1, "Failure & running times of devices")
plot2 <- create_boxplot(fdf2, "Lifetime of devices")
plot3 <- create_boxplot(fdf3, "Lifetimes of Blood cancer patients")

# Arrange the three plots in 1 row, 3 columns
combined_plot <- ggarrange(plot1, plot2, plot3, nrow = 1, ncol = 3)
combined_plot


plot_survival_distributions <- function(x, final_results, legend_pos = c(0.9, 0.7), x_label = "x") {
  x_sorted <- sort(x)
  empirical_survival <- 1 - ecdf(x_sorted)(x_sorted)
  survival_data <- data.frame(x = x_sorted, Empirical = empirical_survival)
  
  
  rs1 = as.numeric(unlist(strsplit(final_results$Params[1], ", ")))
  rs2 = as.numeric(unlist(strsplit(final_results$Params[2], ", ")))
  rs3 = as.numeric(unlist(strsplit(final_results$Params[3], ", ")))
  rs4 = as.numeric(unlist(strsplit(final_results$Params[4], ", ")))
  rs_exp = as.numeric(unlist(strsplit(final_results$Params[5], ", ")))
  rs_nh = as.numeric(unlist(strsplit(final_results$Params[6], ", ")))
  
  
  
  
  survival_LAWD <- function(x) { 1 - (1 / (1 + (rs1[1] * x^rs1[2] + rs1[3] * x^rs1[4])^(-rs1[5]))) }
  survival_LMWD <- function(x) { 1 - (1 / (1 + (rs2[1] * x^rs2[2] + rs2[3] * x)^(-rs2[4]))) }
  survival_LLFRD <- function(x) { 1 - (1 / (1 + (rs3[1] * x^2 + rs3[2] * x)^(-rs3[3]))) }
  survival_LWD <- function(x) { 1 - (1 / (1 + (rs4[1] * x^rs4[2])^(-rs4[3]))) }
  survival_Lexp <- function(x) { 1 - (1 - (1 + (exp(rs_exp[2] * x) - 1)^rs_exp[1] )^(-1)) }
  survival_NH <- function(x) { 1 - (1 - exp(1 - (1 + rs_nh[2] * x)^rs_nh[1])) }
  
  
  x_values <- seq(min(x), max(x), length.out = 1000)
  survival_values <- data.frame(
    x = rep(x_values, 6),
    survival = c(
      survival_LAWD(x_values),
      survival_LMWD(x_values),
      survival_LLFRD(x_values),
      survival_LWD(x_values),
      survival_Lexp(x_values),
      survival_NH(x_values)
    ),
    Distribution = rep(c("LAWD", "LMWD", "LLFRD", "LWD", "Lexp", "NHD"), each = 1000)
  )
  
  survival_values <- survival_values[complete.cases(survival_values), ]
  
  # Plot
  p <- ggplot() +
    geom_line(data = survival_data, aes(x = x, y = Empirical), color = "black", linetype = "solid", size = 0.8) +
    labs(x = x_label, y = "S(x)") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = legend_pos,
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank()
    ) +
    geom_line(data = survival_values, aes(x = x, y = survival, color = Distribution, linetype = Distribution), size = 0.7) +
    scale_color_manual(values = c("LAWD" = "red", "LMWD" = "blue", "LLFRD" = "green", "LWD" = "orange", 
                                  "Lexp" = "purple", "NHD" = "brown")) +
    scale_linetype_manual(values = c("LAWD" = "solid", "LMWD" = "dashed", "LLFRD" = "dotted", "LWD" = "twodash", 
                                     "Lexp" = "dotdash", "NHD" = "longdash"))
  
  print(p)
}


surv1 = plot_survival_distributions(df1, final_results1)
surv2 = plot_survival_distributions(df2, final_results2)
surv3 = plot_survival_distributions(df3, final_results3)


plot_fitted_distributions <- function(x, original_x, c, final_results, y_limit = NULL, x_limit = NULL, bin, legend_pos = c(0.5, 0.78), x_label = "x-axis label") {
  hist_data <- data.frame(x = x)
  p <- ggplot(hist_data, aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = bin, color = "black", fill = "white", alpha = 1) +
    labs(x = x_label, y = "Density") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = legend_pos,
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank()
    )
  
  if (!is.null(y_limit)) {
    p <- p + ylim(y_limit)
  }
  
  if (!is.null(x_limit)) {
    p <- p + xlim(x_limit)
  }
  
  # Extract parameters
  rs1 = as.numeric(unlist(strsplit(final_results$Params[1], ", ")))
  rs2 = as.numeric(unlist(strsplit(final_results$Params[2], ", ")))
  rs3 = as.numeric(unlist(strsplit(final_results$Params[3], ", ")))
  rs4 = as.numeric(unlist(strsplit(final_results$Params[4], ", ")))
  rs_exp = as.numeric(unlist(strsplit(final_results$Params[5], ", ")))
  rs_nh = as.numeric(unlist(strsplit(final_results$Params[6], ", ")))
  
  # Define density functions
  density_LAWD <- function(x) {
    (rs1[5] * (rs1[1] * rs1[2] * x^(rs1[2]-1) + rs1[3] * rs1[4] * x^(rs1[4]-1)) * 
       ((rs1[1] * x^rs1[2] + rs1[3] * x^rs1[4])^(rs1[5]-1))) /
      ((1 + (rs1[1] * x^rs1[2] + rs1[3] * x^rs1[4])^(rs1[5]))^2)
  }
  
  density_LMWD <- function(x) {
    (rs2[4] * (rs2[1] * rs2[2] * x^(rs2[2]-1) + rs2[3]) * 
       ((rs2[1] * x^rs2[2] + rs2[3] * x)^(rs2[4]-1))) /
      ((1 + (rs2[1] * x^rs2[2] + rs2[3] * x)^(rs2[4]))^2)
  }
  
  density_LLFRD <- function(x) {
    (rs3[3] * (rs3[1] * 2 * x + rs3[2]) * 
       ((rs3[1] * x^2 + rs3[2] * x)^(rs3[3]-1))) /
      ((1 + (rs3[1] * x^2 + rs3[2] * x)^(rs3[3]))^2)
  }
  
  density_LWD <- function(x) {
    (rs4[3] * (rs4[1] * rs4[2] * x^(rs4[2]-1)) * 
       ((rs4[1] * x^rs4[2])^(rs4[3]-1))) /
      ((1 + (rs4[1] * x^rs4[2])^(rs4[3]))^2)
  }
  
  density_Lexp <- function(x) {
    rs_exp[1] * rs_exp[2] * (exp(rs_exp[2] * x) * ((exp(rs_exp[2] * x) - 1)^(rs_exp[1]-1))) /
      ((1 + (exp(rs_exp[2] * x) - 1)^rs_exp[1])^2)
  }
  
  density_NH <- function(x) {
    rs_nh[2] * rs_nh[1] * ((1 + rs_nh[2] * x)^(rs_nh[1]-1)) * exp(1 - (1 + rs_nh[2] * x)^rs_nh[1])
  }
  
  # Generate x values
  x_values <- seq(min(x), max(x), length.out = 1000)
  
  # Compute densities
  density_values <- data.frame(
    x = rep(x_values, 6),
    density = c(
      density_LAWD(x_values),
      density_LMWD(x_values),
      density_LLFRD(x_values),
      density_LWD(x_values),
      density_Lexp(x_values),
      density_NH(x_values)
    ),
    Distribution = rep(c("LAWD", "LMWD", "LLFRD", "LWD", "Lexp", "NHD"), each = 1000)
  )
  x_original_mapping <- function(normalized_x) {
    normalized_x * c
  }
  
  # Plot density lines
  p <- p +
    geom_line(data = density_values, aes(x = x, y = density, color = Distribution, linetype = Distribution), size = 0.8) +
    scale_color_manual(values = c("LAWD" = "red", "LMWD" = "blue", "LLFRD" = "green", "LWD" = "orange", 
                                  "Lexp" = "purple", "NHD" = "brown")) +
    scale_linetype_manual(values = c("LAWD" = "solid", "LMWD" = "dashed", "LLFRD" = "dotted", "LWD" = "twodash", 
                                     "Lexp" = "dotdash", "NHD" = "longdash")) +
    scale_x_continuous(
      breaks = seq(0, max(x), length.out=10), # Adjust break positions as needed
      labels = scales::label_number(scale = c) # Show original x values
    )
  
  print(p)
}






# Example usage with x_limit:
p1 = plot_fitted_distributions(df1, fdf1, 100,  final_results1, bin = 20,  y_limit = c(0,1), legend_pos = c(0.8, 0.78), x_label = "")
p2 = plot_fitted_distributions(df2, fdf2, 10, final_results2,  bin = 20, y_limit = c(0,0.36), legend_pos = c(0.8, 0.78), x_label = "")
p3 = plot_fitted_distributions(df3, fdf3, 1000, final_results3, bin = 20 , legend_pos = c(0.9, 0.78),x_label = "")


