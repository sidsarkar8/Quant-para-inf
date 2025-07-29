## Importing all the required libraries 
library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)


#----------------------------------------------------
## Including source files
#source("Documents/Academic_Stuffs/Projects/AK/R codes/DW_DKW_CI_F.R")

## GLD quantile function for CSW parametrization
S1 <- function(u, lambda3,lambda4) {
  if (lambda3 == 0 & lambda4 == 0) {
    return(log(u) - log(1 - u))
  } else if (lambda3 == 0 & lambda4 != 0) {
    return(log(u) - (1 / lambda4) * ((1 - u)^lambda4 - 1))
  } else if (lambda3 != 0 && lambda4 == 0) {
    return((1 / lambda3) * (u^lambda3 - 1) - log(1 - u))
  } else {
    return((1 / lambda3) * (u^lambda3 - 1) - (1 / lambda4) * ((1 - u)^lambda4 - 1))
  }
}

S <- function(u, chi, xi){
  if(chi == 0 && xi == 0.5) return(S1(u,0,0))
  else if(chi !=0 && xi == 0.5*(1+chi)){
    alpha = 1/2*(1/2-xi)/sqrt(xi*(1-xi))
    lambda4 = 2*alpha
    return(S1(u,0,lambda4))
    } 
  else if(chi !=0 && xi == 0.5*(1-chi)){
    alpha = 1/2*(1/2-xi)/sqrt(xi*(1-xi))
    lambda3 = 2*alpha
    return(S1(u,lambda3,0))
  }
  else{
    alpha = 1/2*(1/2-xi)/sqrt(xi*(1-xi))
    beta = 1/2*chi/(sqrt(1-chi^2))
    lambda3 = alpha + beta
    lambda4 = alpha - beta
    return(S1(u,lambda3,lambda4))
  }
}


Q_csw <- function(u, mu_tilde, sigma_tilde, chi, xi) {
  numerator   <- S(u, chi, xi) - S(0.5, chi, xi)
  denominator <- S(0.75, chi, xi) - S(0.25, chi, xi)
  return(mu_tilde + sigma_tilde * (numerator / denominator))
}


#----------------------------------------------------
## Simulating from GLD (CSW parametrization)
rgldcsw = function(n,mu_tilde,sigma_tilde,chi,xi){
  if(sigma_tilde <= 0 || chi <= -1 || chi >= 1 || xi <= 0 || xi >= 1) return("Not in parameter space!")
  else{
    u = runif(n)
    samp = Q_csw(u,mu_tilde,sigma_tilde,chi,xi)
    return(samp)
  }
}

#----------------------------------------------------
## CI for quantile and IQR(Inter quantile range)
CI_quantile = function(samp,L,U,alpha,u){
  if(u != 1 || u!=0){
    sorted_samp = sort(samp)
    u1 = which(U >= u)[1]
    l1 = tail(which((1-L) >= 1-u),1)
    return(c(sorted_samp[u1],sorted_samp[l1]))
  }
  else if(u == 1){
    u1 = which(U >= u)[1]
    return(c(sorted_samp[u1],Inf))
  }
  else if(u == 0){
    l1 = tail(which((1-L) >= 1-u),1)
    return(c(-Inf,sorted_samp[l1]))
  }
}

CI_iqr = function(samp,L,U,alpha,q){
  if(q > 0.5) return("IQR is negative. Choose q <= 0.5")
  else if(q < 0.5) {
    CI_u1 = CI_quantile(samp,L,U,alpha,q)
    CI_u2 = CI_quantile(samp,L,U,alpha,(1-q))
    l1 = CI_u1[1]
    u1 = CI_u1[2]
    l2 = CI_u2[1]
    u2 = CI_u2[2]
    l = max((l2-u1),0)
    u = u2 - l1
    return(c(l,u))
  }
  else return(c(0,0))
}

#------------------------------------------------------
## Width and Coverage of CI for mu and sigma of CSW parametrization
n_values = seq(30,100,by=10)
mu_tilde_seq = c(-1,0,1)
sigma_tilde_seq = c(1/2,1,2)

## CI for mu
results_mu <- matrix(NA, nrow = length(mu_tilde_seq) * length(n_values), ncol = 4)
colnames(results_mu) <- c("mu", "n", "coverage", "width")

row_index <- 1

for (n in n_values) {
  LU <- get_bds_dw(0.05, n)  # Compute once per n
  L <- LU$low
  U <- LU$up
  
  for (mu in mu_tilde_seq) {
    coverage_vec <- numeric(100)
    width_vec <- numeric(100)
    
    for (i in 1:100) {
      samp <- rgldcsw(n, mu, 1, 0, 1/2)
      ci <- CI_quantile(samp, L, U, 0.05, 1/2)
      width_vec[i] <- ci[2] - ci[1]
      coverage_vec[i] <- ifelse((mu - ci[1]) * (mu - ci[2]) < 0, 1, 0)
    }
    
    results_mu[row_index, ] <- c(mu, n, mean(coverage_vec), mean(width_vec))
    row_index <- row_index + 1
  }
}

results_mu_df <- as.data.frame(results_mu)
ggplot(results_mu_df, aes(x = n, y = coverage, color = factor(mu))) +
  geom_line(size = 1) +
  labs(
    title = "Average Coverage vs Sample Size",
    x = "Sample Size (n)",
    y = "Coverage",
    color = "Mu"
  ) +
  ylim(0.5, 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )


ggplot(results_mu_df, aes(x = n, y = width, color = factor(mu))) +
  geom_line(size = 1) +
  labs(
    title = "Average Width vs Sample Size",
    x = "Sample Size (n)",
    y = "Interval Width",
    color = "Mu"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )



## CI for sigma
results_sigma <- matrix(NA, nrow = length(sigma_tilde_seq) * length(n_values), ncol = 4)
colnames(results_sigma) <- c("sigma", "n", "coverage", "width")

row_index <- 1

for (n in n_values) {
  LU <- get_bds_dw(0.05, n)  # Compute once per n
  L <- LU$low
  U <- LU$up
  
  for (sigma in sigma_tilde_seq) {
    coverage_vec <- numeric(100)
    width_vec <- numeric(100)
    
    for (i in 1:100) {
      samp <- rgldcsw(n, 0, sigma, 0, 1/2)
      ci <- CI_iqr(samp,L,U,0.05,1/4)
      width_vec[i] <- ci[2] - ci[1]
      coverage_vec[i] <- ifelse((sigma - ci[1]) * (sigma - ci[2]) < 0, 1, 0)
    }
    
    results_sigma[row_index, ] <- c(sigma, n, mean(coverage_vec), mean(width_vec))
    row_index <- row_index + 1
  }
}

results_sigma_df <- as.data.frame(results_sigma)
ggplot(results_sigma_df, aes(x = n, y = coverage, color = factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = "Average Coverage vs Sample Size",
    x = "Sample Size (n)",
    y = "Coverage",
    color = "Sigma"
  ) +
  ylim(0.5, 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )

ggplot(results_sigma_df, aes(x = n, y = width, color = factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = "Average Width vs Sample Size",
    x = "Sample Size (n)",
    y = "Interval Width",
    color = "Sigma"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )



#write.csv(results_mu_df,"Documents/Academic_Stuffs/Projects/AK/gld_mu_ci.csv")
#write.csv(results_sigma_df,"Documents/Academic_Stuffs/Projects/AK/gld_sigma_ci.csv")


#results_mu_df = read.csv("Documents/Academic_Stuffs/Projects/AK/gld_mu_ci.csv")
#results_sigma_df = read.csv("Documents/Academic_Stuffs/Projects/AK/gld_sigma_ci.csv")






