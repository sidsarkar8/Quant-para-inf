## Importing all the required libraries 
library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)

n = 30
LU = get_bds_dw(0.05,n)
L = LU$low
U = LU$up
lambda = -2
samp = rtlambda(n,lambda)
n = length(samp)
sorted_samp <- sort(abs(samp))
CI_indexwise = matrix(0,nrow = n, ncol = 3)
CI_indexwise[,1] = c(1:n)
N = 100

q_abs = function(p,lambda){
  if(abs(lambda) <= 1e-10) return(log((1+p)/(1-p)))
  else return(((1+p)^(lambda)-(1-p)^(lambda))/(lambda*2^(lambda)))
}
for(i in 1:n){
  f1 <- function(lambda) {q_abs(L[i], lambda) - sorted_samp[i]}
  alpha_i <- uniroot(f1, lower = -N, upper = N, tol = 1e-10)$root
  CI_indexwise[i,2] = alpha_i
  f2 <- function(lambda) {q_abs(U[i], lambda) - sorted_samp[i]}
  beta_i <- uniroot(f2, lower = -N, upper = N, tol = 1e-10)$root
  CI_indexwise[i,3] = beta_i
}

CI_indexwise = as.data.frame(CI_indexwise)
colnames(CI_indexwise) = c("i","low","up")
# Compute intersection of all intervals
intersection_low <- max(CI_indexwise$low)
intersection_up <- min(CI_indexwise$up)

# Extend -Inf to a fixed far-left value for display
display_min <- min(CI_indexwise$low[is.finite(CI_indexwise$low)], na.rm = TRUE) - 3

# Modify data to handle -Inf visually
CI_plot <- CI_indexwise %>%
  mutate(
    low_plot = ifelse(is.infinite(low) & low < 0, display_min, low),
    mark_inf = ifelse(is.infinite(low) & low < 0, TRUE, FALSE)
  )

ggplot(CI_plot, aes(y = i, xmin = low_plot, xmax = up)) +
  geom_errorbarh(height = 0.4, color = "blue") +
  geom_vline(xintercept = intersection_low, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = intersection_up, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = lambda, color = "darkgreen", size = 1) +
  
  annotate("rect", xmin = intersection_low, xmax = intersection_up,
           ymin = min(CI_plot$i) - 0.5, ymax = max(CI_plot$i) + 0.5,
           alpha = 0.1, fill = "red") +
  
  geom_point(
    data = subset(CI_plot, mark_inf),
    aes(x = low_plot, y = i),
    shape = 4, color = "red", size = 3, stroke = 1.2
  ) +
  
  labs(
    title = "Confidence interval per equation",
    x = "Interval",
    y = "Index (1 to 30)"
  ) +
  theme_minimal(base_size = 13)+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )







#------------------------------------------------------
## Tukey quantile
q = function(p, lambda) {
  if (abs(lambda) <= 1e-10) # if lambda = 0
    return(log(p / (1 - p)))
  else
    return(1 / lambda * (p^(lambda) - (1 - p)^(lambda)))
}

#--------------------------------------------------------
## Our method of CI for lambda
lambda_ci <- function(L, U, samp, N = 1000){
  n <- length(samp)
  sorted_samp <- sort(samp)
  neg_idx <- which(sorted_samp < 0)   
  pos_idx <- which(sorted_samp > 0)   
  if(sum(L[neg_idx] > 0.5) > 0) return(NULL)
  else if(sum(U[pos_idx] < 0.5) > 0) return(NULL)
  else
  {  neg_L_indices <- neg_idx[ L[neg_idx] < 0.5 ]
  neg_U_indices <- neg_idx[ U[neg_idx] < 0.5 ]
  pos_L_indices <- pos_idx[ L[pos_idx] >= 0.5 ]
  pos_U_indices <- pos_idx[ U[pos_idx] >= 0.5 ]
  alpha_vals_1 <- numeric(0)
  beta_vals_1  <- numeric(0)
  alpha_vals_2 <- numeric(0)
  beta_vals_2  <- numeric(0)
  # Case 1: negative side, L[i] < 1/2 ⇒ Q(L[i],λ) = X_(i)  ⇒ λ ≤ α_i 
  for(i in neg_L_indices){
    f <- function(lambda) q(L[i], lambda) - sorted_samp[i]
    α_i <- uniroot(f, lower = -N, upper = N, tol = 1e-10)$root
    alpha_vals_1 <- c(alpha_vals_1, α_i)
  }
  # Case 2: negative side, U[i] < 1/2 ⇒ Q(U[i],λ) = X_(i)  ⇒ λ ≥ β_i
  for(i in neg_U_indices){
    f <- function(lambda) q(U[i], lambda) - sorted_samp[i]
    β_i <- uniroot(f, lower = -N, upper = N, tol = 1e-10)$root
    beta_vals_1 <- c(beta_vals_1, β_i)
  }
  # Case 3: positive side, L[i] ≥ 1/2 ⇒ Q(L[i],λ) = X_(i)  ⇒ λ ≥ α_i
  for(i in pos_L_indices){
    f <- function(lambda) q(L[i], lambda) - sorted_samp[i]
    α_i <- uniroot(f, lower = -N, upper = N, tol = 1e-10)$root
    alpha_vals_2 <- c(alpha_vals_2, α_i)
  }
  # Case 4: positive side, U[i] ≥ 1/2 ⇒ Q(U[i],λ) = X_(i)  ⇒ λ ≤ β_i
  for(i in pos_U_indices){
    f <- function(lambda) q(U[i], lambda) - sorted_samp[i]
    β_i <- uniroot(f, lower = -N, upper = N, tol = 1e-10)$root
    beta_vals_2 <- c(beta_vals_2, β_i)
  }
  a <- max(beta_vals_1, alpha_vals_2) 
  b <- min(alpha_vals_1, beta_vals_2)
  return(c(a, b))
  }
}



q_abs = function(p,lambda){
  if(abs(lambda) <= 1e-10) return(log((1+p)/(1-p)))
  else return(((1+p)^(lambda)-(1-p)^(lambda))/(lambda*2^(lambda)))
}

lambda_ci_abs = function(L,U,samp,N=1000){
  n = length(samp)
  sorted_samp_abs <- sort(abs(samp))
  alpha = numeric(0)
  beta = numeric(0)
  for(i in 1:n){
    f1 <- function(lambda) {q_abs(L[i], lambda) - sorted_samp_abs[i]}
    alpha_i <- uniroot(f1, lower = -N, upper = N, tol = 1e-10)$root
    alpha = c(alpha,alpha_i)
    f2 <- function(lambda) {q_abs(U[i], lambda) - sorted_samp_abs[i]}
    beta_i <- uniroot(f2, lower = -N, upper = N, tol = 1e-10)$root
    beta = c(beta,beta_i)
  }
  a = max(alpha)
  b = min(beta)
  return(c(a, b))
}



#----------------------------------------------
## Coverage and Width for our method
coverage_lambda_dw = function(lambda,n,R=1000){
  coun = 0
  width = 0
  LU = get_bds_dw(0.05,n)
  L = LU$low
  U = LU$up
  for(r in 1:R){
    samp = rtlambda(n,lambda)
    ci = lambda_ci(L,U,samp)
    a = ci[1]
    b = ci[2]
    if(is.null(a)==TRUE){ 
      coun = coun
      width = width
    }
    else if((a-lambda)*(b-lambda)<0) {coun = coun + 1
    width = width + (b-a)
    }
  }
  return(c(coun/R*100,width/R))
}


coverage_lambda_dw_abs = function(lambda,n,R=1000){
  coun = 0
  width = 0
  LU = get_bds_dw(0.05,n)
  L = LU$low
  U = LU$up
  for(r in 1:R){
    samp = rtlambda(n,lambda)
    ci = lambda_ci_abs(L,U,samp)
    a = ci[1]
    b = ci[2]
    if(is.null(a)==TRUE){ 
      coun = coun
      width = width
    }
    else if((a-lambda)*(b-lambda)<0) {coun = coun + 1
    width = width + (b-a)
    }
  }
  return(c(coun/R*100,width/R))
}





lambda_values <- c(-2,-1,0, 1,2)
n_values <- c(50,100,250,1000)
results1 <- data.frame(lambda = numeric(), n = integer(), coverage = double(), width = double())
results2 <- data.frame(lambda = numeric(), n = integer(), coverage = double(), width = double())
for (lambda in lambda_values) {
  for (n in n_values) {
    res <- coverage_lambda_dw(lambda, n, R = 100)
    results1 <- rbind(results1, data.frame(lambda = lambda, n = n, coverage = res[1], width = res[2]))
  }
}

for (lambda in lambda_values) {
  for (n in n_values) {
    res <- coverage_lambda_dw_abs(lambda, n, R = 100)
    results2 <- rbind(results2, data.frame(lambda = lambda, n = n, coverage = res[1], width = res[2]))
  }
}
results1$source <- "DW"
results2$source <- "DW(abs)"
results_tlambda_dwdwabs <- bind_rows(results1, results2)
results_tlambda_dwdwabs$source <- factor(results_tlambda_dwdwabs$source, levels = c("DW", "DW(abs)"))
results_tlambda_dwdwabs$coverage = results_tlambda_dwdwabs$coverage/100
results_tlambda_dwdwabs_long <- pivot_longer(results_tlambda_dwdwabs, cols = c(coverage, width), names_to = "metric", values_to = "value")
# Relabel 'n' for better column titles
results_tlambda_dwdwabs_long <- results_tlambda_dwdwabs_long %>%
  mutate(n = factor(n, levels = sort(unique(n)), labels = paste0("n = ", sort(unique(n)))))

library(ggh4x)  

ggplot(results_tlambda_dwdwabs_long, aes(x = lambda, y = value, color = source)) +
  geom_line() +
  geom_point() +
  geom_hline(
    data = data.frame(metric = "coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") +
  labs(
    title = "Comparison of methods", 
    x = expression(lambda), 
    y = "Value", 
    color = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) +
  facetted_pos_scales(
    y = list(
      metric == "coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "coverage" ~ scale_y_continuous()  # Default free scale
    )
  )


# 
# n = 30
# LU = get_bds_dw(0.05,n)
# L = LU$low
# U = LU$up
# lambda = 1
# samp = rtlambda(n,lambda)
# suppressWarnings(lambda_ci(L,U,samp,100))
# suppressWarnings(lambda_ci_abs(L,U,samp,100))
# 1/(max(abs(samp))-(1-L[n])/2)
# 1/(max(abs(samp))-(1-U[n])/2)



