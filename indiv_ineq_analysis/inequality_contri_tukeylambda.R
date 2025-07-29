library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)
library(ggpubr)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)


############################################################################################################

## DW CI
KL_dw = function(p,q){
  if(length(p)!= length(q))
  {
    print("p and q length not same")
    stop()
  }
  
  val = numeric(length(p))
  
  for( i in 1:length(p))
  {
    val[i] = p[i]*log(p[i]/q[i]) + (1-p[i])*log((1-p[i])/(1-q[i]))
    if( p[i] == 0 ){val[i] = log(1/(1-q[i]))}
    if( p[i] == 1 ){val[i] = log(1/q[i]) }
    if( q[i] == 0 ){val[i] = Inf }
    if( q[i] == 1 ){val[i] = Inf}
    
    if( p[i] == 0 & q[i] == 0 ){val[i] = 0}
    if( p[i] == 1 & q[i] == 0 ){val[i] = Inf }
    if( p[i] == 0 & q[i] == 1 ){val[i] = Inf }
    if( p[i] == 1 & q[i] == 1 ){val[i] = 0}
  }
  
  return(val)
}

#########################
C_dw = function(t)
{
  return( log( 1 - log(2*t) - log(2*(1-t)) ) )
}

#########################
D_dw = function(t)
{
  return( log(1 + (C_dw(t))^2) )
}

C_v_dw = function(u,t, v)
{
  len_t = length(t)
  if(length(u)!= len_t)
  {
    print("not same length")
    stop
  }
  
  temp = numeric(len_t)
  for( i in 1:len_t)
  {
    if(max(u[i],t[i]) < 1/2){
      temp[i] = C_dw(max(u[i],t[i])) + (v*D_dw(max(u[i],t[i])))
    } else if(min(u[i],t[i]) > 1/2){
      temp[i] = C_dw(min(u[i],t[i])) + (v*D_dw(min(u[i],t[i])))
    } else {
      temp[i] = 0
    }
  }
  return(temp)
}

#########################

dw_quant = function( v = 3/2, alpha = 0.05, n_runs = 10000, n )
{
  sam = numeric(n_runs)
  t_n = (1:n)/(n) 
  t_n_shift = (0:(n-1))/n
  for( i in 1:n_runs)
  {
    U_sort = sort( runif(n))
    
    sam[i] = max( n*KL_dw(t_n_shift, U_sort) - C_v_dw(t_n_shift, U_sort,v = v),  
                  n*KL_dw(t_n, U_sort) - C_v_dw(t_n, U_sort,  v = v) ) 
  }
  
  n_cuttoff = ceiling(n_runs*(1-alpha))
  
  return(sort(sam)[n_cuttoff])
}

fin_func_dw = function(p,q,quant,v = 3/2,n)
{
  return(KL_dw(p,q) - ((C_v_dw(p,q,v = v)+ quant)/n)  )
}

get_bds_dw = function(alpha = 0.05, n, v = 3/2, quant_val = NA)
{
  if(is.na(quant_val))
  { 
    dw_quant_val = dw_quant(alpha = alpha, n = n, v = v)
  }else{
    dw_quant_val = quant_val
  }
  
  bds_dw = data.frame(low = rep(NA, n),
                      up = rep(NA, n))
  
  t_n = (1:n)/n
  #t_n_shift = (0:(n-1))/n
  
  for(i in 1:(n-1))
  {
    ### upper bound
    if( fin_func_dw(p = t_n[i], q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
      temp_up = 1
    }else{
      temp_up = uniroot( fin_func_dw, lower = t_n[i] , upper = 1-10^(-10),
                         p = t_n[i], quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
    }
    
    bds_dw$up[i]  = temp_up
    bds_dw$low[n-i] = 1 - temp_up
  }
  
  #### upper for i = n
  bds_dw$up[n] = 1
  
  ### upper bound for i = 0
  if( fin_func_dw(p = 0, q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
    temp_up_0 = 1
  }else{
    temp_up_0 = uniroot( fin_func_dw, lower = 0 , upper = 1-10^(-10),
                         p = 0, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
  }
  
  ##### upper for i = 0
  bds_dw$low[n] = 1 - temp_up_0
  
  return(bds_dw)
}

############################################################################################################
############################################################################################################


n = 30
LU = get_bds_dw(0.05,n)
L = LU$low
U = LU$up
lambda = -2
samp = rtlambda(n,lambda)

sorted_samp <- sort(samp)
neg_idx <- which(sorted_samp < 0)   
pos_idx <- which(sorted_samp > 0) 
CI_indexwise = matrix(0,nrow = n, ncol = 3)
CI_indexwise[,1] = c(1:n)
N = 100


for(i in neg_idx){
  if(L[i] > 0.5) {CI_indexwise[i,2] = CI_indexwise[i,3] = lambda}
  else if(U[i] >= 0.5) {
    f1 <- function(lambda) qtlambda(L[i], lambda) - sorted_samp[i]
    alpha_i <- uniroot(f1, lower = -N, upper = N, tol = 1e-10)$root
    CI_indexwise[i,2] = -Inf
    CI_indexwise[i,3] = alpha_i
  }
  else{
    f1 <- function(lambda) qtlambda(L[i], lambda) - sorted_samp[i]
    alpha_i <- uniroot(f1, lower = -N, upper = N, tol = 1e-10)$root
    f2 <- function(lambda) qtlambda(U[i], lambda) - sorted_samp[i]
    beta_i <- uniroot(f2, lower = -N, upper = N, tol = 1e-10)$root
    CI_indexwise[i,2] = beta_i
    CI_indexwise[i,3] = alpha_i
  }
}


for(i in pos_idx){
  if(U[i] < 0.5) {CI_indexwise[i,2]=CI_indexwise[i,3] = lambda}
  else if(L[i] <= 0.5) {
    f2 <- function(lambda) qtlambda(U[i], lambda) - sorted_samp[i]
    beta_i <- uniroot(f2, lower = -N, upper = N, tol = 1e-10)$root
    CI_indexwise[i,3] = beta_i
    CI_indexwise[i,2] = -Inf
  }
  else{
    f1 <- function(lambda) qtlambda(L[i], lambda) - sorted_samp[i]
    alpha_i <- uniroot(f1, lower = -N, upper = N, tol = 1e-10)$root
    f2 <- function(lambda) qtlambda(U[i], lambda) - sorted_samp[i]
    beta_i <- uniroot(f2, lower = -N, upper = N, tol = 1e-10)$root
    CI_indexwise[i,2] = alpha_i
    CI_indexwise[i,3] = beta_i
  }
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
  geom_vline(xintercept = intersection_low, linetype = "dashed", color = "red", size = 0.5) +
  geom_vline(xintercept = intersection_up, linetype = "dashed", color = "red", size = 0.5) +
  geom_vline(xintercept = lambda, color = "darkgreen", size = 0.5) +

  annotate("rect", xmin = intersection_low, xmax = intersection_up,
           ymin = min(CI_plot$i) - 0.5, ymax = max(CI_plot$i) + 0.5,
           alpha = 0.1, fill = "red") +

  geom_point(
    data = subset(CI_plot, mark_inf),
    aes(x = low_plot, y = i),
    shape = 4, color = "red", size = 3, stroke = 1.2
  ) +
  
  labs(
    title = bquote( lambda == .(lambda)),
    x = "Interval",
    y = "Index (1 to 30)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

setwd("~/Documents/Projects/GLD_inference")

ggsave( paste("tl_ineq_lam", lambda, ".pdf", sep = ""), width = 7, height = 5, dpi = 300, units = "in" )

