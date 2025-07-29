library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)
library(ggpubr)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)



#----------------------------------------------------
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



#-------------------------
CI_quantile = function(samp,L,U,alpha,u, sort = F){
  
  if(sort){
    sorted_samp = samp
  }else{
    sorted_samp = sort(samp)
  }
  
  if(u < tail(L,1) & u > U[1] ){
    u1 = which(U >= u)[1]
    l1 = tail(which((1-L) >= 1-u),1)
    return(c(sorted_samp[u1],sorted_samp[l1]))
  }
  else if(u >= tail(L,1)){
    u1 = which(U >= u)[1]
    return(c(sorted_samp[u1],Inf))
  }
  else if(L[1] <= u  & u <= U[1]){
    l1 = tail(which((1-L) >= 1-u),1)
    return(c(-Inf,sorted_samp[l1]))
  }
  else if(L[1] > u){
    return(c(-Inf,sorted_samp[1]))
  }
}


CI_iqr_u1u2 = function(samp,L,U,alpha,u1,u2, sort= F){
  
  if(sort){
    sorted_samp = samp
  }else{
    sorted_samp = sort(samp)
  }
  
  if(u1 > u2) return("Please choose u1 < u2")
  else if(u1 < u2){
    ci1 = CI_quantile(sorted_samp,L,U,alpha,u1, T)
    ci2 = CI_quantile(sorted_samp,L,U,alpha,u2, T)
    l1 = ci1[1]
    u1 = ci1[2]
    l2 = ci2[1]
    u2 = ci2[2]
    l = max((l2-u1),0)
    u = u2 - l1
    return(c(l,u))
  }
  else return(c(0,0))
}


s_u1u2 = function(chi,xi,u1,u2){
  n = S(u2,chi,xi) - S(u1,chi,xi)
  d = S(3/4,chi,xi) - S(1/4,chi,xi)
  return(n/d)
}


CI_s_u1u2 = function(samp,L,U,alpha,u1,u2, sort = F){

  if(sort){
    sorted_samp = samp
  }else{
    sorted_samp = sort(samp)
  }
  
  if(u1 > u2) return("Please choose u1 < u2")
  else if(u1 < u2){
    CI_n = CI_iqr_u1u2(sorted_samp,L,U,alpha,u1= u1,u2 = u2, T)
    CI_dn = CI_iqr_u1u2(sorted_samp,L,U,alpha,u1 = 1/4,u2 = 3/4, T)
    u = CI_n[2]/CI_dn[1]
    l = CI_n[1]/CI_dn[2]
    return(c(l,u))
  }
}

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

chixi_CI = function(samp,L,U,cuts_comb){
  s_ci_all = matrix(0,nrow(cuts_comb),2)
  for(i in 1:nrow(cuts_comb)){
    s_ci_all[i,] = CI_s_u1u2(samp, L, U, 0.05, cuts_comb[i,1], cuts_comb[i,2])
  }
  
  chi_seq = seq((-1+1/100),(1-1/100), length = 100)
  xi_seq = seq(1/100, (1-1/100), length = 100)
  
  chixi_admiss = matrix(nrow = 0, ncol = 2)
  cut_num = nrow(cuts_comb)
  for (chi in chi_seq) {
    for (xi in xi_seq) {
      s_all = matrix(0,cut_num,1)
      for(i in 1:cut_num){
        s_all[i,] = s_u1u2(chi,xi,cuts_comb[i,1],cuts_comb[i,2])
      }
      s_ind = matrix(0,cut_num,1)
      for(i in 1:cut_num){
        s_ind[i,] = (s_all[i,]-s_ci_all[i,1]) * (s_all[i,]-s_ci_all[i,2]) < 0
      }
      if (sum(s_ind) == cut_num) {
        chixi_admiss = rbind(chixi_admiss, c(chi, xi))
      }
    }
  }
  
  
  chixi_admiss_df = as.data.frame(chixi_admiss)
  names(chixi_admiss_df) = c("chi", "xi")
  param_grid = expand.grid(
    chi = chi_seq,
    xi = xi_seq
  )
  param_grid = param_grid %>%
    mutate(
      admissible = vapply(
        seq_len(n()),
        function(i) {
          c = chi[i]
          x = xi[i]
          any(abs(chixi_admiss_df$chi - c) < 1e-8 & abs(chixi_admiss_df$xi - x) < 1e-8)
        },
        logical(1)
      )
    )
  
  highlight_point <- data.frame(
    chi = 0,
    xi =  0.3661
  )
  
  p = ggplot(param_grid, aes(x = chi, y = xi)) +
    geom_raster(aes(fill = factor(admissible))) +
    geom_point(
      data = highlight_point,
      aes(x = chi, y = xi),
      shape = 21,
      fill = "gold",
      color = "black",
      size = 5,
      stroke = 1.5
    ) +
    scale_fill_manual(values = c("blue", "red"), guide = "none") +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      aspect.ratio = 1
    )+
    ggtitle("Chi Xi CI")
  
  q = ggplot(cuts_comb, aes(x = cuts_comb[,1], y = cuts_comb[,2])) +
    geom_point(color = "blue", size = 3) +
    xlab(names(cuts_comb)[1]) +
    ylab(names(cuts_comb)[2]) +
    ggtitle("Cuts Comb") +
    xlim(0, 1) +
    ylim(0, 1)+
    theme(
      aspect.ratio = 1
    )
  ggarrange(p, q, ncol = 2, nrow = 1)
  
}

###########################################################

n = 500
LU = get_bds_dw(0.05, n)
L = LU$low
U = LU$up


#samp = rgldcsw(n, 1/2, 1/2, 0, (1/2 - 1 / sqrt(5))) # uniform
samp = rgldcsw(n, 0, 1, 0, 0.3661) # normal 

#both

cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts,cuts)
cuts_comb_2 = c(1-cuts,1-rev(cuts))
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)
#p_both = chixi_CI(samp,cuts_comb)
chixi_CI(samp,L,U,cuts_comb)

#slope negative 1
cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts)
cuts_comb_2 = c(1-cuts)
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)
#p_slopeneg1 = chixi_CI(samp,L,U,cuts_comb)
chixi_CI(samp,L,U,cuts_comb)


#slope 1
cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts)
cuts_comb_2 = c(1-rev(cuts))
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)
#p_slope1 = chixi_CI(samp,L,U,cuts_comb)
chixi_CI(samp,L,U,cuts_comb)

#slope 1 corner
cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts)
cuts_comb_2 = c(1-rev(cuts))
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)[c(1:2,19:20),]

#p_slope1_corner = chixi_CI(samp,L,U,cuts_comb)
chixi_CI(samp,L,U,cuts_comb)
cuts_comb = rbind.data.frame(c(0.002,0.712),c(0.028,0.738),c(0.210,0.920),c(0.236,0.946))
cuts_comb = rbind.data.frame(c(0.002,0.712),c(0.028,0.738),c(0.210,0.920),c(0.236,0.946))

#slope 1 mid
cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts)
cuts_comb_2 = c(1-rev(cuts))
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)[c(2:19),]
#p_slope1_mid =chixi_CI(samp,L,U,cuts_comb)
chixi_CI(samp,L,U,cuts_comb)


#slope 1 up corner
cuts_comb = rbind.data.frame(c(0.002,0.750),c(0.250,0.998))
p_slope1_up_corner = chixi_CI(samp,L,U,cuts_comb)


#slope 1 upper corner
cuts_comb = rbind.data.frame(c(0.002,0.750),c(0.250,0.998))
p_slope1_up_corner = chixi_CI(samp,L,U,cuts_comb)

#slope 1 low corner
cuts_comb = rbind.data.frame(c(0.002,0.874),c(0.126,0.998))
p_slope1_upper_corner = chixi_CI(samp,L,U,cuts_comb)

#slope 1 lower corner
cuts_comb = cbind.data.frame(c(0.002,0.750),c(0.250,0.998))
p_slope1_lower_corner = chixi_CI(samp,L,U,cuts_comb)



cuts_comb = rbind.data.frame(c(10^(-6),0.5),
                             c(0.5, 1 - 10^(-6)),
                             c(10^(-6),0.3),
                             c(0.7, 1 - 10^(-6)),
                             # c(0.7, 1 - 10^(-6)),
                             # c(1 - 10^(-6), 0.7),
                             # c(10^(-6),0.3),
                             # c(0.7, 1 - 10^(-6)),
                             c(0.1,0.3),
                             c(0.7, 1 - 0.1)
)
colnames(cuts_comb) = c("i","j")
chixi_CI(samp,L,U,cuts_comb)

########################################

k = 10
seq1 = floor(seq(1,(n-1),length = k))/n

cuts_comb = NULL

for( i in 2:(k-1)){
  cuts_comb = rbind(cuts_comb, c(seq1[1],seq1[i]))
  cuts_comb = rbind(cuts_comb, c(seq1[i], seq1[k]))
}

for( i in 3:(k-2)){
  cuts_comb = rbind(cuts_comb, c(seq1[2],seq1[i]))
  cuts_comb = rbind(cuts_comb, c(seq1[i], seq1[k-1]))
}

cuts_comb = as.data.frame(cuts_comb)
chixi_CI(samp,L,U,cuts_comb)


########################################
library(dplyr)
library(MASS) 
library(reshape2) 
library(reshape) 


############################################

#samp = rgldcsw(n, 0, 0.25, 0.2844,  0.3583)

k = ceiling(sqrt(n))
seq1 = seq(10^(-6), 1 - 10^(-6),length = k)
cuts_comb = t(combn(seq1, m = 2)) %>% as.data.frame()
chixi_CI(samp,L,U,cuts_comb)


##############################################

cuts = floor(seq(1,(n/2-1),length = 20))/n
cuts_comb_1 = c(cuts,cuts)
cuts_comb_2 = c((1-cuts),(1-rev(cuts)))
cuts_comb_3 = cbind(cuts_comb_1,cuts_comb_2)

cuts_comb_4 = rbind.data.frame(c(10^(-6),0.5),
                               c(0.5, 1 - 10^(-6)),
                               c(10^(-6),0.3),
                               c(0.7, 1 - 10^(-6)),
                               c(0.1,0.3),
                               c(0.7, 1 - 0.1)
)
colnames(cuts_comb_4) = colnames(cuts_comb_3)

cuts_comb = rbind.data.frame(cuts_comb_3,cuts_comb_4)
chixi_CI(samp,L,U,cuts_comb)

############ 

cuts = sort(c(L,U)) 
k = 100
cuts_sparse = cuts[floor(seq(2, length(cuts)-1, length.out = k))]
cuts_comb = t(combn(cuts_sparse, m = 2)) %>% as.data.frame()

chixi_CI(samp,L,U,cuts_comb)
