library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)


#-------------------------

CI_quantile = function(samp,L,U,alpha,u){
  sorted_samp = sort(samp)
  u1 = which(U >= u)[1]
  l1 = tail(which((1-L) >= 1-u),1)
  return(c(sorted_samp[u1],sorted_samp[l1]))
}


CI_iqr_u1u2 = function(samp,L,U,alpha,u1,u2){
  if(u1 > u2) return("Please choose u1 < u2")
  else if(u1 < u2){
    ci1 = CI_quantile(samp,L,U,alpha,u1)
    ci2 = CI_quantile(samp,L,U,alpha,u2)
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


CI_s_u1u2 = function(samp,L,U,alpha,u1,u2){
  if(u1 > u2) return("Please choose u1 < u2")
  else if(u1 < u2){
    CI_n = CI_iqr_u1u2(samp,L,U,alpha,u1,u2)
    CI_dn = CI_iqr_u1u2(samp,L,U,alpha,1/4,3/4)
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

n = 500
LU = get_bds_dw(0.05, n)
L = LU$low
U = LU$up


#samp = rgldcsw(n, 0, 0.25, 0.2844,  0.3583)  # log normal
samp = rgldcsw(n, 0, 1, 0, 0.3661) # normal 
#samp = rgldcsw(n, 1/2, 1/2, 0, (1/2 - 1 / sqrt(5))) # uniform

# cut 1
cuts = floor(seq(1,(n/2-1),length = 5))/n
cuts_comb_1 = c(cuts,cuts)
cuts_comb_2 = c(1-cuts,1-rev(cuts))
cuts_comb = cbind(cuts_comb_1,cuts_comb_2)

# cut 2
# cuts = c(0.1, 0.25, 0.5, 0.75, 0.9)
# cuts_comb_1 = cuts[combn(1:5,2)[1,]]
# cuts_comb_2 = cuts[combn(1:5,2)[2,]]
# cuts_comb = cbind(cuts_comb_1,cuts_comb_2)

s_ci_all = matrix(0,10,2)
for(i in 1:10){
  s_ci_all[i,] = CI_s_u1u2(samp, L, U, 0.05, cuts_comb[i,1], cuts_comb[i,2])
}

chi_seq = seq((-1+1/100),(1-1/100), length = 100)
xi_seq = seq(1/100, (1-1/100), length = 100)

## Store admissible pairs
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

## Build the full grid of sampled (chi, xi)
param_grid = expand.grid(
  chi = chi_seq,
  xi = xi_seq
)

## Mark admissible vs non-admissible
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

# highlight_point <- data.frame(
#   chi = 0.2844,
#   xi =  0.3583
# )

highlight_point <- data.frame(
  chi = 0,
  xi =  0.3661
)

# highlight_point <- data.frame(
#   chi = 0,
#   xi =  (1/2 - 1 / sqrt(5))
# )

ggplot(param_grid, aes(x = chi, y = xi)) +
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
  )







