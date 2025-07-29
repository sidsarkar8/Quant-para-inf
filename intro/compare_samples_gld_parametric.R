
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



n = 10000

grid1 = seq(-3, 10, length.out = 100)


par(mfrow = c(1,2))

###### 

samp_gld1 = rgldcsw(n, 1, 0.28, 0.2844,  0.3583)  # log normal
samp_p1 = exp(rnorm(n, sd = 0.25))

plot(grid1, ecdf(samp_gld1)(grid1), type = "l")
points(grid1, ecdf(samp_p1)(grid1), col = "red", type = "l")

qqplot(samp_gld1, samp_p1)
abline(a = 0 , b = 1)

###### 

samp_gld2 = rgldcsw(n, 0, 1, 0, 0.3661) # normal 
samp_p2 = rnorm(n)

plot(grid1, ecdf(samp_gld2)(grid1), type = "l")
points(grid1, ecdf(samp_p2)(grid1), col = "red", type = "l")

qqplot(samp_gld2, samp_p2)
abline(a = 0 , b = 1)


######
samp_gld3 = rgldcsw(n, 1/2, 1/2, 0, (1/2 - 1 / sqrt(5))) # uniform
samp_p3 = runif(n, min = -1, max = 1 )

plot(grid1, ecdf(samp_gld3)(grid1), type = "l")
points(grid1, ecdf(samp_p3)(grid1), col = "red", type = "l")

qqplot(samp_gld1, samp_p1)
abline(a = 0 , b = 1)

