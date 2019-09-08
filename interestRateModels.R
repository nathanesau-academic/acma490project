## Acma 490 Project
## Nathan Esau
## 
## This script is from the stocins package I made:
##
##    i.e. see https://github.com/nathanesau/stocins
##
## All of these functions are documented in that package
## explaining what they do and examples of how to use the functions.
##
## This script contains functions for the first half of the 490
## course, related to interest rate models.

## Interest Rate Model

iratemodel <- function(params, class)
{
  class(params) = append(class(params),"iratemodel")
  
  if(class == "ou")
  {
    class(params) = append(class(params),"ou")
    
    return(params)
  }
  else if(class == "wiener")
  {
    class(params) = append(class(params),"wiener")
    
    return(params)
  }
  else if(class == "determ")
  {
    class(params) = append(class(params),"determ")
    
    return(params)
  }
  else if(class == "second")
  {
    class(params) = append(class(params),"second")
    
    return(params)
  }
  else if(class == "ar1")
  {
    class(params) = append(class(params),"ar1")  
    
    return(params)
  }
  else if(class == "arma")
  {
    class(params) = append(class(params),"arma")
    
    return(params)
  }
  else 
  {
    stop("unknown class type")
  }
}

iratemodel.convert <- function(from, to, irm, Delta = 1)
{
  if(from == "second")
  {
    if(to == "arma")
    {
      alpha1 <- irm$alpha1
      alpha2 <- irm$alpha2
      sigma2 <- irm$sigma^2
      
      # eigenvalues
      mu1 <- 0.5 * (alpha1 + sqrt(alpha1^2 + 4*alpha2))
      mu2 <- 0.5 * (alpha1 - sqrt(alpha1^2 + 4*alpha2))
      
      lambda1 <- exp(mu1 * Delta)
      lambda2 <- exp(mu2 * Delta)
      
      phi1 <- lambda1 + lambda2
      phi2 <- -lambda1 * lambda2
      
      P = (-mu1*(1+lambda1^2)*(1-lambda2^2) + mu2*(1+lambda2^2)*(1-lambda1^2)) /
        (mu1*lambda1*(1-lambda2^2) - mu2*lambda2*(1-lambda1^2))
      
      theta1 <- Re(polyroot(c(1,P,1)))
      theta1 <- theta1[which(abs(theta1) == min(abs(theta1)))]
      
      sigma2a = sigma2 * (lambda1 - lambda2)^2 / (2*mu1*(mu1^2 - mu2^2) * (lambda1 - theta1)) /
        ((lambda1 - theta1)/(1 - lambda1^2) - (lambda2 - theta1)/(1 - lambda1*lambda2))
      
      armamodel = iratemodel(params = list(phi1 = phi1, phi2 = phi2,
                                           theta1 = theta1,
                                           sigma = sqrt(sigma2a)), "arma")
      
      return(armamodel)
    }
    else
    {
      stop(paste("conversion from", from, "to", to, "not implemented"))
    }
  }
  else if(from == "arma")
  {
    if(to == "second")
    {
      phi1 <- irm$phi1
      phi2 <- irm$phi2
      
      lambda1 <- 0.5 * (phi1 + sqrt(phi1^2 + 4*phi2))
      lambda2 <- 0.5 * (phi1 - sqrt(phi1^2 + 4*phi2))
      
      mu1 <- log(lambda1) / Delta
      mu2 <- log(lambda2) / Delta
      
      alpha1 = mu1 + mu2
      alpha2 = -mu1 * mu2
      
      secondmodel = iratemodel(params = list(alpha1 = alpha1, 
        alpha2 = alpha2, sigma = 0),"second") # sigma arbitrary for now
      
      return(secondmodel)
    }
    else
    {
      stop(paste("conversion from", from, "to", to, "not implemented"))  
    }
  }
  else if(from == "ar1")
  {
    if(to == "ou")
    {
      phi1 = irm$phi1
      sigma2a = irm$sigma^2
      
      alpha = -log(phi1) / Delta
      sigma2 = 2*alpha*sigma2a / (1 - phi1^2)
      
      oumodel = iratemodel(params = list(alpha = alpha, sigma = sqrt(sigma2)),
                        "ou")
    
    return(oumodel)
    }
    else
    {
      stop(paste("conversion from", from, "to", to, "not implemented"))  
    }
  }
  else if(from == "ou")
  {
    if(to == "ar1")
    {
      alpha = irm$alpha
      sigma2 = irm$sigma^2
      
      phi1 = exp(-alpha*Delta)
      sigma2a = sigma2 / (2*alpha) * (1 - phi1^2)
      
      armodel = iratemodel(params = list(coef = phi1, sigma = sqrt(sigma2a)),
                           "ar1")
      
      return(armodel)
    }
    else
    {
      stop(paste("conversion from", from, "to", to, "not implemented"))  
    }
  }
  else
  {
    stop(paste("conversion from", from, "to", to, "not implemented"))
  }    
}

delta.ev <- function(t, irm)
{
  UseMethod("delta.ev", irm)
}

delta.cov <- function(s, t, irm)
{
  UseMethod("delta.cov", irm)
}

delta.var <- function(t, irm)
{
  UseMethod("delta.var", irm)
}

y.ev <- function(t, irm)
{
  UseMethod("y.ev", irm)
}

y.var <- function(t, irm)
{
  UseMethod("y.var", irm)
}

y.cov <- function(s, t, irm)
{
  UseMethod("y.cov", irm)
}

pv.moment <- function(t, moment, irm)
{
  mu = -y.ev(t, irm)
  sigma2 = y.var(t, irm)
  
  exp(mu * moment + 0.5 * sigma2 * moment^2)
}

pv.ev <- function(t, irm)
{
  pv.moment(t, moment = 1, irm)
}

pv.var <- function(t, irm)
{
  pv.moment(t, moment = 2, irm) - pv.moment(t, moment = 1, irm)^2
}

pv.cov <- function(s, t, irm)
{
  mu <- -(y.ev(s, irm) + y.ev(t, irm))
  sigma2 <- y.var(s,irm) + y.var(t,irm) + 2*y.cov(s,t,irm)
  EXY <- exp(mu + 0.5*sigma2)
  EXEY <- pv.ev(s, irm) * pv.ev(t, irm)
  
  EXY - EXEY
}

## Deterministic Interest Rate

delta.ev.determ <- function(t, irm)
{
  irm$delta * t^0 # constant
}

delta.cov.determ <- function(s, t, irm)
{
  0 * s * t # no covariance
}

delta.var.determ <- function(t, irm)
{
  0 * t # no variance
}

y.ev.determ <- function(t, irm)
{
  irm$delta * t
}

y.var.determ <- function(t, irm)
{
  0 * t # no variance
}

y.cov.determ <- function(s, t, irm)
{
  0 * t # no covariance
}

## Wiener Process

delta.ev.wiener <- function(t, irm)
{
  irm$delta * t^0 # doesn't depend on t
}

delta.cov.wiener <- function(s, t, irm)
{
  irm$sigma^2 * min(s, t)
}

delta.var.wiener <- function(t, irm)
{
  irm$sigma^2 * t
}

y.ev.wiener <- function(t, irm)
{
  irm$delta * t
}

y.cov.wiener <- function(s, t, irm)
{
  irm$sigma^2 * (min(s,t)^2 * max(s, t) / 2 - min(s,t)^3 / 6)
}

y.var.wiener <- function(t, irm)
{
  irm$sigma^2 * t^3 / 3
}

## Ornstein-Uhlenbeck Process

delta.ev.ou <- function(t, irm)
{
  irm$delta + exp(-irm$alpha * t) * (irm$delta0 - irm$delta)
}

delta.cov.ou <- function(s, t, irm)
{
  stop("delta.cov not implemented for ou")
}

delta.var.ou <- function(t, irm)
{
  irm$sigma^2 / (2*irm$alpha) * (1 - exp(-2*irm$alpha*t))
}

y.ev.ou <- function(t, irm)
{
  irm$delta*t + (irm$delta0 - irm$delta) * (1 - exp(-irm$alpha*t)) / irm$alpha
}

y.var.ou <- function(t, irm)
{
  irm$sigma^2/irm$alpha^2*t + irm$sigma^2/(2*irm$alpha^3) * 
    (-3 + 4*exp(-irm$alpha*t) - exp(-2*irm$alpha*t))
}

y.cov.ou <- function(s, t, irm)
{
  irm$sigma^2 / irm$alpha^2 * min(s, t) +
    irm$sigma^2 / (2*irm$alpha^3) * (-2 + 2*exp(-irm$alpha*s) + 2*exp(-irm$alpha*t) 
                             - exp(-irm$alpha*(abs(t-s))) - exp(-irm$alpha*(t+s)))
}

## Second Order Stochcastic Differential Equation

delta.ev.second <- function(t, irm)
{
  f <- function(s)
  {
    lambda1 = (irm$alpha1 - sqrt(irm$alpha1^2 + 4*irm$alpha2)) / 2
    lambda2 = (irm$alpha1 + sqrt(irm$alpha1^2 + 4*irm$alpha2)) / 2
    
    A11 = (lambda2 * exp(lambda2 * s) - lambda1 * exp(lambda1 * s)) / (lambda2 - lambda1)
    A12 = lambda1 * lambda2 * (exp(lambda1 * s) - exp(lambda2 * s)) / (lambda2 - lambda1)
    A21 = (exp(lambda2 * s) - exp(lambda1 * s)) / (lambda2 - lambda1)
    A22 = (lambda2 * exp(lambda1 * s) - lambda1 * exp(lambda2 * s)) / (lambda2 - lambda1)
    A = matrix(c(A11, A21, A12, A22), 2, 2)
    
    out = A %*% matrix(c(irm$delta0prime, irm$delta0 - irm$delta), 2, 1)
    out[2,1] + irm$delta
  }
  
  sapply(t, f)
}

delta.cov.second <- function(s, t, irm)
{
  stop("delta.cov is not implemented for second order sde")
}

delta.var.second <- function(t, irm)
{
  f <- function(s)
  {
    lambda1 = (irm$alpha1 - sqrt(irm$alpha1^2 + 4*irm$alpha2)) / 2
    lambda2 = (irm$alpha1 + sqrt(irm$alpha1^2 + 4*irm$alpha2)) / 2
    
    A11 = (lambda2 * exp(lambda2 * s) - lambda1 * exp(lambda1 * s)) / (lambda2 - lambda1)
    A12 = lambda1 * lambda2 * (exp(lambda1 * s) - exp(lambda2 * s)) / (lambda2 - lambda1)
    A21 = (exp(lambda2 * s) - exp(lambda1 * s)) / (lambda2 - lambda1)
    A22 = (lambda2 * exp(lambda1 * s) - lambda1 * exp(lambda2 * s)) / (lambda2 - lambda1)
    A = matrix(c(A11, A21, A12, A22), 2, 2)
    
    I11 = irm$sigma^2 / (lambda2 - lambda1)^2 *
        (-lambda2/2 * (exp(-2*lambda2*s) - 1) +
          2 * lambda1 * lambda2 / (lambda1 + lambda2) * (exp(-(lambda1 + lambda2)*s) - 1) +
          (-lambda1 / 2) * (exp(-2*lambda1*s) - 1))
    
    I12 = irm$sigma^2 / (lambda2 - lambda1)^2 *
      ( -1/2 * (exp(-2*lambda2*s) - 1) +
          (exp(-(lambda1 + lambda2)*s) - 1) +
          (-1/2) * (exp(-2*lambda1*s) - 1)
      )
    
    I21 = I12
    
    I22 = irm$sigma^2 / (lambda2 - lambda1)^2 *
      ( -1/(2*lambda2) * (exp(-2*lambda2*s) - 1) +
          2/(lambda1 + lambda2) * (exp(-(lambda1 + lambda2)*s) - 1) +
          (-1/(2*lambda1)) * (exp(-2*lambda1 * s) - 1)
      )
    
    I = matrix(c(I11, I21, I12, I22), 2, 2)
    
    out = A %*% I %*% t(A)
    out[2,2]
  }
  
  sapply(t, f)
}

y.ev.second <- function(t, irm)
{
  stop("y.ev is not implemented for second order sde")
}

y.var.second <- function(t, irm)
{
  stop("y.var is not implemented for second order sde")
}

y.cov.second <- function(s, t, irm)
{
  stop("y.cov is not implemented for second order sde")
}

## AR(1) process

delta.ev.ar1 <- function(t, irm)
{
  irm$phi1^t * (irm$delta0 - irm$delta) + irm$delta
}

delta.cov.ar1 <- function(s, t, irm)
{
  irm$phi1^(abs(s-t)) * (1 - irm$phi1^(2*min(s,t))) / (1 - irm$phi1^2) * irm$sigma^2
}

delta.var.ar1 <- function(t, irm)
{
  delta.cov.ar1(t, t, irm)
}

y.ev.ar1 <- function(t, irm)
{
  sum(delta.ev(0:(t-1), irm))
}

y.var.ar1 <- function(t, irm)
{
  total <- 0
  for(s in 0:(t-1))
  {
    for(r in 0:(t-1))
    {
      total <- total + delta.cov(s,r,irm)
    }
  }
  total
}

y.cov.ar1 <- function(s, t, irm)
{
  total <- 0
  for(m in 0:(s-1))
  {
    for(r in 0:(t-1))
    {
      total <- total + delta.cov(m,r,irm)
    }
  }
  total
}

## ARMA(2,1) Process

delta.ev.arma <- function(t, irm)
{
  stop("delta.ev not implemented for arma")
}

delta.cov.arma <- function(s, t, irm)
{
  stop("delta.cov not implemented for arma")
}

delta.var.arma <- function(t, irm)
{
  stop("delta.var not implemented for arma")
}

y.ev.arma <- function(t, irm)
{
  stop("y.ev not implemented for arma")
}

y.var.arma <- function(t, irm)
{
  stop("y.var not implemented for arma")
}

y.cov.arma <- function(s, t, irm)
{
  stop("y.cov not implemented for arma")
}

## Annuity Functions

ann.ev <- function(n,irm)
{
  total = 0
  
  for(t in 1:n)
  {
    total = total + pv.ev(t,irm)
  }
  
  total
}

ann.var <- function(n, irm)
{
  total = 0
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      total = total + pv.cov(i, j, irm)
    }
  }
  total
}