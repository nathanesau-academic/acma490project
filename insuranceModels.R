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
## This script contains functions for the second half of the 490
## course, related to life insurance.

## Mortality Assumptions

mortassumptions <- function(params)
{
  # calculate omega from the table
  params$omega = tail(get(params$table), 1)$x + 1
  
  # append class name
  class(params) = append(class(params), "mortassumptions")

  return(params)
}

kpx <- function(k, mort)
{
  UseMethod("kpx", mort)
}

kdeferredqx <- function(k, mort)
{
  UseMethod("kdeferredqx", mort)
}

kpx.mortassumptions <- function(k, mort)
{
  f <- function(k)
  {
    ifelse(k > mort$omega - mort$x, 0, 
    ifelse(k > 0, prod(1 - get(mort$table)$qx[seq(mort$x, mort$x+k-1, 1) + 1]), 1))
  }
  
  sapply(k, f)
}

kdeferredqx.mortassumptions <- function(k, mort)
{
  f <- function(k)
  {
    ifelse(mort$x + k > mort$omega - 1, 0,
    kpx(k, mort) * get(mort$table)$qx[mort$x+k+1])
  }
  
  sapply(k, f)
}

## Insurance Product

insurance <- function(params, class, subclass = NULL)
{
  class(params) = append(class(params),"insurance")
  
  if(class == "isingle")
  {
    class(params) = append(class(params), "isingle")
    
    if(subclass == "term")
    {
      params$e = 0 # make sure e exists in term for class ``igroup``
      
      class(params) = append(class(params), "termsingle")
      
      return(params)
    }
    else if(subclass == "endow")
    {
      params$e = ifelse(is.null(params$e), 0, params$e)
      
      class(params) = append(class(params), "endowsingle")
      
      return(params)
    }
    else
    {
      stop("unknown subclass type in isingle") 
    }
  }
  else if(class == "iport")
  {
    class(params) = append(class(params), "iport")
    
    if(subclass == "term")
    {
      class(params) = append(class(params), "termport")
      
      return(params)
    }
    else if(subclass == "endow")
    {
      class(params) = append(class(params), "endowport")
      
      return(params)
    }
    else
    {
      stop("unknown subclass type in iport")
    }
  }
  else if(class == "igroup")
  {
    c = sum(unlist(lapply(params, function(x) x$c))) 
    n = max(unlist(lapply(params, function(x) x$single$n)))
    
    params$c = c
    params$n = n
    
    class(params) = append(class(params), "igroup")
    
    return(params)
  }
  else
  {
    stop("unknown class type")
  }
}

z.moment <- function(moment, ins, mort, irm)
{
  UseMethod("z.moment", ins)
}

z.ev <- function(ins, mort, irm)
{
  z.moment(1, ins, mort, irm)
}

z.sd <- function(ins, mort, irm)
{
  sqrt(z.moment(2, ins, mort, irm) - z.moment(1, ins, mort, irm)^2)
}

z.sk <- function(ins, mort, irm)
{
  u1 = z.moment(1, ins, mort, irm)
  u2 = z.moment(2, ins, mort, irm)
  u3 = z.moment(3, ins, mort, irm)

  num = u3 - 3*u2*u1 + 2*u1^3
  den = z.sd(ins, mort, irm)^3

  num / den
}

z.insrisk <- function(ins, mort, irm)
{
  UseMethod("z.insrisk", ins)
}

z.invrisk <- function(ins, mort, irm)
{
  UseMethod("z.invrisk", ins)
}

z.pdf <- function(z, ins, mort, irm)
{
  UseMethod("z.pdf", ins)
}

## Single Insurance Product

z.moment.isingle <- function(moment, ins, mort, irm)
{
  UseMethod("z.moment.isingle", ins)
}

z.insrisk.isingle <- function(ins, mort, irm)
{
  UseMethod("z.insrisk.isingle", ins)
}

z.invrisk.isingle <- function(ins, mort, irm)
{
  UseMethod("z.invrisk.isingle", ins)
}

z.pdf.isingle <- function(z, ins, mort, irm)
{
  UseMethod("z.pdf.isingle", ins)
}

## Single Term Insurance Product

z.moment.isingle.termsingle <- function(moment, ins, mort, irm)
{
  total = 0
  
  for(k in seq(0, ins$n - 1, 1))
  {
    total = total + pv.moment(k+1, moment, irm) * ins$d^moment * kdeferredqx(k, mort)
  }
  
  total
}

z.ev.two.isingle.termsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      mu = -(y.ev(i+1,irm) + y.ev(j+1,irm))
      sigma2 = y.var(i+1,irm) + y.var(j+1,irm) + 2*y.cov(i+1,j+1,irm)
      EPV = exp(mu + 0.5* sigma2) # EPV not accounting for mortality
      
      total = total + EPV * kdeferredqx(i, mort) * 
        kdeferredqx(j, mort) * ins$d^2
    }
  }
  
  total
}

z.ev.three.isingle.termsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      for(k in seq(0, ins$n-1, 1))
      {
        mu = -(y.ev(i+1,irm) + y.ev(j+1,irm) + y.ev(k+1,irm))
        sigma2 = y.var(i+1,irm) + y.var(j+1,irm) + y.var(k+1,irm) +
          2 * y.cov(i+1, j+1, irm) + 2 * y.cov(i+1, k+1, irm) + 2 * y.cov(j+1, k+1, irm)
        
        EPV <- exp(mu + 0.5* sigma2) # EPV not accounting for mortality
        
        total = total + EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) *
          kdeferredqx(k, mort) * ins$d^3
      }
    }
  }
  
  total
}

z.ev.twoone.isingle.termsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      mu <- -(2*y.ev(i+1,irm) + y.ev(j+1,irm))
      sigma2 <- 4*y.var(i+1,irm) + y.var(j+1,irm) + 4*y.cov(i+1,j+1,irm)
      EPV <- exp(mu + 0.5* sigma2) # EPV not accounting for mortality
      
      total = total + EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) * ins$d^2 * ins$d
    }
  }
  
  total
}

z.insrisk.isingle.termsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(k in seq(0,ins$n-1,1))
  {
    for(j in seq(0,ins$n-1,1))
    {
      if(k==j)
      {
        total = total + pv.cov(k+1,j+1,irm) * kdeferredqx(k,mort)
      }
    }
  }
  
  total
}

z.invrisk.isingle.termsingle <- function(ins, mort, irm)
{
  second = 0
  first = 0
  
  for(k in seq(0,ins$n-1,1))
  {
    second = second + pv.moment(k+1,1,irm)^2 * kdeferredqx(k, mort)
    first = first + pv.moment(k+1,1,irm) * kdeferredqx(k, mort)
  }
  
  second - first^2
}

z.pdf.isingle.termsingle <- function(z, ins, mort, irm)
{
  f <- function(z)
  {
    if(z < 0)
    {
      return(0)
    }
    else if(z == 0)
    {
      return(kpx(k = ins$n, mort))
    }
    else
    {
      if(ins$n > 0)
      {
        k = seq(0, ins$n - 1, 1)
        kqx = kdeferredqx(k, mort)
          
        sum(kqx * dlnorm(z, -y.ev(k + 1, irm), sqrt(y.var(k + 1, irm))))
      }
      else
      {
        return (0)
      }
    }
  }
    
  sapply(z, f)
}

## Single endowment insurance product

z.moment.isingle.endowsingle <- function(moment, ins, mort, irm)
{
  total = 0
  
  for(k in seq(0, ins$n - 1, 1))
  {
    total = total + pv.moment(k+1, moment, irm) * ins$d^moment * kdeferredqx(k, mort)
  }
  
  total + kpx(ins$n, mort) * pv.moment(ins$n, moment, irm) * ins$e^moment
}

z.ev.two.isingle.endowsingle <- function(ins, mort, irm) 
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      mu <- -(y.ev(i+1, irm) + y.ev(j+1, irm))
      sigma2 <- y.var(i+1,irm) + y.var(j+1, irm) + 2*y.cov(i+1,j+1,irm)
      EPV <- exp(mu + 0.5* sigma2) # EPV not accounting for mortality
      
      total = total + EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) * ins$d^2
    }
    
    mu <- -(y.ev(i+1, irm) + y.ev(ins$n, irm))
    sigma2 <- y.var(i+1, irm) + y.var(ins$n, irm) + 2*y.cov(i+1, ins$n, irm)
    EPV <- exp(mu + 0.5*sigma2)
    
    total = total + 2 * EPV * kdeferredqx(i, mort) * kpx(ins$n, mort) * ins$d * ins$e
  }
  
  total + pv.moment(ins$n, 2, irm) * kpx(ins$n, mort)^2 * ins$e^2
}

z.ev.three.isingle.endowsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      for(k in seq(0, ins$n-1, 1))
      {
        mu = -(y.ev(i+1, irm) + y.ev(j+1, irm) + y.ev(k+1, irm))
        sigma2 = y.var(i+1, irm) + y.var(j+1, irm) + y.var(k+1, irm) +
          2 * y.cov(i+1, j+1, irm) + 2 * y.cov(i+1, k+1, irm) + 2 * y.cov(j+1, k+1, irm)
        
        EPV <- exp(mu + 0.5* sigma2) # EPV not accounting for mortality
        
        total = total + EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) * 
          kdeferredqx(k, mort) * ins$d^3
      }
      
      mu <- -(y.ev(i+1, irm) + y.ev(j+1, irm) + y.ev(ins$n, irm))
      sigma2 <- y.var(i+1, irm) + y.var(j+1, irm) + y.var(ins$n, irm) +
        2 * y.cov(i+1, j+1, irm) + 2 * y.cov(i+1, ins$n, irm) + 2 * y.cov(j+1, ins$n, irm)
      EPV <- exp(mu + 0.5* sigma2)
      
      total = total + 3 * EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) *
        kpx(ins$n, mort) * ins$d^2 * ins$e
    }
    
    mu <- -(y.ev(i+1, irm) + 2*y.ev(ins$n, irm))
    sigma2 <- y.var(i+1, irm) + 4*y.var(ins$n, irm) + 4*y.cov(i+1, ins$n, irm)
    EPV <- exp(mu + 0.5* sigma2)
    
    total = total + 3 * EPV * kdeferredqx(i, mort) * kpx(ins$n, mort)^2 * ins$d * ins$e^2
  }
  
  total + pv.moment(ins$n, 3, irm) * kpx(ins$n, mort)^3 * ins$e^3
}

z.ev.twoone.isingle.endowsingle <- function(ins, mort, irm)
{
  total = 0
  
  for(i in seq(0, ins$n-1, 1))
  {
    for(j in seq(0, ins$n-1, 1))
    {
      mu <- -(2*y.ev(i+1, irm) + y.ev(j+1, irm))
      sigma2 <- 4*y.var(i+1, irm) + y.var(j+1, irm) + 4*y.cov(i+1,j+1, irm)
      EPV <- exp(mu + 0.5* sigma2) # EPV not accounting for mortality
      
      total = total + EPV * kdeferredqx(i, mort) * kdeferredqx(j, mort) * ins$d^2 * ins$d
    }
    
    mu <- -(2*y.ev(i+1, irm) + y.ev(ins$n, irm))
    sigma2 <- 4*y.var(i+1, irm) + y.var(ins$n, irm) + 4*y.cov(i+1, ins$n, irm)
    EPV <- exp(mu + 0.5 * sigma2)
    
    total = total + EPV * kdeferredqx(i, mort) * kpx(ins$n, mort) * ins$e * ins$d^2
    
    mu <- -(y.ev(i+1, irm) + 2*y.ev(ins$n, irm))
    sigma2 <- y.var(i+1, irm) + 4*y.var(ins$n, irm) + 4*y.cov(i+1, ins$n, irm)
    EPV <- exp(mu + 0.5 * sigma2)
    
    total = total + EPV * kdeferredqx(i, mort) * kpx(ins$n, mort) * ins$e^2 * ins$d
  }
  
  total + pv.moment(ins$n, 3, irm) * kpx(ins$n, mort)^2 * ins$e^3
}

z.insrisk.isingle.endowsingle <- function(ins, mort, irm)
{
  stop("z.insrisk is not implemented for endowsingle class")
}

z.invrisk.isingle.endowsingle <- function(ins, mort, irm)
{
  stop("z.invrisk is not implemented for endowsingle class")
}

z.pdf.isingle.endowsingle <- function(z, ins, mort, irm)
{
  f <- function(z)
  {
    if(z <= 0)
    {
      return(0)
    }
    else
    {
      if(ins$n > 0)
      {
        k = seq(0, ins$n - 1, 1)
        kqx = kdeferredqx(k, mort)
        
        sum(kqx * dlnorm(z, -y.ev(k + 1, irm), sqrt(y.var(k + 1, irm)))) + 
        kpx(ins$n, mort) * dlnorm(z, -y.ev(ins$n, irm), sqrt(y.var(ins$n, irm)))
      }
      else
      {
        ifelse(n == 0 && z == 1, 1, 0)
      }
    }
  }
  
  sapply(z, f)
}

## Insurance Portfolio (identical policies)

z.moment.iport <- function(moment, ins, mort, irm)
{
  UseMethod("z.moment.iport", ins)
}

z.insrisk.iport <- function(ins, mort, irm)
{
  stop("z.insrisk not implemented for iport classes")
}

z.invrisk.iport <- function(ins, mort, irm)
{
  stop("z.invrisk not implemented for iport classes")
}

z.pdf.iport <- function(z, ins, mort, irm)
{
  stop("z.pdf not implemented for iport classes")
}

## Term Insurance Portfolio

z.moment.iport.termport <- function(moment, ins, mort, irm)
{
  epv = 0
  
  if(moment == 1)
  {
    epv = z.moment(1, ins$single, mort, irm) * ins$c
  }
  else if(moment == 2)
  {
    epv = ins$c * (ins$c - 1) * z.ev.two.isingle.termsingle(ins$single, mort, irm) + 
      ins$c * z.moment(2, ins$single, mort, irm)
  }
  else if(moment == 3)
  {
    epv = ins$c * (ins$c - 1) * (ins$c - 2) * z.ev.three.isingle.termsingle(ins$single, mort, irm) +
      3 * ins$c * (ins$c - 1) * z.ev.twoone.isingle.termsingle(ins$single, mort, irm) +
      ins$c * z.moment(3, ins$single, mort, irm)
  }
  else
  {
    stop("moment > 3 not implemented for termport")
  }
  
  epv
}

## Endowment Insurance Portfolio

z.moment.iport.endowport <- function(moment, ins, mort, irm)
{
  epv = 0
  
  if(moment == 1)
  {
    epv = z.moment(1, ins$single, mort, irm) * ins$c
  }
  else if(moment == 2)
  {
    epv = ins$c * (ins$c - 1) * z.ev.two.isingle.endowsingle(ins$single, mort, irm) + 
      ins$c * z.moment(2, ins$single, mort, irm)
  }
  else if(moment == 3)
  {
    epv = ins$c * (ins$c - 1) * (ins$c - 2) * z.ev.three.isingle.endowsingle(ins$single, mort, irm) +
      3 * ins$c * (ins$c - 1) * z.ev.twoone.isingle.endowsingle(ins$single, mort, irm) +
      ins$c * z.moment(3, ins$single, mort, irm)
  }
  else
  {
    stop("moment > 3 not implemented for endowport")
  }
  
  epv
}

## Group Insurance Portfolio

z.ev.two.igroup <- function(ind1, ind2, ins, mort, irm)
{
  if(length(ins) != length(mort) + 2)
  {
    stop("length(ins) != length(mort)")
  }
  
  total = 0
  
  for(i in seq(0, ins[[ind1]]$single$n-1, 1))
  {
    for(j in seq(0, ins[[ind2]]$single$n-1, 1))
    {
      mu = -(y.ev(i+1, irm) + y.ev(j+1, irm))
      sigma2 = y.var(i+1, irm) + y.var(j+1, irm) + 2*y.cov(i+1,j+1, irm)
      EPV = exp(mu + 0.5 * sigma2) # EPV not accounting for mortality
      
      total = total + EPV * kdeferredqx(i, mort[[ind1]]) *
        kdeferredqx(j, mort[[ind2]]) * ins[[ind1]]$single$d * ins[[ind2]]$single$d
    }
  }
  
  for(i in seq(0, ins[[ind1]]$single$n-1, 1))
  {
    mu <- -(y.ev(i+1, irm) + y.ev(ins[[ind2]]$single$n, irm))
    sigma2 <- y.var(i+1, irm) + y.var(ins[[ind2]]$single$n, irm) + 
      2*y.cov(i+1, ins[[ind2]]$single$n, irm)
    EPV <- exp(mu + 0.5*sigma2)
    
    total = total + EPV * kdeferredqx(i, mort[[ind1]]) *
      kpx(ins[[ind2]]$single$n, mort[[ind2]]) * ins[[ind1]]$single$d * ins[[ind2]]$single$e
  }
  
  for(i in seq(0, ins[[ind2]]$single$n-1, 1))
  {
    mu <- -(y.ev(i+1, irm) + y.ev(ins[[ind1]]$single$n, irm))
    sigma2 <- y.var(i+1, irm) + y.var(ins[[ind1]]$single$n, irm) + 
      2*y.cov(i+1, ins[[ind1]]$single$n, irm)
    EPV <- exp(mu + 0.5*sigma2)
    
    total = total + EPV * kpx(ins[[ind1]]$single$n, mort[[ind1]]) *
      kdeferredqx(i, mort[[ind2]]) * ins[[ind2]]$single$d * ins[[ind1]]$single$e
  }
  
  mu = -(y.ev(ins[[ind1]]$single$n, irm) + y.ev(ins[[ind2]]$single$n, irm))
  sigma2 = y.var(ins[[ind1]]$single$n, irm) + y.var(ins[[ind2]]$single$n, irm) + 
    2*y.cov(ins[[ind1]]$single$n, ins[[ind2]]$single$n, irm)
  EPV = exp(mu + 0.5*sigma2)
  
  total + EPV * kpx(ins[[ind1]]$single$n, mort[[ind1]]) *
    kpx(ins[[ind2]]$single$n, mort[[ind2]]) * ins[[ind1]]$single$e * ins[[ind2]]$single$e
}

z.moment.igroup <- function(moment, ins, mort, irm)
{
  if(length(ins) != length(mort) + 2)
  {
    stop("length(ins) != length(mort)")
  }
  
  m = length(mort) # number of groups
  
  epv = 0
  
  if(moment == 1)
  {
    for(i in 1:m)
    {
      epv = epv + z.moment(1, ins[[i]], mort[[i]], irm)
    }
  }
  else if(moment == 2)
  {
    for(i in 1:m)
    {
      epv = epv + ins[[i]]$c * z.moment(2, ins[[i]]$single, mort[[i]], irm)
    }
    
    for(i in 1:m)
    {
      epv = epv + ins[[i]]$c * (ins[[i]]$c - 1) * 
        z.ev.two.isingle.endowsingle(ins[[i]]$single, mort[[i]], irm)
    }
    
    for(i in 1:m)
    {
      for(r in i:m)
      {
        if(i != r)
        {
          epv = epv + 2 * ins[[i]]$c * ins[[r]]$c *
            z.ev.two.igroup(i, r, ins, mort, irm)
        }
      }
    }
  }
  
  epv
}

cashflow.ev <- function(r, ins, mort, irm)
{
  if(length(ins) != length(mort) + 2)
  {
    stop("length(ins) != length(mort)")
  }
  
  m = length(mort) # number of groups
  
  epv = 0
  
  I1 <- function(n, r)
  {
    ifelse(n >= r, 1, 0)
  }
  
  I2 <- function(n, r)
  {
    ifelse(n == r, 1, 0)
  }
  
  for(i in 1:m)
  {
    epv = epv + ins[[i]]$single$d * ins[[i]]$c * kdeferredqx(r - 1, mort[[i]]) *
      I1(ins[[i]]$single$n, r)
  }
  
  for(i in 1:m)
  {
    epv = epv + ins[[i]]$single$e * ins[[i]]$c * kpx(ins[[i]]$single$n, mort[[i]]) *
      I2(ins[[i]]$single$n, r)
  }
  
  epv
}

cashflow.cov <- function(s, r, ins, mort, irm)
{
  if(length(ins) != length(mort) + 2)
  {
    stop("length(ins) != length(mort)")
  }
  
  m = length(mort) # number of groups
  
  total = 0
  
  I1 <- function(n, r)
  {
    ifelse(n >= r, 1, 0)
  }
  
  I2 <- function(n, r)
  {
    ifelse(n == r, 1, 0)
  }
  
  if(s == r) # variance
  {
    for(i in 1:m)
    {
      total = total + ins[[i]]$single$d^2 * ins[[i]]$c *
        kdeferredqx(r - 1, mort[[i]]) *
        (1 - kdeferredqx(r - 1, mort[[i]])) * I1(ins[[i]]$single$n, r)
    }
    
    for(i in 1:m)
    {
      total = total + ins[[i]]$single$e^2 * ins[[i]]$c *
        kpx(ins[[i]]$single$n, mort[[i]]) *
        (1 - kpx(ins[[i]]$single$n, mort[[i]])) * I2(ins[[i]]$single$n, r)
    }
    
    for(i in 1:m)
    {
      total = total - 2 * ins[[i]]$single$d * ins[[i]]$single$e * ins[[i]]$c *
        kdeferredqx(r - 1, mort[[i]]) *
        kpx(r, mort[[i]]) * I2(ins[[i]]$single$n, r)
    }
  }
  else
  {
    for(i in 1:m)
    {
      total = total - ins[[i]]$single$d^2 * ins[[i]]$c *
        kdeferredqx(min(s,r) - 1, mort[[i]]) *
        kdeferredqx(max(s,r) - 1, mort[[i]]) * I1(ins[[i]]$single$n, max(s,r))
    }
    
    for(i in 1:m)
    {
      total = total - ins[[i]]$single$d * ins[[i]]$single$e * ins[[i]]$c *
        kdeferredqx(min(s,r) - 1, mort[[i]]) *
        kpx(ins[[i]]$single$n, mort[[i]]) * I2(ins[[i]]$single$n, max(s,r))
    }
  }
  
  total
}

z.insrisk.igroup <- function(ins, mort, irm)
{
  total = 0
  
  for(r in 1:ins$n)
  {
    for(s in 1:ins$n)
    {
      mu <- -(y.ev(r, irm) + y.ev(s, irm))
      sigma2 <- y.var(r, irm) + y.var(s, irm) + 2*y.cov(r, s, irm)
      EPV <- exp(mu + 0.5*sigma2)
      
      total = total +
        EPV * cashflow.cov(r, s, ins, mort, irm)
    }
  }
  
  total
}

z.invrisk.igroup <- function(ins, mort, irm)
{
  total = 0
  
  for(r in 1:ins$n)
  {
    for(s in 1:ins$n)
    {
      total = total + cashflow.ev(r, ins, mort, irm) *
        cashflow.ev(s, ins, mort, irm) * pv.cov(r, s, irm)
    }
  }
  
  total
}