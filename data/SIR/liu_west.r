library(tlpl)
library(smcUtils)

logit = function(p) log(p/(1-p))
expit = function(r) 1/(1+1/exp(r))

liu_west = function(y, sckm, n.particles, delta, prior,...)
{
  # Liu-West parameter
  a = (3*delta-1)/(2*delta)
  h = sqrt(1-a^2)

  # Indices
  p.i = 1:sckm$r
  r.i = sckm$r + (1:sckm$r)

  # Set up data structures
  n = nrow(y)
  weights = lweights = matrix(NA, n.particles, n+1)
  X       = array(NA, dim=c(sckm$s, n.particles, n+1))
  theta   = array(NA, dim=c(2*sckm$r, n.particles, n+1))

  ## Initialize
  X[,,1] = sckm$X
  for (i in p.i) 
  {
    theta[i,,1]        = logit(rbeta( n.particles, prior$prob$a[i], prior$prob$b[i]))
    theta[i+sckm$r,,1] = log(  rgamma(n.particles, prior$rate$a[i], prior$rate$b[i]))
  }
  weights[,1] = 1/n.particles
  lweights[,1] = log(weights[,1])
 
  # Variables used in algorithm
  w      = numeric(n.particles)
  exp.nr = matrix(NA,sckm$r,n.particles)
  m      = matrix(NA, 2*sckm$r, n.particles)


  for (i in 1:n)
  {
    # Calculate sample mean and covariance of weighted particles
    #mn = apply(theta[,,i], 1, wtd.mean, normwt=T, weights=weights[,i])
    #cv = cov(t(theta[,,i]))

    cw = cov.wt(t(theta[,,i]), weights[,i])
    mn = cw$center
    cv = cw$cov
    
    ch = tryCatch(t(chol(cv, pivot=FALSE)),
                  warning = function(w) { print(w) },
                  error   = function(e) { print(e); NULL })
    # If cholesky fails, use APF but do not jitter particles


    # Shrink each particle toward grand mean
    if (is.null(ch)) # Cholesky failed
    {
      m = theta[,,i]
    } else 
    {                # Cholesky succeeded
      for (j in 1:n.particles)
      {
        m[,j] = a*theta[,j,i] + (1-a)*mn
      }
    }

    # Calculate APF weights
    for (j in 1:n.particles) 
    {
      sckm$X = X[,j,i]
      sckm$theta = exp(theta[r.i,j,i])
      exp.nr[,j] = ceiling(hazard(sckm)$h) # Expected number of reactions/transitions
      w[j] = lweights[j,i] + sum(dbinom(y[i,], exp.nr[,j], expit(m[p.i,j]), log=T))
    }

    # Resample particles
    w = renormalize(w, log=T)

    # If all particles have zero probability 
    if (all(is.nan(w))) 
    {
        warning(paste("All particles have zero probability at time ",i,". Aborting run.",sep=""))
        return(list(weights=weights, X=X, 
                    p = expit(theta[p.i,,]), 
                    r = exp(  theta[r.i,,])))
    }

    
    res = resample(w, log=F, normalized=T,...)
    kk = res$indices

    # Propagate
    for (j in 1:n.particles)
    {
      k = kk[j]

      if (is.null(ch))
      {
        theta[,j,i+1] = m[,k] # = theta[,k,i]
      } else 
      {
        theta[,j,i+1] = m[,k] + h*ch%*%rnorm(2*sckm$r)
      }

      sckm$theta = exp(theta[r.i,j,i+1])
      sckm$X = X[,k,i]
      tmp = tau_leap(sckm, engine="C")
      X[,j,i+1] = tmp$X[2,]
      lweights[j,i+1] = log(res$weights[k])+ # In case particles were not resampled
                        sum(dbinom(y[i,], tmp$nr,     expit(theta[p.i,j,i+1]), log=T)-
                            dbinom(y[i,], exp.nr[,k], expit(    m[p.i,k    ]), log=T))
    }
    weights[,i+1] = renormalize(lweights[,i+1], log=T)
  }

  return(list(weights=weights, X=X, 
              p = expit(theta[p.i,,]), 
              r = exp(  theta[r.i,,])))
}


