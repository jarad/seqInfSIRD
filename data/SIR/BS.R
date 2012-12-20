library(smcUtils)
library(tlpl)

load("sims.RData")
source("filter-settings.R")

weights = array(NA, dim=c(n.particles, n+1, n.sims))
X = array(NA, dim=c(sckm$s, n.particles, n+1, n.sims))
p = r = array(NA, dim=c(sckm$r, n.particles, n+1, n.sims))

for (i in 1:n.sims)
{
  d = data[[i]]

  # Draw prior
  weights[,1,i] = 1/n.particles
  X[,,1,i] = sckm$X
  for (j in 1:sckm$r)
  {
    p[j,,1,i] = rbeta( n.particles, prior$prob$a[j], prior$prob$b[j])
    r[j,,1,i] = rgamma(n.particles, prior$rate$a[j], prior$rate$b[j])
  }

  for (j in 1:n )
  {
    # Resample
    if (j>1) 
    {
      kk = resample(weights[,j,i], method=resampling.function)$indices
    } else {
      kk = 1:n.particles
    }

    # Propagate
    for (k in 1:n.particles)
    {
      sckm$X     = X[,kk[k],j,i]
      sckm$theta = r[,kk[k],j,i]
      tmp = tau_leap(sckm)
      X[,k,j+1,i] = tmp$X[2,]
      r[,k,j+1,i] = r[,kk[k],j,i]
      p[,k,j+1,i] = p[,kk[k],j,i]
    }

    # Weights
    for (k in 1:n.particles)
    {
      weights[k,j+1,i] = sum(dbinom(d[j,], tmp$nr, p[,k,j+1,i], log=T))
    }
    weights[,j+1,i] = renormalize.weights(weights[,j+1,i], log=T)
  }
}

rm(n, tmp, kk)
save(weights,X,p,r, file="BS.RData")

#q("no")

