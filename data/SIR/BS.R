library(smcUtils)

load("sims.RData")
source("filter-settings.R")

n = max(d$time)+1
weights = array(NA, dim=c(n.particles, n, n.sims))
X = array(NA, dim=c(sckm$s, n.particles, n, n.sims))
p = r = array(NA, dim=c(sckm$r, n.particles, n))

for (i in 1:n.sims)
{
  # Draw prior
  weights[,1,i] = 1/n.particles
  X[,,1,i] = sckm$X
  for (j in 1:sckm$r)
  {
    p[j,,1,i] = rbeta(n.particles, prior$prob$a[j], prior$prob$b[j])
    r[j,,1,i] = rbeta(n.particles, prior$prob$a[j], prior$prob$b[j])
  }

  for (j in (d$time+1))
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
    }

    # Weights
    for (k in 1:n.particles)
    {
      weights[k,j,i] = sum(dbinom(data[j,7:8], tmp$nr, p[,kk[k],j,i], log=T))
    }
    weights[,j,i] = renormalize.weights(weights[,j,i], log=T)

  }
}

rm(n, tmp, kk)
save(weights,X,p,r, file="BS.RData")

#q("no")

