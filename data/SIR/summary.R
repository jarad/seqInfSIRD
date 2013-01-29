# Read in truth, means, and quantiles and 
# calculate MSE, MAD, MAPE, coverage
#
load("sims.RData")

source("quantile-settings.R")
N.METHODS = 4
load("PLq.RData")
load("LW90q.RData"); lw90q = lwq
load("LW95q.RData"); lw95q = lwq
load("LW99q.RData"); lw99q = lwq
method = c("PL","LW90","LW95","LW99")

states = sckm$states
rxns   = sckm$rxns
ns     = sckm$s
nr     = sckm$r
comp = c(states, paste("p:",rxns), paste("r:",rxns))


nn = c(length(plq), ns+2*nr, n+1, N.METHODS)

# MSE
which.md = which(probs==0.5)
SE = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    SE[i,j,,1] = (sims[[i]]$X[,j] - plq  [[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,2] = (sims[[i]]$X[,j] - lw90q[[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,3] = (sims[[i]]$X[,j] - lw95q[[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,4] = (sims[[i]]$X[,j] - lw99q[[i]]$X.quantiles[j,which.md,])^2
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + j
    SE[i,k,,1] = (sims[[i]]$probs[j] - plq  [[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,2] = (sims[[i]]$probs[j] - lw90q[[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,3] = (sims[[i]]$probs[j] - lw95q[[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,4] = (sims[[i]]$probs[j] - lw99q[[i]]$p.quantiles[j,which.md,])^2
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + sckm$r + j
    SE[i,k,,1] = (sims[[i]]$rates[j] - plq  [[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,2] = (sims[[i]]$rates[j] - lw90q[[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,3] = (sims[[i]]$rates[j] - lw95q[[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,4] = (sims[[i]]$rates[j] - lw99q[[i]]$r.quantiles[j,which.md,])^2
  }
}

MSE = apply(SE, 2:4, mean, na.rm=T)

for (i in 1:nn[2])
{
  plot(0:n, MSE[i,,1], type="l", ylim=range(MSE[i,,], na.rm=T), main=comp[i])
  for (j in 2:4) lines(0:n, MSE[i,,j], col=j, lty=j)
  legend("topright",method,col=1:4, lty=1:4)
  readline("")
}


# Coverage
which.ci = c(1,5)
CV = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    for (l in 1:(n+1))
    {
      CV[i,j,l,1] = sum(sims[[i]]$X[l,j] > plq  [[i]]$X.quantiles[j,which.ci,l])==1
      CV[i,j,l,2] = sum(sims[[i]]$X[l,j] > lw90q[[i]]$X.quantiles[j,which.ci,l])==1
      CV[i,j,l,3] = sum(sims[[i]]$X[l,j] > lw95q[[i]]$X.quantiles[j,which.ci,l])==1
      CV[i,j,l,4] = sum(sims[[i]]$X[l,j] > lw99q[[i]]$X.quantiles[j,which.ci,l])==1
    }
  }
  for (j in 1:sckm$r)
  {
    for (l in 1:(n+1))
    {
      k = sckm$s + j
      CV[i,k,l,1] = sum(sims[[i]]$probs[j] > plq  [[i]]$p.quantiles[j,which.ci,l])==1
      CV[i,k,l,2] = sum(sims[[i]]$probs[j] > lw90q[[i]]$p.quantiles[j,which.ci,l])==1
      CV[i,k,l,3] = sum(sims[[i]]$probs[j] > lw95q[[i]]$p.quantiles[j,which.ci,l])==1
      CV[i,k,l,4] = sum(sims[[i]]$probs[j] > lw99q[[i]]$p.quantiles[j,which.ci,l])==1
    }
  }
  for (j in 1:sckm$r)
  {
    for (l in 1:(n+1))
    {
      k = sckm$s + sckm$r + j
      CV[i,k,l,1] = sum(sims[[i]]$rates[j] > plq  [[i]]$r.quantiles[j,which.ci,l])==1
      CV[i,k,l,2] = sum(sims[[i]]$rates[j] > lw90q[[i]]$r.quantiles[j,which.ci,l])==1
      CV[i,k,l,3] = sum(sims[[i]]$rates[j] > lw95q[[i]]$r.quantiles[j,which.ci,l])==1
      CV[i,k,l,4] = sum(sims[[i]]$rates[j] > lw99q[[i]]$r.quantiles[j,which.ci,l])==1
    }
  }
}

MCV = apply(CV, 2:4, mean, na.rm=T)

for (i in 1:nn[2])
{
  plot(0:n, MCV[i,,1], type="l", ylim=range(MCV[i,,], na.rm=T), main=comp[i])
  for (j in 2:4) lines(0:n, MCV[i,,j], col=j, lty=j)
  legend("topright",method,col=1:4, lty=1:4)
  readline("")
}






