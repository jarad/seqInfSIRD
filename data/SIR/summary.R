# Read in truth, means, and quantiles and 
# calculate MSE, MAD, MAPE, coverage
#
load("sims.RData")
rm(data)

source("quantile-settings.R")

load("PLq.RData")
load("LW90q.RData"); lw90q = lwq
load("LW95q.RData"); lw95q = lwq
load("LW99q.RData"); lw99q = lwq

# Build data frames
i.m = which(probs==0.5)

states = expand.grid( time=0:n, sim=1:n.sims, method=c("PL","LW90","LW95","LW99"))

trueS = trueI = trueR = NULL
for (i in 1:n.sims)
{
  trueS = c(trueS, truth[[i]]$S)
  trueI = c(trueI, truth[[i]]$I)
  trueR = c(trueR, truth[[i]]$R)
}
states$trueS = rep(trueS, nlevels(states$method))
states$trueI = rep(trueI, nlevels(states$method))
states$trueR = rep(trueR, nlevels(states$method))


