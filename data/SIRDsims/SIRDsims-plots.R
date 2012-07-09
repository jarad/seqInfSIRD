# SIRD simulation plots

load("SIRDsims-mcmc-density.RData")
load("SIRDsims.RData")

truth = NULL
for (i in 1:length(n)) {
  truth = c(truth,sims[[i]]$x$S[n[i]], sims[[i]]$x$I[n[i]], thetas[i,1:2])
}

par(mfrow=c(2,2))
for (i in 1:length(kd)) {
  plot(kd[[i]], main=key[i])
  abline(v=truth[i], col="red")
}

