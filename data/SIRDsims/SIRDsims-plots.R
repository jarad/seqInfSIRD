# SIRD simulation plots

load("SIRDsims-mcmc-density.RData")

par(mfrow=c(2,2))
for (i in 1:length(kd)) {
  plot(kd[[i]], main=key[i])
}

