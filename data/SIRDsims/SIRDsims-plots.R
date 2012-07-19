# SIRD simulation plots

load("SIRDsims-mcmc-density.RData")
load("SIRDsims.RData")
load("SIRDsims-smcLW-density.RData") 
load("SIRDsims-smcPL-density.RData") 
load("SIRDsims-smcSV-density.RData") 

truth = NULL
for (i in 1:length(n)) {
  truth = c(truth,sims[[i]]$x$S[n[i]], sims[[i]]$x$I[n[i]], thetas[i,1:2])
}

par(mfrow=c(2,2))
for (i in 1:length(kd)) {
  plot(kd[[i]], main=key[i], 
       xlim=range(kd[[i]]$x, kd.LW[[i]]$x, kd.SV[[i]]$x, kd.PL[[i]]$x),
       ylim=range(kd[[i]]$y, kd.LW[[i]]$y, kd.SV[[i]]$y, kd.PL[[i]]$y))
  lines(kd.LW[[i]], col=2)
  lines(kd.SV[[i]], col=3)
  lines(kd.PL[[i]], col=4)
  abline(v=truth[i], col="black")
  legend("bottomright", c("MCMC","LW","SV","PL"), col=1:4, lwd=1, bg="white")
  if (i%%4==0) readline("hit <enter>:")
}

