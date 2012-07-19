# SIRD simulation plots

load("SIRDsims-mcmc-density.RData")
load("SIRDsims.RData")
load("SIRDsims-smcLW-density.RData") 
load("SIRDsims-smcPL-density.RData") 
load("SIRDsims-smcSV-density.RData") 

mcmc.quantile <- read.csv("../../data/SIRDsims/SIRDsims-mcmc-quantiles.csv")
pl.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcPL-quantiles.csv")
lw.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcLW-quantiles.csv")
sv.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcSV-quantiles.csv")


truth = NULL
for (i in 1:length(n)) {
  truth = c(truth,sims[[i]]$x$S[n[i]], sims[[i]]$x$I[n[i]], thetas[i,1:2])
}

par(mfrow=c(2,2),mar=c(4,4,2,1), oma=c(0,0,0,1))
# Make 4 panels, each panel shows the 4 different algorithm outputs (in terms of terminal densities)

stat.ndx <- list();
stat.ndx[[1]] <- 1:3;  stat.ndx[[2]] <- 7:9; 
stat.ndx[[3]] <- 4:6;  stat.ndx[[4]] <- 10:12;
for (i in 1:length(kd)) {
  k <- (i-1) %% 4 + 1
  j <- floor((i-1)/4)+1
  plot(kd[[i]], main=key[i], 
       xlim=range(kd[[i]]$x, kd.LW[[i]]$x, kd.SV[[i]]$x, kd.PL[[i]]$x),
       ylim=range(kd[[i]]$y, kd.LW[[i]]$y, kd.SV[[i]]$y, kd.PL[[i]]$y))
  lines(kd.LW[[i]], col=2)
  lines(kd.SV[[i]], col=3)
  lines(kd.PL[[i]], col=4)
  points(mcmc.quantile[j,stat.ndx[[k]] ],rep(0,3),pch=4,col=1, cex=1.5)
  points(lw.quantile[j,stat.ndx[[k]] ],rep(0,3),pch=5,col=2, cex=1.5)
  points(sv.quantile[j,stat.ndx[[k]] ],rep(0,3),pch=6,col=3, cex=1.5)
  points(pl.quantile[j,stat.ndx[[k]] ],rep(0,3),pch=7,col=4, cex=1.5)
  
  abline(v=truth[i], col="black")
  if (k==2)
    legend("right", c("MCMC","LW","SV","PL"), col=1:4, lwd=1, bg="white")
  if (i%%4==0) readline("hit <enter>:")
}

