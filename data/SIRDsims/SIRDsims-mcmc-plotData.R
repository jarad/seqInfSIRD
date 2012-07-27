load("SIRDsims-mcmc-script.RData")

# Generate kernel density estimates for 
# last time point of S, I, S->I, and I->R
n.sims = length(res)

kd = list()
key = paste("simulation ",rep(1:n.sims,each=4),", ",c("S","I","S->I","I->R"),sep="")
ii = 1
for (i in 1:n.sims) {
  kd[[ii]] = density(unlist(res[[i]][,2+n[i]])) # S
  ii = ii+1
  dat = unlist(res[[i]][,2+2*n[i]])
  x = seq(range(dat)[1],range(dat)[2])
  y = rep(NA,length(x))
  for (j in 1:length(x)) y[j] = mean(dat==x[j])
  kd[[ii]] = list(x=x,y=y) # I
  ii = ii+1
  kd[[ii]] = density(unlist(res[[i]][,1])) # S->I
  ii = ii+1
  kd[[ii]] = density(unlist(res[[i]][,2])) # I->R
  ii = ii+1
}

save(kd,key,file="SIRDsims-mcmc-density.RData")

# Save median and 95% credible intervals

sims.MCMC = matrix(NA,n.sims,3*4)
for (i in 1:n.sims) {
  cols = c(2+n[i],1,2+2*n[i],2)
  tmp = apply(as.matrix(res[[i]][,cols]), 2, function(x) quantile(x,prob=c(.5,.025,.975)))
  sims.MCMC[i,] = as.vector(tmp[])
}

colnames(sims.MCMC) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(sims.MCMC,"SIRDsims-mcmc-quantiles.csv",row.names=F)


################################################
# Sequential MCMC analysis
################################################
rm(list=ls())
load("SIRDsims-mcmc-sequential.RData")

sims.MCMC = matrix(NA,n,3*4)
for (i in 2:n) {
  cols = c(2+i,1,2+2*i,2)
  tmp = apply(as.matrix(res[[i]][,cols]), 2, function(x) quantile(x,prob=c(.5,.025,.975)))
  sims.MCMC[i,] = as.vector(tmp[])
}

colnames(sims.MCMC) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(sims.MCMC,"SIRDsims-mcmc-seq-quants.csv",row.names=F)

