require(rjags)
load("SIRDsims.RData")

# determines which simulation should be used for sequential analysis
sim.id = which.max(n)
sim = sims[[sim.id]]
n = n[sim.id]

N.RXNS = 2
N.STATES = 3

res = list()
run.time = matrix(NA,n,3)
for (i in 2:n) {
  cat("Time point",i,"of",n, "\n")
  X0 = as.numeric(sim$x[1,1:N.STATES])
  dat   = list(X0 = X0, 
               y  = as.matrix(sim$y[1:i,1:N.RXNS]), 
               p=probs[sim.id,1:N.RXNS], 
               n=i, 
               N=sum(X0), a=prior$theta$a[1:N.RXNS], b=prior$theta$b[1:N.RXNS])
  inits = list(dx=as.matrix(sim$dx[1:i,1:N.RXNS]),
               theta=thetas[sim.id,1:N.RXNS])

  ptr = proc.time()[1:3]

  mod   = jags.model("../../code/working/SIR.txt", data=dat, 
                    inits=inits, n.adapt=1e5)
  res[[i]]   = coda.samples(mod, c("theta","x"), 1e6, thin=10) 

  run.time[i,] = proc.time()[1:3]-ptr

  save.image("SIRDsims-mcmc-sequential.RData")
}

#res.quantiles = apply(as.matrix(res), 2, function(x) quantile(x,c(.025,.975)))

#write.csv(res.quantiles,"SIRDsims-mcmc-quantiles.csv",row.names=T)

save.image("SIRDsims-mcmc-sequential.RData")

q("no")

