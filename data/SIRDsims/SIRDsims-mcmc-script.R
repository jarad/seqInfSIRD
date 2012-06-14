require(rjags)
load("SIRDsims.RData")

i = 1
cols = 1:4
dat   = list(X0=X0[cols], y=as.matrix(sims[[i]]$y[,cols]), p=probs[i,cols], 
             n=nrow(sims[[i]]$y), 
             N=sum(X0[cols]), a=prior$theta$a[cols], b=prior$gamma$b[cols])
inits = list(dx=as.matrix(sims[[i]]$dx[,cols]),
             theta=gammas[i,cols])

mod   = jags.model("../../code/working/SIR.txt", data=dat, 
                  inits=inits, n.adapt=1e3)
res   = coda.samples(mod, c("theta","x"), 1e3, thin=10) 


ptr = proc.time()

mod   = jags.model("../../code/working/SIRD.txt", data=dat, 
                  inits=inits, n.adapt=1e5)
res   = coda.samples(mod, c("gamma","x"), 1e5, thin=10) 

run.time = proc.time()-ptr

save.image("SIRDsims-mcmc-script.RData")

q("no")

