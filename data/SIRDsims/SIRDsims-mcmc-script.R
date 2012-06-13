load("SIRDsims.RData")

i = 1
cols = 1:4
dat   = list(X0=X0[cols], y=as.matrix(sims[[i]]$y[,cols]), p=probs[i,cols], 
             n=nrow(sims[[i]]$y), 
             N=sum(X0[cols]), a=prior$gamma$a[cols], b=prior$gamma$b[cols])
inits = list(dx=as.matrix(sims[[i]]$dx[,cols]),
             gamma=gammas[i,cols])

mod   = jags.model("../../code/working/SIR.txt", data=dat, 
                  inits=inits, n.adapt=1e3)
res   = coda.samples(mod, c("gamma","x"), 1e3, thin=10) 



mod   = jags.model("../../code/working/SIRD.txt", data=dat, inits=inits, n.adapt=1e6)
update(mod,1e5)
res   = coda.samples(mod, c("gamma","x"), 1e5, thin=100) 
