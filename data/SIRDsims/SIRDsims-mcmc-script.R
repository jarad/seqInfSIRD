require(rjags)
load("SIRDsims.RData")

i = 1
cols = 1:4
dat   = list(X0=X0[cols], y=as.matrix(sims[[i]]$y[,cols]), p=probs[i,cols], 
             n=nrow(sims[[i]]$y), 
             N=sum(X0[cols]), a=prior$theta$a[cols], b=prior$theta$b[cols])
inits = list(dx=as.matrix(sims[[i]]$dx[,cols]),
             theta=thetas[i,cols])

ptr = proc.time()

mod   = jags.model("../../code/working/SIRD.txt", data=dat, 
                  inits=inits, n.adapt=1e6)
res   = coda.samples(mod, c("theta","x"), 1e6, thin=100) 

run.time = proc.time()-ptr

res.quantiles = apply(as.matrix(res), 2, function(x) quantile(x,c(.025,.975)))

write.csv(res.quantiles,"SIRDsims-mcmc-quantiles.csv",row.names=T)

save.image("SIRDsims-mcmc-script.RData")

q("no")

