require(rjags)
load("SIRDsims.RData")

N.RXNS = 2
N.STATES = 3

i = 1
  X0 = as.numeric(sims[[i]]$x[1,1:N.STATES])
  dat   = list(X0 = X0, 
               y  = as.matrix(sims[[i]]$y[1:n[i],1:N.RXNS]), 
               p=probs[i,1:N.RXNS], 
               n=n[i], 
               N=sum(X0), a=prior$theta$a[1:N.RXNS], b=prior$theta$b[1:N.RXNS])
  inits = list(dx=as.matrix(sims[[i]]$dx[1:n[i],1:N.RXNS]),
               theta=thetas[i,1:N.RXNS])

ptr = proc.time()

mod   = jags.model("../../code/working/SIR.txt", data=dat, 
                  inits=inits, n.adapt=1e6)
res   = coda.samples(mod, c("theta","x"), 1e6, thin=10) 

run.time = proc.time()-ptr

res.quantiles = apply(as.matrix(res), 2, function(x) quantile(x,c(.025,.975)))

write.csv(res.quantiles,"SIRDsims-mcmc-quantiles.csv",row.names=T)

save.image("SIRDsims-mcmc-script.RData")

q("no")

