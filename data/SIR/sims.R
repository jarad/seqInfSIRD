library(tlpl)
library(plyr)

# Generate data
## Set up SIR model
sckm = list()
sckm$s = 3 # species (S,I,R)
sckm$r = 2 # reactions (S->I, I->R)
#                   S -> I    I -> R
sckm$Pre  = rbind( c(1,1,0), c(0,1,0))
sckm$Post = rbind( c(0,2,0), c(0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(16000,10,0)
N = sum(sckm$X)
sckm$lmult = log(c(1/N,1))

## Simulate data
n.sims = 10
set.seed(20121218)
n = 50

prior = list(prob=list(a=rep(1,sckm$r), b=rep(1,sckm$r)),
             rate=list(a=c(1,.5)*10, b=rep(10,sckm$r)))

# A function to produce a single simulation
sim.f = function()
{
  try({
    rates = rgamma(sckm$r, prior$rate$a, prior$rate$b)
    sckm$theta = rates
    out = tau_leap(sckm, n)
    out$rates = rates

    out$probs = rbeta( sckm$r, prior$prob$a, prior$prob$b)
    out$y = t(rbind(rbinom(n, out$nr[,1], out$p[1]), 
                    rbinom(n, out$nr[,2], out$p[2])))
    }, silent=T)

  if (exists("out")) return(out) else { return(NA) }
}

# Apply the function n.sims times and return as a list
sims = rlply(n.sims, sim.f, .progress="time")

# If a simulation had an error, resimulate
for (i in 1:n.sims)
{
  while(!is.list(sims[[i]])) sims[[i]] = sim.f()
}

save.image("sims.RData")

