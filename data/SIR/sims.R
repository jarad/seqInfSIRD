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
sckm$X = c(16000,100,0)
N = sum(sckm$X)
sckm$lmult = log(c(1/N,1))
sckm$states = c("S","I","R")
sckm$rxns = c("S->I","I->R")

## Simulate data
n.sims = 30
set.seed(20121218)
n = 60

prior = list(prob=list(a=rep(5,sckm$r), b=rep(95,sckm$r)),
             rate=list(a=c(.5,.25)*10, b=rep(10,sckm$r)))

# A function to produce a single simulation
sim.f = function()
{
  failed = TRUE
  
  try({

    somey = FALSE
    while(!somey) 
    {
      rates = rgamma( sckm$r, prior$rate$a*10, prior$rate$b*10)
      sckm$theta = rates
      out = tau_leap(sckm, n)
      out$rates = rates

      out$probs = rbeta( sckm$r, prior$prob$a*10, prior$prob$b*10)

      out$y = t(rbind(rbinom(n, out$nr[,1], out$p[1]), 
                      rbinom(n, out$nr[,2], out$p[2])))
      somey = sum(out$y[1:5,1]>0)
    }

    }, silent=T)

  if (exists("out")) return(out) else { return(NA) }
}

# Apply the function n.sims times and return as a list
sims = rlply(n.sims, sim.f, 
             .progress = progress_text(style=ifelse(interactive(), 3, 1)))

# If a simulation had an error, resimulate
for (i in 1:n.sims)
{
  while(!is.list(sims[[i]])) sims[[i]] = sim.f()
}

save.image("sims.RData")

q(ifelse(interactive(),"ask","no"))

