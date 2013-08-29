library(tlpl)
library(plyr)


# Generate data
## Set up SIR model
sckm = sckm("sir", X=c(16000,100,0))

## Simulate data
source("settings.R")
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
      sckm$X = as.numeric(rmultinom(1,N,sckm$X/N))
      sckm$theta = rgamma( sckm$r, prior$rate$a*10, prior$rate$b*10)

      out = tau_leap(sckm, n)
      out$sckm = sckm
      out$rates = sckm$theta
      out$probs = rbeta( sckm$r, prior$prob$a*10, prior$prob$b*10)
      out$y = cbind(rbinom(n, out$nr[,1], out$p[1]), 
                    rbinom(n, out$nr[,2], out$p[2]))

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

