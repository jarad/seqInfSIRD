library(tlpl)

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

sim = list()
param = list()

for (i in 1:n.sims)
{
  sckm$theta = rgamma(sckm$r, prior$rate$a, prior$rate$b)
  ### True states and transitions
  tl = tau_leap(sckm, n)

  ### Sample transitions
  p = rbeta(sckm$r, prior$prob$a, prior$prob$b)
  y = t(rbind(rbinom(n, tl$nr[,1], p[1]), rbinom(n, tl$nr[,2], p[2])))
  
  param[[i]] = list(rate=sckm$theta, prob=p)

  sim[[i]] = data.frame(sim=i, time=0:50, 
                  S=tl$X[,1], I=tl$X[,2], R=tl$X[,3],
                  StoI=c(tl$nr[,1], NA), ItoR=c(tl$nr[,2], NA),
                  yStoI=c(y[,1], NA), yItoR=c(y[,2], NA))
}

rm(i,tl,y,p)
save.image("SIR-sims.RData")

#q("no")


