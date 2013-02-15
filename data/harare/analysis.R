library(tlpl)

# Generate data
## Set up SEIR model
sckm = list()
sckm$s = 4 # species (S,E,I,R)
sckm$r = 3 # reactions (S->E, E->I, I->R)
#                    S -> E      E -> I      I -> R
sckm$Pre  = rbind( c(1,0,1,0), c(0,1,0,0), c(0,0,1,0))
sckm$Post = rbind( c(0,1,1,0), c(0,0,1,0), c(0,0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(16000,0,10,0)
N = sum(sckm$X)
sckm$lmult = log(c(1/N,1,1))
sckm$states = c("S","E","I","R")
sckm$rxns = c("S->E","E->I","I->R")

# Read data
d = read.csv("harareClean.csv")
d$new = diff(c(0,d$total))
n = nrow(d)
y = matrix(0, nrow=3, ncol=n)
y[2,] = d$new

library(tlpl)
sckm$theta = rep(0,sckm$r)

prior = list(prob=list(a=rep(1,5,1), b=c(1e5,95,1e5)),
             rate=list(a=c(10,1e5,5), b=c(10,1e5,10)),
             X = sckm$X)

res = tlpl(list(y=y, tau=1), sckm=sckm, prior=prior, n.particles=10000, verbose=1)
q   = tlpl_quantile(res, c(.025,.975), verbose=1)


par(mfrow=c(3,max(sckm$s,sckm$r)))
# States 
for (j in 1:sckm$s)  {
  plot(0:n,0:n, type="n", ylim=range(q$X.quantiles[j,,]), xlim=c(0,n), main=sckm$states[j])
  lines(0:n, q$X.quantiles[j,1,], col=2)
  lines(0:n, q$X.quantiles[j,2,], col=2)
}

# Probabilities
for (j in 1:sckm$r)
{
  plot(0,0, type="n", ylim=range(q$p.quantiles[j,,]), xlim=c(0,n), main=paste("p:",sckm$rxns[j]))   
  lines(0:n, q$p.quantiles[j,1,], col=2)
  lines(0:n, q$p.quantiles[j,2,], col=2)
}

plot(0,0,type="n", axes=F, xlab="", ylab="")

# Rates
for (j in 1:sckm$r)
{

  plot(0,0, type="n", ylim=range(q$r.quantiles[j,,]), xlim=c(0,n), main=paste("r:",sckm$rxns[j]))
  lines(0:n, q$r.quantiles[j,1,], col=2)
  lines(0:n, q$r.quantiles[j,2,], col=2)
}

  

