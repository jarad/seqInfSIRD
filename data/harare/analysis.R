library(tlpl)

# Generate data
## Set up SEIR model

sys = sckm("seir", X=c(16000,10,10,0))

# Read data
d = read.csv("harareClean.csv")
d$new = d$diff2
n = nrow(d)


# Simulate data
sys$theta = c(1,1,.5)
out = tau_leap(sys,60)
y_tmp = rbinom(n, out$nr[,2], .02)
#plot(y_tmp)



y = matrix(0, nrow=3, ncol=n)
y[2,] = d$new
n = 20
y = y[,1:n]
sys$theta = rep(0,sys$r)

prior = list(prob=list(a=rep(1,2,1), b=c(1e5,998,1e5)),
             rate=list(a=c(10,1e5,5)/2, b=c(10,1e5,10)/2),
             X = sys$X)

res = tlpl(list(y=y, tau=1), sys=sys, prior=prior, n.particles=1e4, verbose=1)
q   = tlpl_quantile(res, c(.025,.975), verbose=1)


par(mfrow=c(2,max(sys$s,sys$r)))
# States 
for (j in 1:sys$s)  {
  plot(0:n,0:n, type="n", ylim=range(q$X.quantiles[j,,]), xlim=c(0,n), main=sys$states[j], ylab="Number", xlab="Week")
  lines(0:n, q$X.quantiles[j,1,], col=2)
  lines(0:n, q$X.quantiles[j,2,], col=2)
}

# Probabilities
for (j in 2)
{
  plot(0,0, type="n", ylim=range(q$p.quantiles[j,,]), xlim=c(0,n), main=paste("p:",sys$rxns[j]), ylab="", xlab="Week")   
  lines(0:n, q$p.quantiles[j,1,], col=2)
  lines(0:n, q$p.quantiles[j,2,], col=2)
}

plot(0,0,type="n", axes=F, xlab="", ylab="")

# Rates
for (j in c(1,3))
{
  plot(0,0, type="n", ylim=range(q$r.quantiles[j,,]), xlim=c(0,n), main=paste("r:",sys$rxns[j]), ylab="", xlab="Week")
  lines(0:n, q$r.quantiles[j,1,], col=2)
  lines(0:n, q$r.quantiles[j,2,], col=2)
}

  

