library(tlpl)

# Generate data
## Set up SEIR model

cc = 1
sys = sckm("seir", X=c(16000,10,10,0)/cc)
N = sum(sys$X)

# Read data
d = read.csv("harareClean.csv")[-61,]
n = nrow(d)

n = 26



y = matrix(0, nrow=n, ncol=3)
y[,2] = d$diff2[1:n]
sys$theta = rep(0,sys$r)

np = 1e4
prior = list(prob=list(a=c(1,10,1)/cc, b=c(1e5,990,1e5)/cc),
             rate=list(a=rep(1,3), b=rep(1,3)),
             X = rmultinom(np, sum(sys$X), sys$X/sum(sys$X)))
for (i in 1:ncol(prior$X)) 
{
  N = rbinom(1,1.5e6, 0.01)
  prior$X[,i] = rmultinom(1,N, c(.998,.001,.001,0))
}


res = tlpl(list(y=y, tau=1), sckm=sys, prior=prior, n.particles=np, verbose=1,
           method="stratified")
q   = tlpl_quantile(res, c(.025,.975), verbose=1)


prior$rate = list(a=c(1,1e5,1), b=c(1,1e5,1))
res2 = tlpl(list(y=y, tau=1), sckm=sys, prior=prior, n.particles=np, verbose=1,
           method="stratified")
q2   = tlpl_quantile(res2, c(.025,.975), verbose=1)



pdf("harare-fit.pdf", width=10)
par(mfrow=c(2,max(sys$s,sys$r)))
# States 
for (j in 1:sys$s)  {
  plot(0:n,0:n, type="n", ylim=c(0,1), xlim=c(0,n), main=sys$states[j], ylab="Proportion", xlab="Week")
  lines(0:n, q$X.quantiles[j,1,]/N, col=2, lwd=2)
  lines(0:n, q$X.quantiles[j,2,]/N, col=2, lwd=2)
  lines(0:n, q2$X.quantiles[j,1,]/N, col=4, lwd=2)
  lines(0:n, q2$X.quantiles[j,2,]/N, col=4, lwd=2)
}

# Probabilities
for (j in 2)
{
  plot(0,0, type="n", ylim=range(q$p.quantiles[j,,]), xlim=c(0,n), main=paste("p:",sys$reactions[j]), ylab="", xlab="Week")   
  lines(0:n, q$p.quantiles[j,1,], col=2, lwd=2)
  lines(0:n, q$p.quantiles[j,2,], col=2, lwd=2)
  lines(0:n, q2$p.quantiles[j,1,], col=4, lwd=2)
  lines(0:n, q2$p.quantiles[j,2,], col=4, lwd=2)
}

#plot(0,0,type="n", axes=F, xlab="", ylab="")

# Rates
for (j in 1:3)
{
  plot(0,0, type="n", ylim=range(q$r.quantiles[j,,]), xlim=c(0,n), main=paste("r:",sys$reactions[j]), ylab="", xlab="Week")
  lines(0:n, q$r.quantiles[j,1,], col=2, lwd=2)
  lines(0:n, q$r.quantiles[j,2,], col=2, lwd=2)
  lines(0:n, q2$r.quantiles[j,1,], col=4, lwd=2)
  lines(0:n, q2$r.quantiles[j,2,], col=4, lwd=2)
}
dev.off()
  

