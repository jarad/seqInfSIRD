sckm = list()
sckm$s = 3 # species (S,I,R)
sckm$r = 2 # reactions (S->I, I->R)
#                   S -> I    I -> R
sckm$Pre  = rbind( c(1,1,0), c(0,1,0))
sckm$Post = rbind( c(0,2,0), c(0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(100,2,0)
sckm$theta = c(2/100,1)

set.seed(1)
n = 15
y = tau.leap(sckm, n)


cat("Running sequential inference...\n")
prior = tlpl.prior(sckm$X, 1e4, 1, sckm$theta*10, 10, sckm$r)
z = tlpl(list(y=y$nr, tau=1), sckm, prior=prior, n.particles=100)

qs = tlpl_quantile(z)

# Make figures
ld = 2
clrs = c("green","red","blue")
plot(y$X[,1], type="l", ylim=c(0,sckm$X[1]), lwd=ld, col=clrs[1], xlab="Time", ylab="Number")
lines(y$X[,2], lwd=ld, col=clrs[2])
lines(y$X[,3], lwd=ld, col=clrs[3])
legend("right", c("S","I","R"), col=clrs, lwd=ld)


## States
xx = 0:n
par(mfrow=c(1,3))

### Susceptible
plot( xx, qs$X[,1,1], type="l", ylim=c(0, max(sckm$X)), main="Susceptible", ylab="Count", xlab="Time")
lines(xx, qs$X[,1,2], lwd=2)
lines(xx, qs$X[,1,3])

### Infecteds
plot( xx, qs$X[,2,1], type="l", ylim=c(0, max(sckm$X)), main="Infected", ylab="Count", xlab="Time")
lines(xx, qs$X[,2,2], lwd=2)
lines(xx, qs$X[,2,3])

### Recovered
plot( xx, qs$X[,3,1], type="l", ylim=c(0, max(sckm$X)), main="Recovered", ylab="Count", xlab="Time")
lines(xx, qs$X[,3,2], lwd=2)
lines(xx, qs$X[,3,3])


## Sampling probabilities and reaction rates
par(mfrow=c(2,2))

### S->I probability
plot( xx, qs$p[,1,1], type="l", ylim=range(qs$p[,1,]), main="S -> I", ylab="Probability", xlab="Time")
lines(xx, qs$p[,1,2], lwd=2)
lines(xx, qs$p[,1,3])
abline(h=1, col="red")

### I->R probability
plot( xx, qs$p[,2,1], type="l", ylim=range(qs$p[,2,]), main="I -> R", ylab="Probability", xlab="Time")
lines(xx, qs$p[,2,2], lwd=2)
lines(xx, qs$p[,2,3])
abline(h=1, col="red")


### S->I rate
plot( xx, qs$r[,1,1], type="l", ylim=range(qs$r[,1,]), main="S -> I", ylab="Rate", xlab="Time")
lines(xx, qs$r[,1,2], lwd=2)
lines(xx, qs$r[,1,3])
abline(h=sckm$theta[1], col="red")


### I->R rate
plot( xx, qs$r[,2,1], type="l", ylim=range(qs$r[,2,]), main="I -> R", ylab="Rate", xlab="Time")
lines(xx, qs$r[,2,2], lwd=2)
lines(xx, qs$r[,2,3])
abline(h=sckm$theta[2], col="red")






