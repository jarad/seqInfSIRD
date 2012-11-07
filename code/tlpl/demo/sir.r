
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
sckm$theta = c(0.5/N,0.25)

## Simulate data
set.seed(2)
n = 50

### True states and transitions
tl = tau.leap(sckm, n)

### Sample transitions
p = c(0.5,0.5) # Sample probabilities for S->I and I->R respectively
y = cbind(rbinom(n, tl$nr[,1], p[1]), rbinom(n, tl$nr[,2], p[2]))

# Perform inference
cat("Running sequential inference...\n")
sckm$theta[1] = sckm$theta[1]*N
prior = tlpl.prior(sckm$X, 1e1, 1e1, sckm$theta*2, 2, sckm$r)
z = tlpl(list(y=y, tau=1), sckm, prior=prior, n.particles=1e4, mult=c(1/N,1), verbose=T)
qs = tlpl_quantile(z)

# Make figures
ld = 2
clrs = c("green","red","blue")

## Data

plot( tl$X[,1], type="l", ylim=c(0,sckm$X[1]), lwd=ld, col=clrs[1], xlab="Time", ylab="Number")
lines(tl$X[,2], lwd=ld, col=clrs[2])
lines(tl$X[,3], lwd=ld, col=clrs[3])
legend("right", c("S","I","R"), col=clrs, lwd=ld)


## States
xx = 0:n
par(mfrow=c(1,3))

### Susceptible
plot( xx, qs$X[,1,1], type="l", ylim=c(0, max(sckm$X)), main="Susceptible", ylab="Count", xlab="Time")
lines(xx, qs$X[,1,2], lwd=2)
lines(xx, qs$X[,1,3])
lines(xx, tl$X[,1], col="red")

### Infecteds
plot( xx, qs$X[,2,1], type="l", ylim=c(0, max(sckm$X)), main="Infected", ylab="Count", xlab="Time")
lines(xx, qs$X[,2,2], lwd=2)
lines(xx, qs$X[,2,3])
lines(xx, tl$X[,2], col="red")


### Recovered
plot( xx, qs$X[,3,1], type="l", ylim=c(0, max(sckm$X)), main="Recovered", ylab="Count", xlab="Time")
lines(xx, qs$X[,3,2], lwd=2)
lines(xx, qs$X[,3,3])
lines(xx, tl$X[,3], col="red")


## Sampling probabilities and reaction rates
par(mfrow=c(2,2))

### S->I probability
plot( xx, qs$p[,1,1], type="l", ylim=range(qs$p[,1,]), main="S -> I", ylab="Probability", xlab="Time")
lines(xx, qs$p[,1,2], lwd=2)
lines(xx, qs$p[,1,3])
abline(h=p[1], col="red")

### I->R probability
plot( xx, qs$p[,2,1], type="l", ylim=range(qs$p[,2,]), main="I -> R", ylab="Probability", xlab="Time")
lines(xx, qs$p[,2,2], lwd=2)
lines(xx, qs$p[,2,3])
abline(h=p[2], col="red")


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






