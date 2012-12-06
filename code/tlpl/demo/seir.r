plots = FALSE
ask   = TRUE

# Generate data
## Set up SIR model
sckm = list()
sckm$s = 4 # species (S,E,I,R)
sckm$r = 3 # reactions (S->I, I->R)
#                    S -> E      E -> I      I -> R
sckm$Pre  = rbind( c(1,0,1,0), c(0,1,0,0), c(0,0,1,0))
sckm$Post = rbind( c(0,1,1,0), c(0,0,1,0), c(0,0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(16000,0,100,0)
N = sum(sckm$X)
sckm$theta = c(0.5,0.25,0.25)
sckm$lmult = log(c(1/N,1,1))

## Simulate data
set.seed(2)
n = 100

### True states and transitions
tl = tau.leap(sckm, n)

### Sample transitions
p = c(0,.05,0) # Sample probabilities for S->I, E->I, and I->R respectively
y = rbind(rbinom(n, tl$nr[,1], p[1]), 
          rbinom(n, tl$nr[,2], p[2]),
          rbinom(n, tl$nr[,3], p[3]))


clrs = c("seagreen","pink","red","blue")
ld=3
par(mfrow=c(1,3))
plot( tl$X[,1], type="l", ylim=c(0,sckm$X[1]), lwd=ld, col=clrs[1], 
      xlab="Time", ylab="Number", main="Truth")
lines(tl$X[,2], lwd=ld, col=clrs[2])
lines(tl$X[,3], lwd=ld, col=clrs[3])
lines(tl$X[,4], lwd=ld, col=clrs[4])
legend("topright", c("S","E","I","R"), col=clrs, lwd=ld)




plot( y[1,], type="p", pch=19, ylim=range(y), lwd=ld, col=clrs[2], 
      xlab="Time", ylab="Number", main="Observations")
points(y[2,], lwd=ld, col=clrs[3], pch=19)
points(y[3,], lwd=ld, col=clrs[4], pch=19)
legend("topright", c("S->E","E->I","I->R"), col=clrs[2:4], lwd=ld)



### Cumulative transitions
y2 = t(apply(y,1,cumsum))
plot( y2[1,], type="l", ylim=range(y2), lwd=ld, col=clrs[2], 
      xlab="Time", ylab="Number", main="Cumulative Observations")
lines(y2[2,], lwd=ld, col=clrs[3])
lines(y2[3,], lwd=ld, col=clrs[4])
legend("topleft", c("S->E","E->I","I->R"), col=clrs[2:4], lwd=ld)


if (plots) dev.copy2pdf(file="example2-data.pdf")

if (ask) readline("Hit <enter> to continue:")


# Perform inference
cat("\nRunning sequential inference...\n")
prior = tlpl.prior(sckm$X, 1e6*(p+.0001), 1e6*(1-(p+.0001)), sckm$theta*1e6, 1e6, sckm$r)
z = tlpl(list(y=y, tau=1), sckm, prior=prior, n.particles=1e4, engine="C", verbose=1)

cat("\nCalculating quantiles...\n")
qs = tlpl_quantile(z)




# Make figures
ld = 2

## States
xx = 0:n
par(mfrow=c(2,2))

### Susceptible
plot( xx, qs$X[1,1,], type="l", ylim=c(0, max(sckm$X)), main="Susceptible", ylab="Count", xlab="Time")
lines(xx, qs$X[1,2,], lwd=ld)
lines(xx, qs$X[1,3,])
lines(xx, tl$X[,1], col="red", lwd=ld)

### Exposed
plot( xx, qs$X[2,1,], type="l", ylim=c(0, max(sckm$X)), main="Exposed", ylab="Count", xlab="Time")
lines(xx, qs$X[2,2,], lwd=ld)
lines(xx, qs$X[2,3,])
lines(xx, tl$X[,2], col="red", lwd=ld)

### Infected
plot( xx, qs$X[3,1,], type="l", ylim=c(0, max(sckm$X)), main="Infected", ylab="Count", xlab="Time")
lines(xx, qs$X[3,2,], lwd=ld)
lines(xx, qs$X[3,3,])
lines(xx, tl$X[,3], col="red", lwd=ld)


### Recovered
plot( xx, qs$X[4,1,], type="l", ylim=c(0, max(sckm$X)), main="Recovered", ylab="Count", xlab="Time")
lines(xx, qs$X[4,2,], lwd=ld)
lines(xx, qs$X[4,3,])
lines(xx, tl$X[,4], col="red", lwd=ld)


if (plots) dev.copy2pdf(file="example2-states.pdf")

if (ask) readline("Hit <enter> to continue:")


## Sampling probabilities and reaction rates
par(mfrow=c(2,3))

### S->E probability
plot( xx, qs$p[1,1,], type="l", ylim=range(qs$p[1,,]), 
      main="S -> E", ylab="Probability", xlab="Time")
lines(xx, qs$p[1,2,], lwd=ld)
lines(xx, qs$p[1,3,])
abline(h=p[1], col="red", lwd=ld)

### E->I probability
plot( xx, qs$p[2,1,], type="l", ylim=range(qs$p[2,,]), main="E -> I", ylab="Probability", xlab="Time")
lines(xx, qs$p[2,2,], lwd=ld)
lines(xx, qs$p[2,3,])
abline(h=p[2], col="red", lwd=ld)

### I->R probability
plot( xx, qs$p[3,1,], type="l", ylim=range(qs$p[3,,]), main="I -> R", ylab="Probability", xlab="Time")
lines(xx, qs$p[3,2,], lwd=ld)
lines(xx, qs$p[3,3,])
abline(h=p[3], col="red", lwd=ld)





### S->E rate
plot( xx, qs$r[1,1,], type="l", ylim=range(qs$r[1,,]), main="S -> E", ylab="Rate", xlab="Time")
lines(xx, qs$r[1,2,], lwd=ld)
lines(xx, qs$r[1,3,])
abline(h=sckm$theta[1], col="red", lwd=ld)


### E->I rate
plot( xx, qs$r[2,1,], type="l", ylim=range(qs$r[2,,]), main="E -> I", ylab="Rate", xlab="Time")
lines(xx, qs$r[2,2,], lwd=ld)
lines(xx, qs$r[2,3,])
abline(h=sckm$theta[2], col="red", lwd=ld)


### I->R rate
plot( xx, qs$r[3,1,], type="l", ylim=range(qs$r[3,,]), main="I -> R", ylab="Rate", xlab="Time")
lines(xx, qs$r[3,2,], lwd=ld)
lines(xx, qs$r[3,3,])
abline(h=sckm$theta[3], col="red", lwd=ld)



if (plots) dev.copy2pdf(file="example2-parameters.pdf")




# Predictions
tt = 30
np = 1000
z2 = list(X = z$X[,1:np,tt])
z2$hyper = list()
z2$hyper$rate = list(a = z$hyper$rate$a[,1:np,tt], 
                     b = z$hyper$rate$b[,1:np,tt])
z2$hyper$prob = list(a = z$hyper$prob$a[,1:np,tt], 
                     b = z$hyper$prob$b[,1:np,tt])

z3 = tlpl_predict(sckm, n-tt, z2, verbose=1)

Xq = tlpl_quantile(z3, which="x")



## States
par(mfrow=c(2,2))
x1 = xx[1:tt]
x2 = xx[tt:n]

### Susceptible
plot( x1, qs$X[1,1,1:tt], type="l", ylim=c(0, max(sckm$X)), xlim=c(0,n), main="Susceptible", ylab="Count", xlab="Time")
lines(x1, qs$X[1,2,1:tt], lwd=ld)
lines(x1, qs$X[1,3,1:tt])

lines(x2, Xq$X.quantiles[1,1,], col="blue")
lines(x2, Xq$X.quantiles[1,2,], lwd=ld, col="blue")
lines(x2, Xq$X.quantiles[1,3,], col="blue")

lines(xx, tl$X[,1], col="red", lwd=ld)


### Exposed
plot( x1, qs$X[2,1,1:tt], type="l", ylim=c(0, max(sckm$X)), xlim=c(0,n), main="Exposed", ylab="Count", xlab="Time")
lines(x1, qs$X[2,2,1:tt], lwd=ld)
lines(x1, qs$X[2,3,1:tt])

lines(x2, Xq$X.quantiles[2,1,], col="blue")
lines(x2, Xq$X.quantiles[2,2,], lwd=ld, col="blue")
lines(x2, Xq$X.quantiles[2,3,], col="blue")

lines(xx, tl$X[,2], col="red", lwd=ld)


### Infected
plot( x1, qs$X[3,1,1:tt], type="l", ylim=c(0, max(sckm$X)), xlim=c(0,n), main="Infected", ylab="Count", xlab="Time")
lines(x1, qs$X[3,2,1:tt], lwd=ld)
lines(x1, qs$X[3,3,1:tt])

lines(x2, Xq$X.quantiles[3,1,], col="blue")
lines(x2, Xq$X.quantiles[3,2,], lwd=ld, col="blue")
lines(x2, Xq$X.quantiles[3,3,], col="blue")

lines(xx, tl$X[,3], col="red", lwd=ld)


### Recovered
plot( x1, qs$X[4,1,1:tt], type="l", ylim=c(0, max(sckm$X)), xlim=c(0,n), main="Recovered", ylab="Count", xlab="Time")
lines(x1, qs$X[4,2,1:tt], lwd=ld)
lines(x1, qs$X[4,3,1:tt])

lines(x2, Xq$X.quantiles[4,1,], col="blue")
lines(x2, Xq$X.quantiles[4,2,], lwd=ld, col="blue")
lines(x2, Xq$X.quantiles[4,3,], col="blue")

lines(xx, tl$X[,4], col="red", lwd=ld)



if (plots) dev.copy2pdf(file="example2-prediction")





# Predicting data

# need to get particle estimate of number of reactions before time 'tt'
# this should be returned from tlpl
# currently cheating
ycusum = apply(z3$y[2,,],1,cumsum)
y2q = apply(ycusum, 1, function(x) quantile(x, c(.025,.5,.975)))

plot(x2[-1], y2q[1,], type="l", ylim=c(0,max(y2q)), xlim=c(0,n), main="Cumulative E->I", ylab="Count", xlab="Time")
lines(x2[-1], y2q[2,], lwd=2)
lines(x2[-1], y2q[3,])

lines(xx[-1], y2[2,]-y2[2,tt], col="red", lwd=ld)

