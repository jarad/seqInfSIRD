# Particle Gibbs samping for random walk with noise
library(smcUtils)

# Simulate data
n <- 20
mu <- 0
x <- rnorm(n,mu)
y <- rnorm(n,x)

# Evolution variance is unknown
#   prior is Unif(0,10) on standard deviation
draw.mu <- function(x) {
  rnorm(1, mean(x), 1/sqrt(length(x)))
}

draw.x <- function(y,x,mu,J=length(y)*5) {
  w <- xp <- A <- matrix(NA,J,length(x))
  xp[,1] <- 0 # intial x is zero
  w[,1] <- 1/J
  
  # Deterministically copy x
  xp[1,] <- x  
  A[1,] <- 1

  ks <- rep(NA,J)
  for (i in 2:n) {
    w[1,i] <-  dnorm(y[i], xp[1,i], log=T)
    ks[-1] <- resample(w[,i-1],J-1,"stratified")$indices
    for (j in 2:J) {
      A[j,i] <- ks[j]
      xp[j,i] <- rnorm(1,mu)
      w[j,i] <- dnorm(y[i], xp[j,i], log=T)
    }
    w[,i] <- renormalize.weights(w[,i], log=T)
  }
  
  # Sample a path
  xc <- NA*x
  k <- sample(J,1,prob=w[,i])
  for (i in n:1) {
    xc[i] <- xp[k,i]
    k <- A[k,i]
  }
  return(xc)
}



## MCMC
n.reps <- 1e3
x.reps <- matrix(NA,n.reps,n); x.reps[1,] <- x 
mu.reps <- rep(NA,n.reps);     mu.reps[1] <- 0

cat("  |",rep(" ",20),"|100%\n  |",sep="")
for (i in 2:n.reps) {
  if (i%%(n.reps/20)==0) cat("*")
  mu.reps[i] <- draw.mu( x.reps[i-1,])
  x.reps[i,] <- draw.x(y,x.reps[i-1,], mu.reps[i])
}
cat("|\n")

hist(mu.reps,100,freq=F)
curve(dnorm(x, mean(y), 1/sqrt(n)), col="red", add=T)

