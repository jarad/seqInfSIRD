# Particle Gibbs samping for random walk with noise

# Simulate data
n <- 20
y <- x <- rep(NA,n)
x[1] <- 0
y[1] <- rnorm(1,x[1])
for (i in 2:n) {
  x[i] <- rnorm(1,x[i-1])
  y[i] <- rnorm(1,x[i])
}

# Evolution variance is unknown
#   prior is Unif(0,10) on standard deviation
draw.sigma <- function(x) {
  e <- diff(x)
  repeat {
    s <- sqrt(1/rgamma(1, length(e)/2, sum(e^2)/2))
    if (s<10) return(s)
  }
}

draw.x <- function(y,x,s,J=length(y)*5) {
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
      xp[j,i] <- rnorm(1,xp[A[j,i],i-1],s)
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
n.reps <- 1e4
x.reps <- matrix(NA,n.reps,n); x.reps[1,] <- x 
s.reps <- rep(NA,n.reps);      s.reps[1] <- 1

cat("  |",rep(" ",20),"|100%\n  |",sep="")
for (i in 2:n.reps) {
  f (i%%(n.reps/20)==0) cat("*")
  s.reps[i] <- draw.sigma(x.reps[i-1,])
  x.reps[i,] <- draw.x(y,  x.reps[i-1,], s.reps[i])
}
cat("|\n")



