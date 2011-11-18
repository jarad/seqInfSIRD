# This will implement particle Gibbs sampling for 
#   S -> I -> R
# 
#####################################################################
# 
library(smcUtils)


p <- c(.5,.5)
g.true <- c(.3,.1)
dt <- 1

# Generate data
S0 <- 98
I0 <- 2
R0 <- 0
N <- S0+I0+R0

n <- 20

X <- matrix(NA,n+1,3, dimnames=list(1:(n+1),c("S","I","R")))
X[1,] <- c(S0,I0,R0)

y <- matrix(NA,n,2) 

oneStep <- function(X,g) {
  newI <- rpois(1, g[1]*X[1]*X[2]/N*dt)
  newR <- rpois(1, g[2]*X[2])

  newI <- min(X[1],      newI)
  newR <- min(X[2]+newI, newR)
 
 
  X[1] <- X[1]-newI
  X[2] <- X[2]+newI-newR
  X[3] <- X[3]     +newR

  y <- rep(NA,2)
  y[1] <- rbinom(1,newI,p[1])
  y[2] <- rbinom(1,newR,p[2])

  return(list(X=X,y=y))
}

for (i in 1:n) {
  tmp <- oneStep(X[i,],g.true)
  X[i+1,] <- tmp$X
  y[i,]   <- tmp$y
}


# MCMC
draw.g <- function(y,X, prior=matrix(1,2,2)) {
  newI <- -diff(X[,1])
  newR <-  diff(X[,3])

  n <- nrow(X)
  EnewI <- X[-n,1]*X[-n,2]/N
  EnewR <- X[-n,2]

  return(
    rgamma(2, prior[,1]+c(sum( newI), sum( newR)),
                 prior[,2]+c(sum(EnewI), sum(EnewR)))
  )
}

obsEqn <- function(y,X) {
  tmp <- sum(dbinom(y,X,p,log=T))
  if (is.nan(tmp)) stop()
  return(tmp)
}

draw.X <- function(y,X,g,J=nrow(y)*5) {
  n <- nrow(X)
  d <- ncol(X)
  Xp <- array(NA, dim=c(J,n,d))
  for (i in 1:d) Xp[,1,i] <- X[1,i]
  w <- matrix(1/J,J,n)
  A <- matrix(NA,J,n) # ancestors

  for (i in 2:n) {
    # Deterministically copy current particle
    Xp[1,i,] <- X[i,]
    w[1,i]   <- obsEqn(y[i-1,],Xp[1,i,])
    A[1,i] <- 1

    ks <- resample(w[,i-1], J, "stratified")$indices # one extra draw
    for (j in 2:J) {
      A[j,i] <- ks[j]
      tmp <- oneStep(Xp[ks[j],i-1,],g)
      Xp[j,i,] <- tmp$X
      w[j,i]   <- obsEqn(y[i-1,],Xp[j,i,])
    }
    w[,i] <- renormalize.weights(w[,i],log=T)
  }  
  k <- sample(J,1,prob=w[,n])
  Xc <- NA*X
  for (i in n:1) {
    Xc[i,] <- Xp[k,i,]
    k <- A[k,i]
  }
  return(Xc)
}


# Actual MCMC
n.reps <- 1e3
g.reps <- matrix(NA,n.reps,length(g.true))
X.reps <- array(NA, dim=c(n.reps, nrow(X), ncol(X)))
g.reps[1,]  <- g.true
X.reps[1,,] <- X


cat("  |",rep(" ",20),"|100%\n  |",sep="")
for (i in 2:n.reps) {
  if (i%%(n.reps/20)==0) cat("*")
  g.reps[i,] <- draw.g(y,X.reps[i-1,,])
  X.reps[i,,] <- draw.X(y,X.reps[i-1,,],g.reps[i,])
}
cat("|\n")

