if (!is.element("rjags",installed.packages()[,1])) install.packages("rjags")
library(rjags)


######### Functions ########
fix.updates <- function(x,dx,p) {
  dx[1] <- min(c(dx[1],x[1]))
  if ( dx[2]+dx[4] > x[2]+dx[1]) {
    dx[2] <- rbinom(1, x[2]+dx[1], p)
    dx[4] <- x[2]+dx[1]-dx[2]
  }
  dx[4] <- x[2]+dx[1]-dx[2] 
  return(dx)
}

# Simulate data
gamma <- rgamma(4,1); gamma <- c(.2,.1,.05,.01)
p <- rep(0.5,4)
N <- 100
I0 <- 2

n <- 20
x <- y <- dx <- matrix(NA,n,4)
X0 <- c(N-I0,I0,0,0)
dx[1,1] <- rpois(1, gamma[1]*X0[1]*X0[2]/N)
dx[1,2] <- rpois(1, gamma[2]*X0[2])
dx[1,3] <- rpois(1, gamma[3]*X0[1])
dx[1,4] <- rpois(1, gamma[4]*X0[2])
dx[1, ] <- fix.updates(X0,dx[1,],gamma[2]*X0[2]/(gamma[2]*X0[2]+gamma[4]*X0[2]))

x[1,1] <- X0[1]-dx[1,1]        -dx[1,3]
x[1,2] <- X0[2]+dx[1,1]-dx[1,2]        -dx[1,4]
x[1,3] <- X0[3]                +dx[1,3]
x[1,4] <- X0[4]                        +dx[1,4]

for (j in 1:4) {
  y[1,j] <- rbinom(1,dx[1,j],p[j])
}

for (i in 2:n) {
  dx[i,1] <- rpois(1, gamma[1]*x[i-1,1]*x[i-1,2]/N)
  dx[i,2] <- rpois(1, gamma[2]*x[i-1,2])
  dx[i,3] <- rpois(1, gamma[3]*x[i-1,1])
  dx[i,4] <- rpois(1, gamma[4]*x[i-1,2])
  dx[i,]  <- fix.updates(x[i-1,], dx[i,], 
                         gamma[2]*x[i-1,2]/(gamma[2]*x[i-1,2]+gamma[4]*x[i-1,2]))

  x[i,1] <- x[i-1,1]-dx[i,1]        -dx[i,3]
  x[i,2] <- x[i-1,2]+dx[i,1]-dx[i,2]        -dx[i,4]
  x[i,3] <- x[i-1,3]                +dx[i,3]
  x[i,4] <- x[i-1,4]                        +dx[i,4]

  for (j in 1:4) {
    y[i,j] <- rbinom(1,dx[i,j],p[j])
  }
}

plot(x[,2])


# Inference
mod <- jags.model("SIRD.txt", list(X0=X0, y=y, p=p, n=nrow(y), N=sum(X0)),
                  list(dx=y)) 
res <- coda.samples(mod, c("gamma","x"), 1e3) 
