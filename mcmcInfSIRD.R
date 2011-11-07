if (!is.element("rjags",installed.packages()[,1])) install.packages("rjags")
library(rjags)

# Simulate data
gamma <- rgamma(2,1); gamma <- c(.2,.1)
p <- rep(0.5,2)
N <- 100
I0 <- 2

n <- 20
x <- y <- dx <- matrix(NA,n,2)
X0 <- c(N-I0,I0)
dx[1,1] <- rpois(1, gamma[1]*X0[1]*X0[2]/N)
dx[1,2] <- rpois(1, gamma[2]*X0[2])
x[1,1] <- X0[1]-dx[1,1]
x[1,2] <- X0[2]+dx[1,1]-dx[1,2]
for (j in 1:2) {
  y[1,j] <- rbinom(1,dx[1,j],p[j])
}

for (i in 2:n) {
  dx[i,1] <- rpois(1, gamma[1]*x[i-1,1]*x[i-1,2]/N)
  dx[i,2] <- rpois(1, gamma[2]*x[i-1,2])
  x[i,1] <- x[i-1,1]-dx[i,1]
  x[i,2] <- x[i-1,2]+dx[i,1]-dx[i,2]
  for (j in 1:2) {
    y[i,j] <- rbinom(1,dx[i,j],p[j])
  }
}

plot(x[,2])


# Inference
mod <- jags.model("SIRD.txt", list(X0=X0, y=y, p=p, n=n, N=N),
                  list(dx=y)) 
res <- coda.samples(mod, c("gamma","x"), 1e3) 
