# Particle Gibbs samping for random walk with noise
library(smcUtils)

# Simulate data
n <- 20
mu <- 0
x <- rnorm(n,mu)
y <- rnorm(n,x)

source("pgs.R")

out <- pgs(y, 
           function(y,x,theta) dnorm(y,x,log=T), 
           function(x,theta) rnorm(length(x), theta),
           function(x) rnorm(1, mean(x), 1/sqrt(length(x))),
           function(n,theta) rep(0,n),
           list(x=x*0, theta=mu), n.iter=1e3)


hist(out$theta,100,freq=F)
curve(dnorm(x, mean(y), 1/sqrt(n)), col="red", add=T)

