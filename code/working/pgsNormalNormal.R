# Particle Gibbs samping for random walk with noise
library(smcUtils)

# Simulate data
n <- 20
mu <- 0
x <- rnorm(n,mu)
y <- rnorm(n,x)

source("pgs.R")

out <- pgs(y, 
           function(y,x,theta) dnorm(y,x,log=T),             # y ~ N(x,1)
           function(x,mu) rnorm(length(x), mu),        # x ~ N(mu,1)
           function(x) rnorm(1, mean(x), 1/sqrt(length(x))), # mu|x 
           function(n,mu) rnorm(n,mu),                 
           list(x=x*0, theta=mu), n.iter=1e3)


hist(out$theta,100,freq=F)
curve(dnorm(x, mean(y), sqrt(2/n)), col="red", add=T)

