if (!is.element("rjags",installed.packages()[,1])) install.packages("rjags")
library(rjags)

# Simulate data
gamma <- 0.01
p <- 0.5
N <- 100

n <- 20
x <- y <-dx <- rep(NA,n)
for (i in 1:n) {
  dx[i] = rpois(1,gamma*ifelse(i==1,N,x[i-1]))
  x[i]  = ifelse(i==1,N,x[i-1])-dx[i]
  y[i]  = rbinom(1,dx[i],p)
}

plot(x)


# Inference
mod <- jags.model("SIRD.txt", list(S0=N, y=y, p=p, n=n),
                  list(dx=y)) 
res <- coda.samples(mod, c("gamma","x"), 1e3) 
