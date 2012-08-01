# This file no longer works as a test file as of 2012/6/13 - JBN


source("inference.R")

n = list()
n$reps = 9
n$states = 4
n$steps = 3
X0 = c(16000,2,0,0)
N = sum(X0)
X = matrix(X0,nrow=n$states,ncol=n$reps)
theta = matrix(rgamma(N.RXNS*n$reps,1),N.RXNS,n$reps)
prop = matrix(rbeta(N.RXNS*n$reps,1,1),N.RXNS,n$reps)
hyper = array(rgamma(N.RXNS*n$reps*2,1),dim=c(N.RXNS,2,n$reps))

# Check hazard


# Check simulate
set.seed(1)
sims = simulate.R(X, theta, n.steps=3)
set.seed(1)
stopifnot(all.equal(sims,simulate.C(X, theta, n.steps=n$steps)))

# Check inference.one.step
for (i in 1:n$reps) {
  inf.R = inference.one.step.R(sims$X[,1,i],sims$dX[,2,i],N,prop[,i],hyper[,,i])
  inf.C = inference.one.step.C(sims$X[,1,i],sims$dX[,2,i],N,prop[,i],hyper[,,i])
  stopifnot(all.equal(inf.C,inf.R))
}





# Check one.step
set.seed(1)
os = one.step.R(X,hyper,theta,prop)
set.seed(1)
stopifnot(all.equal(os,os.C <- one.step.C(X,hyper,theta,prop)))















########################################################
# Old stuff
########################################################

#source("plot-functions.R")

# Zimbabwe measles data
#n.deaths <- 631
#n.cases  <- 13783
#n.confirmed <- 693 

#I0 <- 1
#N <- 20000
#theta <- c(1, .4, 0, n.deaths/n.cases*.4)


#X <- matrix(c(N-I0, I0, 0, 0), nrow=4, ncol=2)
#out <- simulate.R(X, theta, sum(X), 50)

#inf <- inference.R(out$X[,,1], out$dX[,,1], N, prop=rep(1,4))
#par(mfrow=c(2,3))
#plot.simulation(out$X[2,,])
#plot(0,0,type='n', axes=F, xlab='', ylab='')
#for (i in 1:N.RXNS) plot.inference(inf[,i,])

# Multiple starting states and hyper parameters
#n.reps <- 1e4
#I0     <- rpois(n.reps, 2)+1
#N      <- 20000
#X      <- cbind(N-I0,I0,0,0)
#hyper  <- array( rgamma(n.reps*N.RXNS*2,1), dim=c(n.reps, N.RXNS, 2))
#theta  <- matrix(rgamma(n.reps*N.RXNS  ,1),       n.reps, N.RXNS)
#out    <- one.step.R(X, hyper, theta)


# Timing comparison
#n.times <- 1e2 # this will take about 3 mins
#C.times <- R.times <- rep(NA,n.times)
#for (i in 1:n.times) {
#  n.reps <- 1e4
#  I0     <- rpois(n.reps, 2)+1
#  N      <- 20000
#  X      <- cbind(N-I0,I0,0,0)
#  hyper  <- array( rgamma(n.reps*N.RXNS*2,1), dim=c(n.reps, N.RXNS, 2))
#  theta  <- matrix(rgamma(n.reps*N.RXNS  ,1),       n.reps, N.RXNS)
#  prop = matrix(1,N.RXNS,n.reps)
#  R.times[i] <- system.time(one.step.R(X,hyper,theta,prop=prop))[3]
#  C.times[i] <- system.time(one.step.C(X,hyper,theta))[3]
#}
#c(mean(R.times), sd(R.times))
#c(mean(C.times), sd(C.times))
#mean(R.times)/mean(C.times)
#sd(R.times/C.times)

                      
