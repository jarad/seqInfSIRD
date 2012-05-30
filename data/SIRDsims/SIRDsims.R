source("../../code/working/SIRD.R")

#### Simulation parameter settings
# Poisson transition parameters (gamma)
prior.g = rbind(c(1.5*10,10), # S -> I
                c( .5*10,10), # I -> R
                c(.01*10,10), # S -> R
                c(.05*10,10)) # I -> D
                
# Binomial observation parameters (beta)           
prior.p = rbind(c(1,1), # S -> I
                c(1,1), # I -> R
                c(1,1), # S -> R
                c(1,1)) # I -> D
      
# Number of steps              
n = 52 # weekly for a year

# Initial state
X0 = c(16000,2,0,0)

# Number of simulations
n.sims = 2^4 

# Whether the draws should be random from the prior or gridded
random = TRUE 

#### Random draws
if (random) {
  gammas = probs = matrix(NA,n.sims,4)
  for (i in 1:4) {
    gammas[,i] = rgamma(n.sims, prior.g[i,1], prior.g[i,2])	
     probs[,i] =  rbeta(n.sims, prior.p[i,1], prior.p[i,2])	
  }
}
gammas[,3] = 0

#### Gridded draws
if (!random) {
  n.each.dim = floor(n.sims^.25)
  qts = (1:n.each.dim)/(n.each.dim+1)	
  gammas = probs = matrix(NA,n.sims,4)
  for (i in 1:4) {
    gammas[,i] = qgamma(qts, prior.g[i,1], prior.g[i,2])	
     probs[,i] =  qbeta(qts, prior.p[i,1], prior.p[i,2])	
  }
}

# Need to expand.grid for the gridded draws

##################### Generate simulations ########################
set.seed(1)
sims = list()
for (i in 1:n.sims) sims[[i]] = SIRDsim(X0,gammas[i,],probs[i,],n)

save(sims,gammas,probs, file="SIRDsims.RData")

q("no")
