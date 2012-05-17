source("../../code/working/SIRD.R")

#### Simulation parameter settings
# Poisson transition parameters (gamma)
prior.g = rbind(c(1,1) # S -> I
                c(1,1) # I -> R
                c(1,1) # S -> R
                c(1,1) # I -> D)
                
# Binomial observation parameters (beta)           
prior.p = rbind(c(1,1) # S -> I
                c(1,1) # I -> R
                c(1,1) # S -> R
                c(1,1) # I -> D)
      
# Number of steps              
n = 52 # weekly for a year

# Initial state
X0 = c(16000,2,0,0)

# Number of simulations

n.sims = 2^4



#### Gridded draws
n.sims

