source("../../code/working/SIRD.R")

set.seed(1)
      
# Number of steps              
n = 52 # weekly for a year

# Initial state
X0 = c(16000,2,0,0)


# Parameters
ItoRrate = rgamma(2, .5*100,100)
R        = runif( 2,1,3)

probs    = c(.1,.5)
parms    = expand.grid(ItoRrate=ItoRrate,R=R,StoIprob=probs, ItoRprob=probs)

parms$StoIrate = apply(parms[,1:2],1,prod)

gammas = cbind(parms$StoIrate, parms$ItoRrate, 0, 0)
probs  = cbind(parms$StoIprob, parms$ItoRprob, 0, 0)


##################### Generate simulations ########################
sims = list()
for (i in 1:nrow(gammas)) sims[[i]] = SIRDsim(X0,gammas[i,],probs[i,],n)

save(sims,gammas,probs, file="SIRDsims.RData")

q("no")
