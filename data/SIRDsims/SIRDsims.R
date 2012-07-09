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

thetas = cbind(parms$StoIrate, parms$ItoRrate, 0, 0)
probs  = cbind(parms$StoIprob, parms$ItoRprob, 0, 0)


##################### Generate simulations ########################
sims = list()
for (i in 1:nrow(thetas)) sims[[i]] = SIRDsim(X0,thetas[i,],probs[i,],n)

prior = list()
prior$theta = list()
prior$theta$a = c(150, 50, 0, 0)
prior$theta$b = c(100,100,10,10)

# Find the end of the outbreak
# defined as 4 consecutive weeks of zero S->I observed
for (i in 1:nrow(thetas)) {
  zeros = which(sims[[i]]$y$StoI==0)
  diff3 = diff(zeros,3)
  n[i] = min(zeros[match(3,diff3)]+3,52)
}

save(sims,thetas,probs,n,prior, file="SIRDsims.RData")

#q("no")
