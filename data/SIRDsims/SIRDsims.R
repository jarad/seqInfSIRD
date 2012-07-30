source("../../code/working/SIRD.R")


set.seed(1)
      
# Number of steps              
N = 52 # weekly for a year

# Initial state
X0 = c(16000,10,0,0)

# 
N.STATES = 3
N.RXNS   = 2

# Parameters
ItoRrate = c(.9,1.1)
R        = c(1.5,2.0,2.5)

probs    = c(.1,.5)
parms    = expand.grid(ItoRrate=ItoRrate,R=R,StoIprob=probs, ItoRprob=probs)

parms$StoIrate = apply(parms[,1:2],1,prod)

thetas = cbind(parms$StoIrate, parms$ItoRrate, 0, 0)
probs  = cbind(parms$StoIprob, parms$ItoRprob, 0, 0)


##################### Generate simulations ########################
sims = list()
for (i in 1:nrow(thetas)) sims[[i]] = SIRDsim(X0,thetas[i,],probs[i,],N)

prior = list()
prior$theta = list()
prior$theta$a = c(200, 100,  0,  0)
prior$theta$b = c(100, 100, 10, 10)
prior$p$a     = c(1, 1, 1, 1)
prior$p$b     = c(1, 1, 1, 1)

# Find the end of the outbreak
# defined as 4 consecutive weeks of zero S->I observed
n = numeric(length(sims))
for (i in 1:nrow(thetas)) {
  zeros = which(sims[[i]]$y$StoI==0)
  ymax  = which.max(sims[[i]]$y$StoI)
  zeros = zeros[zeros>ymax]
  diff3 = diff(zeros,3)
  n[i] = min(zeros[match(3,diff3)]+3,N)
}

save(sims,thetas,probs,n,prior, file="SIRDsims.RData")

q("no")

