# Script to test C code

source("inference.r")

random.system = function(s=rpois(1,5)+1, r=rpois(1,5)+1) {
  d = s*r
  Pre  = matrix(rpois(d,.5), r, s)
  Post = matrix(rpois(d,.5), r, s)
  stoich = t(Post-Pre)
  X = rpois(s,5)
  theta = rgamma(r,1)
  sys = list(s=s,r=r,Pre=Pre,Post=Post,stoich=stoich,X=X,theta=theta)
  sys$h = hazard(sys)$h
  
  return(sys)
}


random.sir = function() {
  s = 3
  r = 2
  Pre  = rbind(c(1,1,0),c(0,1,0))
  Post = rbind(c(0,2,0),c(0,0,1))
  stoich = t(Post-Pre)
  X = c(rpois(1,1000),rpois(1,5),0)
  theta = rgamma(r,10,10)/c(sum(X),1) # S->I scaled by N
  sys = list(s=s,r=r,Pre=Pre,Post=Post,stoich=stoich,X=X,theta=theta)
  sys$h = hazard(sys)$h
  
  return(sys)
}


for (i in 1:100) 
{
  sys = random.system()

  stopifnot(all.equal(hazard.part(sys, engine="R"), hazard.part(sys, engine="C")))

  stopifnot(all.equal(hazard(sys, engine="R"), hazard(sys, engine="C")))

  sys$h = hazard(sys)$h
  a = .Random.seed; b = sim.poisson(sys, engine="R")
  .Random.seed = a; c = sim.poisson(sys, engine="C")
  stopifnot(all.equal(b,c))

  r = sim.poisson(sys)
  stopifnot(all.equal(update.species(sys,r, engine="R"), update.species(sys,r,engine="C"))) 
}


for (i in 1:100) 
{
  sir = random.sir()

  a = .Random.seed; r = sim.one.step(sir, engine="R")
  .Random.seed = a; c = sim.one.step(sir, engine="C")
  stopifnot(all.equal(r,c))
}



