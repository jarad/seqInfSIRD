# Script to test C code

source("inference.r")

random.system = function(s=rpois(1,5)+1, r=rpois(1,5)+1) {
  d = s*r
  Pre  = matrix(rpois(d,.5), r, s)
  Post = matrix(rpois(d,.5), r, s)
  stoich = t(Post-Pre)
  X = rpois(s,5)
  theta = rgamma(r,1)
  
  return(list(s=s,r=r,Pre=Pre,Post=Post,stoich=stoich,X=X,theta=theta))
}



for (i in 1:100) 
{
  sys = random.system()

  stopifnot(all.equal(hazard.part(sys, engine="R"), hazard.part(sys, engine="C")))

  stopifnot(all.equal(hazard(sys, engine="R"), hazard(sys, engine="C")))

  a = .Random.seed; b = sim.poisson(sys, engine="R")
  .Random.seed = a; c = sim.poisson(sys, engine="C")
  stopifnot(all.equal(b,c))

  r = rpois(sys$r, hazard(sys))
  stopifnot(all.equal(update.species(sys,r, engine="R"), update.species(sys,r,engine="C"))) 

  a = .Random.seed; b = sim.one.step(sys, engine="R")
  .Random.seed = a; c = sim.one.step(sys, engine="C")
  stopifnot(all.equal(b,c))
}

