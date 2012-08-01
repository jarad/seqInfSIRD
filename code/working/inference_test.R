# Script to test C code

source("inference.r")

random.system = function(s=rpois(1,5)+1, r=rpois(1,5)+1) {
  d = s*r
  Pre  = matrix(rpois(d,.5), r, s)
  Post = matrix(rpois(d,.5), r, s)
  stoich = t(Post-Pre)
  
  return(list(s=s,r=r,Pre=Pre,Post=Post,stoich=stoich))
}



for (i in 1:10) 
{
  sys = random.system()

  sys$X = rpois(sys$s,5)
  stopifnot(all.equal(hazard.part(sys, engine="R"), hazard.part(sys, engine="C")))

  sys$theta=rgamma(sys$r,1)
  stopifnot(all.equal(hazard(sys, engine="R"), hazard(sys, engine="C")))
}

