# Script to test C code

source("inference.r")

random.system = function(s=rpois(1,5)+1, r=rpois(1,5)+1) {
  d = s*r
  Pre  = matrix(rbinom(d,1,.2), r, s)
  Post = matrix(rbinom(d,1,.2), r, s)
  stoich = t(Post-Pre)
  
  return(list(s=s,r=r,Pre=Pre,Post=Post,stoich=stoich))
}

sys = random.system()


sys$X = rpois(sys$s,5)
stopifnot(all.equal(hazard.part(sys, engine="R"), hazard.part(sys, engine="C")))

sys$theta=rgamma(sys$r,1)
stopifnot(all.equal(hazard(sys, engine="R"), hazard(sys, engine="C")))


