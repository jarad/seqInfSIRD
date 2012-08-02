# Script to test C code

source("inference.r")

n.reps = 100
  


for (i in 1:n.reps) 
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


for (i in 1:n.reps) 
{
  sir = random.sir()

  a = .Random.seed; r = sim.one.step(sir, engine="R")
  .Random.seed = a; c = sim.one.step(sir, engine="C")
  stopifnot(all.equal(r,c))

  y = rbinom(sir$r,r$dX,sir$p)
  a = .Random.seed; r = sim.one.step(sir, y, engine="R")
  .Random.seed = a; c = sim.one.step(sir, y, engine="C")
  stopifnot(all.equal(r,c))

  dX = sim.one.step(sir, y, engine="R")$dX
  stopifnot(all.equal(inf.one.step(sir,y,dX, engine="R"), inf.one.step(sir,y,dX, engine="C")))
}



