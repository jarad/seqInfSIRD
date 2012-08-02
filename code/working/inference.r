# Functions for discrete-time compartment models

dyn.load("inference.so")


# Model creation functions
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

# Utility functions

engine.error = function() stop("'engine' must be 'C' or 'R'")



# Model functions

hazard.part = function(sys, engine="R") 
{ 
    if (engine=="R") 
    {
      h = rep(1,sys$r)
      for (i in 1:sys$r) 
      {
        for (j in 1:sys$s) 
        {
          h[i] = h[i]*choose(sys$X[j], sys$Pre[i,j]) 
        }
      }          
      return(h)
    } else if (engine=="C") 
    {
      out = .C("hazard_part", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), h=integer(sys$r))
      return(out$h)
    } else 
    {
        engine.error()
    }
}

hazard = function(sys, engine="R") 
{ 
    if (engine=="R") 
    {
      hp = hazard.part(sys)
      h = numeric(sys$r)
      for (i in 1:sys$r) h[i] = hp[i]*sys$theta[i]
      return(list(h=h,hp=hp))
    } else if (engine=="C") 
    {
      out = .C("hazard", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), as.double(sys$theta), 
               hp=integer(sys$r), h=double(sys$r))
      return(list(h=out$h,hp=out$hp))
    } else 
    {
        engine.error()
    }
}

sim.poisson = function(sys, engine="R")
{
    if (engine=="R")
    {
        return(rpois(sys$r,sys$h))
    } else if (engine=="C")
    {
        out = .C("sim_poisson",
                 as.integer(sys$r), as.double(sys$h), r=integer(sys$r))
        return(out$r)
    } else 
    {
        engine.error()
    }
}

update.species = function(sys, nRxns, engine="R") 
{
    if (engine=="R")
    {
        return(as.numeric(sys$X+sys$stoich%*%nRxns))
    } else if (engine=="C")
    {
        out = .C("update_species",
                 as.integer(sys$s), as.integer(sys$r), X=as.integer(sys$X),
                 as.integer(sys$stoich), as.integer(nRxns))
        return(out$X)
    } else 
    {
        engine.error()
    }
}

sim.one.step = function(sys, engine="R")
{
    if (engine=="R")
    {
        r = sim.poisson(sys)
        whileCount = 0
        X = update.species(sys,r,engine="R")
        while ( any(X<0) ) 
        {
            r = sim.poisson(sys)
            X = update.species(sys,r,engine="R")
            whileCount = whileCount+1
            if (whileCount>1000) stop("R:sim.one.step: Too many unsuccessful simulation iterations.")
        }
        return(list(X=X,dX=r))    
    } else if (engine=="C")
    {
        out = .C("sim_one_step", 
                 as.integer(sys$s), as.integer(sys$r), X=as.integer(sys$X),
                 as.integer(sys$stoich), r=integer(sys$r), as.double(sys$h))
        return(list(X=out$X,dX=out$r))
    } else 
    {
        engine.error()
    }
}



