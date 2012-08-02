# Functions for discrete-time compartment models

dyn.load("inference.so")

engine.error = function() stop("'engine' must be 'C' or 'R'")

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
      h = hazard.part(sys)
      for (i in 1:sys$r) h[i] = h[i]*sys$theta[i]
      return(h)
    } else if (engine=="C") 
    {
      out = .C("hazard", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), as.double(sys$theta), h=double(sys$r))
      return(out$h)
    } else 
    {
        engine.error()
    }
}


sim.poisson = function(sys, engine="R")
{
    h = hazard(sys)
    if (engine=="R")
    {
        return(rpois(sys$r,h))
    } else if (engine=="C")
    {
        out = .C("sim_poisson",
                 as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
                 as.integer(t(sys$Pre)), as.double(h), r=integer(sys$r))
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
                 as.integer(t(sys$Pre)), as.integer(sys$stoich), as.double(sys$theta),
                 r=integer(sys$r))
        return(list(X=out$X,dX=out$r))
    } else 
    {
        engine.error()
    }
}

