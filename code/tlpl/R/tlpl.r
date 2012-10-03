library(smcUtils)


# Calculates the predictive likelihood for each particle
calc.pred.like = function(y, tau, sys, part, log=T, engine="R") 
{
    stopifnot(length(y) == sys$r)
    engine = pmatch(engine, c("R","C"))
    check.system(sys)

    switch(engine,
    {
        # R implementation
        ph = part$p * hazard.part(sys) * tau
        prob = ph/(part$hyper$rate$a+ph) 
        nz = which(ph>0)
        logPL = sum(dnbinom(y[nz],part$hyper$rate$b[nz],prob[nz],log=T))
        return(ifelse(log,logPL,exp(logPL)))
    },
    {
        # C implementation
        out = .C("calc_log_pred_like_R",
                 as.integer(y), as.double(tau), 
                 as.integer(sys$s), as.integer(sys$r), as.integer(t(sys$Pre)), as.integer(t(sys$Post)),
                 as.integer(part$X), 
                 as.double(part$hyper$prob$a), as.double(part$hyper$prob$b), 
                 as.double(part$hyper$rate$a), as.double(part$hyper$rate$b), 
                 as.double(part$p), as.double(part$rate),
                 logPL=double(1))

        return(ifelse(log, out$logPL, exp(out$logPL)))
    },
    {
        stop(paste("No engine matching ",engine,".",sep=""))
    })
}









tlpl = function(data, sckm, swarm=NULL, prior=NULL, n.particles=NULL, engine="R")
{
    nr = sckm$r
 
    # Check data, sckm, swarm
    n = nrow(data$y)
    stopifnot(length(data$tau)==n)
    stopifnot(ncol(data$y)==nr) # Only observations on transitions are currently implemented

    check.system(sckm)
    if (!is.null(swarm)) check.swarm(swarm)

    # Create swarm
    if (is.null(swarm)) 
    {
        if (is.null(prior)) stop("Both prior and swarm cannot be NULL.")
        np = n.particles

        swarm = list()
        swarm$n.particles = np
        swarm$weights = rep(1/np, np)
        swarm$normalized = TRUE
        swarm$log.weights = FALSE
        swarm$x = matrix( ,s,np)
        swarm$hyper = list()
        swarm$hyper$prob = list()
        swarm$hyper$prob$a = matrix(prior$prob$a, r, np)
        swarm$hyper$prob$b = matrix(prior$prob$b, r, np)
        swarm$hyper$rate = list()
        swarm$hyper$rate$a = matrix(prior$rate$a, r, np)
        swarm$hyper$rate$b = matrix(prior$rate$b, r, np)

        
    } else { # !is.null(swarm)
        np = swarm$n.particles
    }

    # Create output 
    out = list()
    out$x = array(NA, dim=c(n+1,s,np))
    out$hyper = list()
    out$hyper$prob = list()
    out$hyper$prob$a = array(NA, dim=c(n+1,r,np))
    out$hyper$prob$b = array(NA, dim=c(n+1,r,np))
    out$hyper$rate = list()
    out$hyper$rate$a = array(NA, dim=c(n+1,r,np))
    out$hyper$rate$b = array(NA, dim=c(n+1,r,np))

    # Fill output with initial values
    out$x[1,] = swarm$x
    out$hyper$prob$a[1,,] = swarm$hyper$prob$a
    out$hyper$prob$b[1,,] = swarm$hyper$prob$b
    out$hyper$rate$a[1,,] = swarm$hyper$rate$a
    out$hyper$rate$b[1,,] = swarm$hyper$rate$b

    # Create variables used throughout
    part = list()
    part$hyper = list()
    part$hyper$prob = list()
    part$hyper$rate = list()
    w = rep(NA, np)
    newswarm = swarm


    # Run through all data points
    for (i in 1:n) 
    {  
        # Sample observation probability
        swarm$p = matrix(rbeta(r*np, swarm$hyper$prob$a, swarm$hyper$prob$b), r, np)

        # Calculate all particle weights
        for (j in 1:n.particles) 
        {           
            # Extract a particle
            part$p = swarm$p[,j]
            part$hyper$prob$a = swarm$hyper$prob$a[,j]
            part$hyper$prob$b = swarm$hyper$prob$b[,j]            
            part$hyper$rate$a = swarm$hyper$rate$a[,j]
            part$hyper$rate$b = swarm$hyper$rate$b[,j]
            
            # Calculate particle's weight
            w[j] = calc.pred.like(data$y[i,], data$tau[i], sckm, part)
        }

        # Resample particles
        w = renormalize.weights(w, log=swarm$log.weights)
        k = resample(w)$indices
        for (j in 1:n.particles)
        {
            part$p = swarm$p[,k[j]]
            part$hyper$prob$a = swarm$hyper$prob$a[,k[j]]
            part$hyper$prob$b = swarm$hyper$prob$b[,k[j]]            
            part$hyper$rate$a = swarm$hyper$rate$a[,k[j]]
            part$hyper$rate$b = swarm$hyper$rate$b[,k[j]]

            # Sample trajectory
            newswarm$x[,j] = tau.leap(sckm, tau=tau[i])[2,]

            # Update states
            # Update sufficient statistics

        }
        swarm = newswarm

        # Fill output with current values
        out$x[i+1,] = swarm$x
        out$hyper$prob$a[i+1,,] = swarm$hyper$prob$a
        out$hyper$prob$b[i+1,,] = swarm$hyper$prob$b
        out$hyper$rate$a[i+1,,] = swarm$hyper$rate$a
        out$hyper$rate$b[i+1,,] = swarm$hyper$rate$b
    }

    return(out)
}




