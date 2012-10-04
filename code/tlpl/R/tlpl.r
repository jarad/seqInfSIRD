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
    hp = numeric(nr)

    # Run through all data points
    for (i in 1:n) 
    {  
        y = data$[i,]
        tau = data$tau[i]

        # Sample observation probability
        swarm$p = matrix(rbeta(r*np, swarm$hyper$prob$a, swarm$hyper$prob$b), r, np)

        # Calculate all particle weights
        for (j in 1:n.particles) 
        {           
            for (k in 1:nr) hp[k] = exp(sum(lchoose(swarm$x[,j], sckm$Pre[k,])))

            ph = swarm$p[,j] * hp * tau
            prob = ph/(swarm$hyper$rate$a[,j]+ph) 
            nz = which(ph>0)
            w[j] = sum(dnbinom(y[nz],swarm$hyper$rate$b[nz,j],prob[nz],log=T))
        }

        # Resample particles
        w = renormalize.weights(w, log=T)
        rs = resample(w)$indices

        for (j in 1:n.particles)
        {
            kk = rs[j] # new particle id

            # Calculate mean for unknown transitions
            lambda = rgamma(nr, swarm$hyper$rate$a[,kk], swarm$hyper$rate$b[,kk])
            for (k in 1:nr) hp[k] = exp(sum(lchoose(swarm$x[,kk], sckm$Pre[k,])))
            mn = (1-swarm$p[,kk])* lambda * tau * hp

            # Sample transitions and update state
            z = rpois(nr, mn) # unobserved transitions
            n.rxns = y + z    # total transitions
            newswarm$x[,j] = swarm$x[,kk] + sckm$stoich %*% n.rxns 

            # Update sufficient statistics
            newswarm$hyper$prob$a[,j] = swarm$hyper$prob$a[,kk] + y
            newswarm$hyper$prob$b[,j] = swarm$hyper$prob$b[,kk] + z
            newswarm$hyper$rate$a[,j] = swarm$hyper$rate$a[,kk] + n.rxns
            newswarm$hyper$rate$b[,j] = swarm$hyper$rate$b[,kk] + hp * tau
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




