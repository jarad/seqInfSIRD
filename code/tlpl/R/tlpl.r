tlpl = function(data, sckm, swarm=NULL, prior=NULL, n.particles=NULL, engine="R",...)
{
    nr = sckm$r
    ns = sckm$s 

    # Check data, sckm, swarm
    n = nrow(data$y)
    if (length(data$tau)==1) data$tau = rep(data$tau,n)
    stopifnot(length(data$tau)==n)
    stopifnot(ncol(data$y)==nr) # Only observations on transitions are currently implemented

    check.system(sckm)
    if (!is.null(swarm)) check.swarm(swarm)

    # Create swarm
    if (is.null(swarm)) 
    {
        if (is.null(n.particles)) 
        {
            # Determine number of particles based on number of reactions/species
            # and number of time points
            # n.particles = 2 set for testing purposes
            n.particles = 2
        }

        if (is.null(prior)) 
        {
            prior = list()
            prior$X = rep(1, ns)
            prior$prob = list()
            prior$prob$a = rep(1,nr)
            prior$prob$b = rep(1,nr)
            prior$rate = list()
            prior$rate$a = rep(1,nr)
            prior$rate$b = rep(1,nr)
        }
        np = n.particles

        swarm = list()
        swarm$n.particles = np
        swarm$weights = rep(1/np, np)
        swarm$normalized = TRUE
        swarm$log.weights = FALSE
        swarm$X = matrix(prior$X, ns, np)
        swarm$hyper = list()
        swarm$hyper$prob = list()
        swarm$hyper$prob$a = matrix(prior$prob$a, nr, np)
        swarm$hyper$prob$b = matrix(prior$prob$b, nr, np)
        swarm$hyper$rate = list()
        swarm$hyper$rate$a = matrix(prior$rate$a, nr, np)
        swarm$hyper$rate$b = matrix(prior$rate$b, nr, np)

        
    } else { # !is.null(swarm)
        np = swarm$n.particles
    }

    # Create output 
    out = list()
    out$X = array(NA, dim=c(n+1,ns,np))
    out$hyper = list()
    out$hyper$prob = list()
    out$hyper$prob$a = array(NA, dim=c(n+1,nr,np))
    out$hyper$prob$b = array(NA, dim=c(n+1,nr,np))
    out$hyper$rate = list()
    out$hyper$rate$a = array(NA, dim=c(n+1,nr,np))
    out$hyper$rate$b = array(NA, dim=c(n+1,nr,np))

    # Fill output with initial values
    out$X[1,,] = swarm$X
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
        y = data$y[i,]
        tau = data$tau[i]

        # Sample observation probability
        swarm$p = matrix(rbeta(nr*np, swarm$hyper$prob$a, swarm$hyper$prob$b), nr, np)

        # Calculate all particle weights
        for (j in 1:n.particles) 
        {           
            for (k in 1:nr) hp[k] = exp(sum(lchoose(swarm$X[,j], sckm$Pre[k,])))

            ph = swarm$p[,j] * hp * tau
            prob = ph/(swarm$hyper$rate$a[,j]+ph) 
            nz = which(ph>0)
            w[j] = sum(dnbinom(y[nz],swarm$hyper$rate$b[nz,j],prob[nz],log=T))
        }

        # Resample particles
        w = renormalize.weights(w, log=T)
        rs = resample(w,...)$indices

        for (j in 1:n.particles)
        {
            kk = rs[j] # new particle id

            # Calculate mean for unobserved transitions
            lambda = rgamma(nr, swarm$hyper$rate$a[,kk], swarm$hyper$rate$b[,kk])
            for (k in 1:nr) hp[k] = exp(sum(lchoose(swarm$X[,kk], sckm$Pre[k,])))
            mn = (1-swarm$p[,kk])* lambda * tau * hp

            # Sample transitions and update state
            z = rpois(nr, mn) # unobserved transitions
            n.rxns = y + z    # total transitions
            newswarm$X[,j] = swarm$X[,kk] + sckm$stoich %*% n.rxns 

            # Update sufficient statistics
            newswarm$hyper$prob$a[,j] = swarm$hyper$prob$a[,kk] + y
            newswarm$hyper$prob$b[,j] = swarm$hyper$prob$b[,kk] + z
            newswarm$hyper$rate$a[,j] = swarm$hyper$rate$a[,kk] + n.rxns
            newswarm$hyper$rate$b[,j] = swarm$hyper$rate$b[,kk] + hp * tau
        }
        swarm = newswarm

        # Fill output with current values
        out$X[i+1,,] = swarm$X
        out$hyper$prob$a[i+1,,] = swarm$hyper$prob$a
        out$hyper$prob$b[i+1,,] = swarm$hyper$prob$b
        out$hyper$rate$a[i+1,,] = swarm$hyper$rate$a
        out$hyper$rate$b[i+1,,] = swarm$hyper$rate$b
    }

    return(out)
}




