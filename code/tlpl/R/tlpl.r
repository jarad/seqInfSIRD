# A convenience function to make the prior
tlpl.prior = function(X, p.a, p.b, r.a, r.b, nr) 
{
    if (length(p.a)==1) p.a = rep(p.a, nr)
    if (length(p.b)==1) p.b = rep(p.b, nr)
    if (length(r.a)==1) r.a = rep(r.a, nr)
    if (length(r.b)==1) r.b = rep(r.b, nr)

    stopifnot(length(p.a)==nr, length(p.b)==nr, length(r.a)==nr, length(r.b)==nr)

    prior = list()
    prior$X = sckm$X
    prior$prob = list()
    prior$prob$a = p.a
    prior$prob$b = p.b
    prior$rate = list()
    prior$rate$a = r.a
    prior$rate$b = r.b

    return(prior)
}





tlpl = function(data, sckm, swarm=NULL, prior=NULL, n.particles=NULL, 
                mult=rep(1, sckm$r), engine="R", verbose=F, ...)
{
    # mult is a constant multiplier for each reaction

    nr = sckm$r
    ns = sckm$s 

    # Check data, sckm, swarm
    n = nrow(data$y)
    if (length(data$tau)==1) data$tau = rep(data$tau,n)
    stopifnot(length(data$tau)==n)
    stopifnot(ncol(data$y)==nr) # Only observations on transitions are currently implemented

    check.system(sckm)
    if (!is.null(swarm)) check.swarm(swarm)

    stopifnot(length(mult)==sckm$r)
    lmult = log(mult)

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
            prior = tlpl.prior(sckm$X, 1, 1, 1, 1, nr)
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
    hp = matrix(NA, nr, np)

    # Run through all data points
    for (i in 1:n) 
    {  
        if (verbose) cat(paste("Time point ",i,", ",round(i/n*100), "% completed.\n", sep=''))
        y = data$y[i,]
        tau = data$tau[i]

        # Sample observation probability
        swarm$p = matrix(rbeta(nr*np, swarm$hyper$prob$a, swarm$hyper$prob$b), nr, np)

        # Calculate all particle weights
        for (j in 1:n.particles) 
        {          
            for (k in 1:nr) hp[k,j] = exp(sum(lchoose(swarm$X[,j], sckm$Pre[k,]))+lmult[k])

            ph = swarm$p[,j] * hp[,j] * tau
            prob = ph/(swarm$hyper$rate$b[,j]+ph) 
            nz = which(ph>0) # otherwise NaNs produced
            w[j] = sum(dnbinom(y[nz], swarm$hyper$rate$a[nz,j], 1-prob[nz], log=T))

            # If particle outbreak is over but data indicates continuing outbreak,
            # particle weight becomes 0 ( log(weight)=-Inf )
            if (any(ph==0)) { if (any(y[ph==0]!=0)) w[j] = -Inf }
        }

        # Resample particles
        w = renormalize.weights(w, log=T)
        rs = resample(w,...)$indices

        for (j in 1:np)
        {               
            if (verbose && (j%%100)==0) 
                cat(paste("  Particle ",j,", ",round(j/np*100), "% completed.\n", sep=''))

            any.negative = T
            while (any.negative) {
                # To ensure a new particle is resampled
                # Clearly reasonable for multinomial resampling, but what about the rest? 
                kk = resample(w,1, method="multinomial")$indices # new particle id

                # Calculate mean for unobserved transitions
                lambda = rgamma(nr, swarm$hyper$rate$a[,kk], swarm$hyper$rate$b[,kk])
                mn = (1-swarm$p[,kk])* lambda * tau * hp[,kk]

                # Sample transitions and update state
                z = rpois(nr, mn) # unobserved transitions
                n.rxns = y + z    # total transitions
                newswarm$X[,j] = swarm$X[,kk] + sckm$stoich %*% n.rxns

                # Check to see if any state is negative 
                any.negative = any(newswarm$X[,j]<0)

                if (any.negative && verbose) {
                    cat(paste("Particle",kk,"failed.\n"))
                    cat("Probabilities: ")
                    for (k in 1:nr) cat(paste(swarm$p[k,kk]," "))
                    cat("\nRates: ")
                    for (k in 1:nr) cat(paste(lambda[k]," "))
                    cat("\nStates: ")
                    for (k in 1:ns) cat(paste(swarm$X[k,kk]," "))
                    cat("\nData: ")
                    for (k in 1:nr) cat(paste(y[k]," "))
                    cat("\nHazard parts: ")
                    for (k in 1:nr) cat(paste(hp[k,kk]," "))
                    cat(paste("\nWeight:",w[kk],"\n"))
                }
            }

            # Update sufficient statistics
            newswarm$hyper$prob$a[,j] = swarm$hyper$prob$a[,kk] + y
            newswarm$hyper$prob$b[,j] = swarm$hyper$prob$b[,kk] + z
            newswarm$hyper$rate$a[,j] = swarm$hyper$rate$a[,kk] + n.rxns
            newswarm$hyper$rate$b[,j] = swarm$hyper$rate$b[,kk] + hp[,kk] * tau
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




