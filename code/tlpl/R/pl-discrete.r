
predictive.likelihood = function(y, sys, part, log=T, engine="R") 
{
    engine = pmatch(engine, c("R","C"))
    check.system(sys)

    switch(engine,
    {
        # R implementation
        stop("Not implemented.")
        ph = part$p * hazard.part(sys)
        #prob = ph/(part$hyper+ph) 
        logPL = sum(dnbinom(y,gamma,prob,log=T))
        return(ifelse(log,logPL,exp(logpL)))
    },
    {
        # C implementation
        out = .C("calculate_log_predictive_likelihood_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(t(sys$Pre)), 
                 as.integer(y), as.double(tau), 
                 as.integer(part$X), as.double(part$p), as.double(part$hyper),
                 logPL=double(1))
        return(ifelse(log, out$logPL, exp(out$logPL)))
    },
    {
        stop(paste("No engine matching ",engine,".",sep=""))
    })
}

