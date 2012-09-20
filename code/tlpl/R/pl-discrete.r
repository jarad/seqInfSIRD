
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

