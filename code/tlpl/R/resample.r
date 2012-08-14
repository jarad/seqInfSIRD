
###############################################################
# Utility functions
###############################################################


is.increasing = function(v, engine="C") 
{
    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        return(!is.unsorted(v))
    },
    {
        # C implementation
        out = .C("is_increasing_wrap", 
                 as.integer(length(v)), 
                 as.double(v), 
                 increasing=integer(1))
        return(as.logical(out$increasing))
    })
}


cusum = function(v, engine="C")
{
    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        return(cumsum(v))
    },
    {
        # C implementation
        out = .C("cumulative_sum_wrap", 
                 as.integer(length(v)), 
                 cusum=as.double(v))
        return(out$cusum)
    })
}


rep2id = function(rep, engine="C") 
{
    engine=pmatch(engine, c("R","C"))
    sum = sum(rep)

    switch(engine,
    {
        # R implementation
        id = integer(sum)
        current.index = integer(1)
        for (i in 1:length(rep)) 
        {
            if (rep[i] != 0) 
            {
                id[current.index + 1:rep[i]] = i
                current.index = current.index + rep[i]
            }
        }
        return(id)
    },
    {
        # C implementation
        out = .C("rep2id_wrap", 
                 as.integer(rep),
                 as.integer(sum), 
                 id=integer(sum))
        return(out$id)
    })  
}


inverse.cdf.weights = function(weights, uniforms, engine="C")
{
    check.weights(weights)
    check.weights(uniforms)
    
    engine=pmatch(engine, c("R","C"))
    if (is.unsorted(uniforms)) uniforms = sort(uniforms)
    n.samples = length(uniforms)

    switch(engine,
    {
        # R implementation
        ids       = integer(n.samples)
        cusum     = cumsum(weights)
        index     = 1
        for (i in 1:n.samples) 
        {
            found = FALSE
            while (!found) 
            {
                if (uniforms[i] > cusum[index]) 
                {
                    index = index + 1
                }
                else 
                {
                    found = TRUE
                }
            }
            ids[i] = index
        }
        return(ids)
    },
    {
        # C implementation
        out = .C("inverse_cdf_weights_wrap", 
                 as.integer(length(weights)),
                 as.double(weights),
                 as.integer(n.samples),
                 as.double(uniforms), 
                 id=integer(n.samples))
        return(out$id)
    })  
}


renormalize = function(weights, log=F, engine="C")
{   
    if (!log && any(weights<0)) stop("log=F but negative weights exist.")

    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        if (log) weights = exp(weights - max(weights))
        return(weights/sum(weights))
    },
    {
        # C implementation
        n = length(weights)
        out = .C("renormalize_wrap", 
                 as.integer(n),
                 as.integer(log),  
                 weights=as.double(weights))
        return(out$weights)
    })   
}



###############################################################
# Effective sample size functions
###############################################################


ess = function(weights, engine="C") 
{
    check.weights(weights, log=F, normalized=T)
    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        return(1/sum(weights^2))
    },
    {
        # C implementation
        out = .C("ess_wrap", 
                 as.integer(length(weights)),
                 as.double(weights), 
                 ess = double(1))
        return(out$ess)
    })
}


cov2 = function(weights, engine="C") 
{
    check.weights(weights, log=F, normalized=T)
    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        return(var(weights)/mean(weights)^2)
    },
    {
        # C implementation
        out = .C("cov2_wrap", 
                 as.integer(length(weights)),
                 as.double(weights), 
                 cov2 = double(1))
        return(out$cov2) 
    })  
}



entropy = function(weights, engine="C") 
{
    check.weights(weights, log=F, normalized=T)
    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        return(-sum(weights * log2(weights + .Machine$double.eps)))
    },
    {
        # C implementation
        out = .C("entropy_wrap", 
                 as.integer(length(weights)),
                 as.double(weights), 
                 entropy = double(1))
        return(out$entropy)
    })
}




###############################################################
# Resampling functions
###############################################################

stratified.resample = function(weights, n.samples=length(weights), engine="C")
{
    check.weights(weights, log=F, normalized=T)
    stopifnot(n.samples>0)

    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        lbs = seq(0, by=1/n.samples, length=n.samples)
        ubs = lbs+1/n.samples
        u = runif(n.samples, lbs, ubs)
        return(inverse.cdf.weights(weights,u,engine="R"))
    },
    {
        # C implementation
        out = .C("stratified_resample_wrap", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(n.samples),
                 id = integer(n.samples))
        return(out$id)
    })
}


multinomial.resample = function(weights, n.samples=length(weights), engine="C")
{
    check.weights(weights, log=F, normalized=T)
    stopifnot(n.samples>0)

    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        # R implementation
        # sample does not perform the same as inverse.cdf method
        # so to be consistent with C, use inverse.cdf
        #return(sort(sample(n, n.samples, replace = TRUE, 
        #                   prob = weights)))
        u = runif(n.samples)
        return(inverse.cdf.weights(weights,u,engine="R"))
    },
    {
        # C implementation
        out = .C("multinomial_resample_wrap", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(n.samples),
                 id = integer(n.samples))
        return(out$id)
    })
}





residual.resample = function(weights, n.samples, 
   rrf=c("stratified","multinomial","systematic"), 
   engine="C")
{
    rrf = pmatch(rrf, 
                 c("stratified","multinomial","systematic"))
    if (is.na(rrf)) stop("No matching residual resampling function.")

    out = .C("residual_resample_wrap",
             as.integer(length(weights)), as.double(weights), 
             as.integer(n.samples), id=integer(n.samples), as.integer(rrf))

    return(out$id)
}

