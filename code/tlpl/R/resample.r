

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
    stopifnot(length(weights)>0, length(uniforms)>0,
              all(weights>0), all(uniforms>0))
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


residual.resample = function(weights, n.samples, residual.resampling.function=c("stratified","multinomial","systematic"))
{
    rrf = pmatch(residual.resampling.function, c("stratified","multinomial","systematic"))
    if (is.na(rrf)) stop("No matching residual resampling function.")

    out = .C("residual_resample_wrap",
             as.integer(length(weights)), as.double(weights), 
             as.integer(n.samples), id=integer(n.samples), as.integer(rrf))

    return(out$id)
}

