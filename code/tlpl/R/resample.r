

is.increasing = function(v, engine="R") 
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


residual.resample = function(weights, num.samples, residual.resampling.function=c("stratified","multinomial","systematic"))
{
    rrf = pmatch(residual.resampling.function, c("stratified","multinomial","systematic"))
    if (is.na(rrf)) stop("No matching residual resampling function.")

    out = .C("residual_resample_wrap",
             as.integer(length(weights)), as.double(weights), 
             as.integer(num.samples), id=integer(num.samples), as.integer(rrf))

    return(out$id)
}

