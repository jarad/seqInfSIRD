
predictive.likelihood = function(y, sys, part, log=T, engine="R") 
{
    engine = pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
#        ph = part$p * hazard.part
#        prob = 
#        sum(dnbinom(
    },
    {
    })
}

