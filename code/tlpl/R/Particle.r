# Functions for creating particles

# This particle complements the system created by test.system
test.particle = function(i=1)
{
    s = 3
    r = 8
    foundi = FALSE
    switch(i,
    {
        X = c(0,1,2)
        p = rep(.5,r)
        hyper = matrix(1:(s*r),2,3)
        foundi=TRUE        
    })

    # Used in place of a default switch case
    if (!foundi) stop("Test system doesn't exist.")

    return(list(X=X,p=p,hyper=hyper))    
}


