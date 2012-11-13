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
        hyper = list()
        hyper$prob = list(a=rep(1,r), b=rep(1,r))
        hyper$rate = list(a=rep(1,r), b=rep(1,r))
        foundi=TRUE        
    })

    # Used in place of a default switch case
    if (!foundi) stop("Test system doesn't exist.")

    return(list(X=X,p=p,hyper=hyper))    
}


check.swarm = function(swarm)
{
    stopifnot(!is.null(swarm),
              !is.null(swarm$n.particles),
              !is.null(swarm$weights),
              !is.null(swarm$normalized),
              !is.null(swarm$log.weights),
              !is.null(swarm$X),
              !is.null(swarm$hyper),
              !is.null(swarm$hyper$prob),
              !is.null(swarm$hyper$prob$a),
              !is.null(swarm$hyper$prob$b),
              !is.null(swarm$hyper$rate),
              !is.null(swarm$hyper$rate$a),
              !is.null(swarm$hyper$rate$b))

    np = swarm$n.particles
    nr = nrow(swarm$hyper$prob$a)

    stopifnot(np == length(swarm$weights),
              np == ncol(swarm$X),
              np == ncol(swarm$hyper$prob$a),
              np == ncol(swarm$hyper$prob$b),
              np == ncol(swarm$hyper$rate$a),
              np == ncol(swarm$hyper$rate$b),
              nr == nrow(swarm$hyper$prob$a),
              nr == nrow(swarm$hyper$prob$b),
              nr == nrow(swarm$hyper$rate$a),
              nr == nrow(swarm$hyper$rate$b))
}



