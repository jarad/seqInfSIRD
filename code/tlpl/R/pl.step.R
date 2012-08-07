
# A single step of PL in tau-leaping model

pl.step = function(y,tau,mod,part,while.max=1000) 
{
    check.model(mod)

    out = .C("discrete_all_particle_update",
             as.integer(mod$s), as.integer(mod$r), as.integer(t(mod$Pre)), as.integer(mod$stoich),
             as.integer(y), as.double(tau), as.integer(part$n), as.integer(while.max),
             X=as.integer(t(part$X)), hyper=as.double(t(part$hyper)))
    return(list(X    = matrix(out$X, part$n, sys$s, byrow=T), 
                hyper= matrix(out$hyper, part$n, 4*sys$r, byrow=T)))
} 

