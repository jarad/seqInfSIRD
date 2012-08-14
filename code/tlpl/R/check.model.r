check.model = function(sys) {
    stopifnot(sys$s == ncol(  sys$Pre),
              sys$s == ncol(  sys$Post),
              sys$s == nrow(  sys$stoich),
              sys$s == length(sys$X),
              sys$r == nrow(  sys$Pre),
              sys$r == nrow(  sys$Post),
              sys$r == ncol(  sys$stoich),
              sys$r == length(sys$theta),
              all(sys$theta>=0),
              all(sys$X    >=0),
              all(sys$Pre  >=0),
              all(sys$Post >=0))

    if (!all.equal(sys$stoich,t(sys$Post-sys$Pre))) warning("sys$stoich!=t(sys$Post-sys$Pre))")
}


