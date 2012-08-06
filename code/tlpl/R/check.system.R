check.system = function(sys) {
    stopifnot(sys$s == ncol(  sys$Pre),
              sys$s == ncol(  sys$Post),
              sys$s == nrow(  sys$stoich),
              sys$s == length(sys$X),
              sys$r == nrow(  sys$Pre),
              sys$r == nrow(  sys$Post),
              sys$r == ncol(  sys$stoich),
              sys$r == length(sys$theta))
}


