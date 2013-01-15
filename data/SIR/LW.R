# Performs LW for preset tuning parameter
library(plyr)

load("sims.RData")
source("filter-settings.R")
source("liu_west.r")

lw = llply(lapply(sims, function(x) return(x$y)), liu_west,
           sckm=sckm, 
           n.particles=n.particles,
           delta=delta, 
           prior=prior,
           method=resampling.function)

save(lw, file=paste("LW", round(delta*100), ".RData", sep=""))

#q("no")

