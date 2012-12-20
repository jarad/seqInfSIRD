library(tlpl)

load("sims.RData")
source("filter-settings.R")

pl = list()

for (i in 1:n.sims)
{
  d = list(y=t(data[[i]]), tau=1)
  prior$X = sckm$X
  pl[[i]] = tlpl(d, sckm, prior=prior, n.particles=n.particles)
}

save(pl, file="PL.RData")

q("no")

