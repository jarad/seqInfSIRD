library(tlpl)

load("PL.RData")

n.sims = length(pl)
plq = list()

for (i in 1:n.sims)
{
  plq[[i]] = tlpl_quantile(pl[[i]])
}

save(plq, file="PLq.RData")

q("no")

