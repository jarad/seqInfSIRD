library(tlpl)

load("PL.RData")
source("quantile-settings.R")

plq = llply(pl, tlpl_quantile, probs=probs, verbose=0, 
            .progress = progress_text(style=ifelse(interactive(), 3, 1)))

save(plq, file="PLq.RData")

q("no")

