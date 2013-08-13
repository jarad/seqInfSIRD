library(tlpl)
library(plyr)

load("PL.RData")
source("settings.R")

plq = llply(pl, tlpl_quantile, probs=probs, verbose=0, 
            .progress = progress_text(style=ifelse(interactive(), 3, 1)))

save(plq, file="PLq.RData")

q(ifelse(interactive(),"ask","no"))

