
fl = system("ls *sum.csv",T)

d = list()
for (i in 1:length(fl)) d[[i]] = read.csv(fl[i])

library(plyr)
library(reshape2)
d = dcast(rbind.fill(d), parameter + time + method ~ statistic)

d$parameter = factor(d$parameter, c("S","I","R","p: S->I","p: I->R","blank","r: S->I","r: I->R"))

library(ggplot2)
qplot(time, MCV, data=d, geom="line", colour=method, facets=~parameter)




