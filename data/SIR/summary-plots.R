

fl = system("ls *sum.csv",T)

d = list()
for (i in 1:length(fl)) d[[i]] = read.csv(fl[i])

library(plyr)
library(reshape2)
d = dcast(rbind.fill(d), parameter + time + method ~ statistic)

d$parameter = factor(d$parameter, c("S","I","R","p: S->I","p: I->R","blank","r: S->I","r: I->R"))

library(ggplot2)
qplot(time, MCV, data=d, geom="line", colour=method, facets=~parameter)

qplot(time, MSE, data=d, geom="line", colour=method, facets=~parameter)
qplot(time, MSE, data=d[d$parameter%in% levels(d$parameter)[c(4:5,7:8)],], geom="line", colour=method, facets=~parameter)


d$tmp = d$MSE+2*d$seMSE

qplot(time, tmp, data=d[d$parameter%in% levels(d$parameter)[c(4:5,7:8)],], geom="line", colour=method, facets=~parameter)

ggplot(d[d$parameter=="p: S->I",], aes(x = time, y = MSE, group=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = MSE-2*seMSE, ymax = MSE+2*seMSE, fill = "05%-95%"), alpha = .25)


