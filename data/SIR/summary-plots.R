
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



ggplot(d, aes(x = MSE, group = 1)) +
  geom_ribbon(aes(ymin = quantile.0.05, ymax = quantile.0.95, fill = "05%-95%"), alpha = .25) + 
  geom_ribbon(aes(ymin = quantile.0.25, ymax = quantile.0.75, fill = "25%-75%"), alpha = .25) +
  geom_line(aes(y = quantile.0.5)) +
  scale_fill_manual(name = "", values = c("25%-75%" = "red", "05%-95%" = "blue")) 



