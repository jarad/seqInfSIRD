
d = read.csv("summary-statistics.csv")

library(ggplot2)
qplot(time, value, data=d[d$stat=="MSE",], facets=~parameter)




