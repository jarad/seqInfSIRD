

fl = system("ls *sum.csv",T)

d = list()
for (i in 1:length(fl)) d[[i]] = read.csv(fl[i])

library(plyr)
library(reshape2)
d = dcast(rbind.fill(d), parameter + time + method ~ statistic)

d$parameter = factor(d$parameter, c("S","I","R","p: S+I->2I","p: I->R","blank","r: S+I->2I","r: I->R"))

library(ggplot2)

# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





d = d[d$time>0,]

qplot(time, MCV, data=d, geom="line", colour=method, facets=~parameter)

q1 = qplot(time, sqrt(MSE), data=d[d$parameter%in% levels(d$parameter)[1:3],], geom="line", colour=method, facets=~parameter)
q2 = qplot(time, sqrt(MSE), data=d[d$parameter%in% levels(d$parameter)[4:5],], geom="line", colour=method, facets=~parameter)
q3 = qplot(time, sqrt(MSE), data=d[d$parameter%in% levels(d$parameter)[7:8],], geom="line", colour=method, facets=~parameter)


pdf("mse.pdf")
multiplot(q1,q2,q3)
dev.off()



d$tmp = d$MSE+2*d$seMSE

qplot(time, tmp, data=d[d$parameter%in% levels(d$parameter)[c(4:5,7:8)],], geom="line", colour=method, facets=~parameter)





ggplot(d[d$parameter=="p: S+I->2I",], aes(x = time, y = MSE, colour=method)) +
  geom_ribbon(aes(ymin = MSE-2*seMSE/10, ymax = MSE+2*seMSE/10, fill = "05%-95%"), alpha = .25)


