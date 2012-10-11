sckm = list()
sckm$s = 3 # species (S,I,R)
sckm$r = 2 # reactions (S->I, I->R)
#                   S -> I    I -> R
sckm$Pre  = rbind( c(1,1,0), c(0,1,0))
sckm$Post = rbind( c(0,2,0), c(0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(100,2,0)
sckm$theta = c(2/100,1)

set.seed(1)
y = tau.leap(sckm, 15)

ld = 2
clrs = c("green","red","blue")
plot(y[,1], type="l", ylim=c(0,sckm$X[1]), lwd=ld, col=clrs[1], xlab="Time", ylab="Number")
lines(y[,2], lwd=ld, col=clrs[2])
lines(y[,3], lwd=ld, col=clrs[3])
a <- readline("Press <enter> to continue:")


