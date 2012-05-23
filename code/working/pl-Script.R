source("inference.R")
source("SMLib.R");
require(smcUtils)
source("pfSIR.R")

options(error=recover)
require(plotrix)

#base.params <- list(    initP=c(0.25, 0.25, 0.5, 0.5),
#    initX=c(19980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,1.6,80),
#    trueTheta = array(c(0.75, 0.5, 0, 0.02),dim=c(4,1)) )

N.RXNS <- 4

# initP are the sampling proportions: S->I, I->R, S->R, I->D (in that order)
# initX is the initial state of S/I/R/D
# hyperPrior is the Gamma prior on the 4 thetas: (alpha_i,beta_i)_{i=1}^4
# trueTheta are the actual thetas for S->I, I->R, S->R, I->D
base.params <- list(    initP=c(0.1, 0, 0, 0),
    initX=c(15980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,0.6,80),
    trueTheta = array(c(0.8, 0.5, 0, 0.002),dim=c(4,1)) )



# RUN particleSampledSIR
#
tt <- plSIR(5000, 50, LOOPN=5,verbose="CI")
lw <- particleSampledSIR(5000, 50, LOOPN=5,aLW=0.99,verbose="f",trueX=tt$trueX,Y=tt$Y)  # Liu-West

# compare different methods
tt <- plSIR(5000, 52, LOOPN=10,verbose="f")
lw <- particleSampledSIR(5000, 50, LOOPN=10,aLW=0.99,verbose="f",trueX=tt$trueX,Y=tt$Y) # liu-west
str <- particleSampledSIR(5000, 50, LOOPN=10,aLW=1.99,verbose="f",trueX=tt$trueX,Y=tt$Y) # storvik

plot(c(0,50), c(0.75,0.75), type="l",col="red",xlab="Periods t",ylab="Beta Posterior",lwd=2)
plot.ci.mcError(str$stat,1,col1="#ffeedd99", col2="black")
plot.ci.mcError(tt$stat,1,col1="#eeddff99", col2="blue")
plot.ci.mcError(lw$stat,1,col1="#ddffee99", col2="green")

############
# To see a more realistic example
temp <- plSIR(5000, 70, LOOPN=10,verbose="HIST")
plot.ci(temp$stat[1,,],temp$trueX,temp$theta,1)
plot.ci.mcError(temp$stat,1)


