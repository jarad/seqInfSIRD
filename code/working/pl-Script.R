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
base.params <- list(    initP=c(0.05, 0, 0, 0),
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


### Try to handle the actual weekly data set: 59 weeks total
Nweek <- 59
harare <- read.csv("data/harare/harareClean.csv")
harY <- array(0,dim=c(Nweek,4))
harY[,1] <- diff(harare$total)
hararePL <- plSIR(5000,Nweek,LOOPN=1,Y=harY,verbose="HIST")
plot.ci(hararePL$stat[1,,],NULL,hararePL$theta,1)
hararePL <- plSIR(5000,Nweek,LOOPN=5,Y=harY,verbose="CI")
plot(c(0,60),c(0.6,0.6),type="l",col="red",ylab="Beta posterior",xlab="Periods")
plot.ci.mcError(hararePL$stat,1,col1="#ddddddaa", col2="blue")

# The prior on GAMMA has a strong impact on the POSTERIOR of BETA:
base.params$hyperPrior <- c(1,1.2,20,40,0,1000,0,1000)
hararePL3 <- plSIR(5000,Nweek-1,LOOPN=10,Y=harY,verbose="CI")
plot(c(0,60),c(0.6,0.6),type="l",col="red",ylab="Beta posterior",xlab="Periods")
plot.ci.mcError(hararePL$stat,1,col1="#ddddddaa", col2="blue") # close to 0.5

base.params$hyperPrior <- c(2.25,3,1,2,0,1000,0,1000)  # close to 1 due to the wide prior of Gamma
hararePL <- plSIR(5000,Nweek-1,LOOPN=10,Y=harY,verbose="CI")

################################################################
# Generate a scenario
N.week <- 60
scen.params <- base.params
scen.params$initP <- c(1,1,1,1) # use full observations first
set.seed(5)
sampScen <- generate.scenario(scen.params,tauLeap,N.week) 
sampScen$theta
# [1] 0.7159145 0.5000000 0.0000000 0.0020000
plot(sampScen$X[,2]/2,type="l")
lines(sampScen$Y[,1],col="green") 
scen.params$trueTheta <- sampScen$theta 

# now generate various observations by sub-sampling
# case 1: typical good observations
Ymed <- array(0,dim=c(N.week+1,4))
Ymed[,1] <-rbinom(N.week+1,sampScen$Y[,1],0.2)
scen.params$initP <- c(0.2,0,0,0)
scenMed <- plSIR(8000,N.week,LOOPN=1,Y=Ymed,trueX=sampScen$X,verbose="CI",model.params=scen.params)

# typical low observations (like in Zimbabwe)
Ylow <- array(0,dim=c(N.week+1,4))
Ylow[,1] <-rbinom(N.week+1,sampScen$Y[,1],0.05)
scen.params$initP <- c(0.05,0,0,0)
scenLow <- plSIR(8000,N.week,LOOPN=1,Y=Ylow,trueX=sampScen$X,verbose="CI",model.params=scen.params)

# 20%-observations of infecteds transitioning to recovereds as well
Ygam <- Ymed
Ygam[,2] <- rbinom(N.week+1,sampScen$Y[,2],0.2)
scen.params$initP <- c(0.2,0.2,0,0)
scenGamma <- plSIR(8000,N.week,LOOPN=1,Y=Ygam,trueX=sampScen$X,verbose="CI",model.params=scen.params)

# observe 20% of infecteds as well as ALL deaths (very implicit info about current infecteds)
Ydeaths <- Ymed
Ydeaths[,4] <- sampScen$Y[,4]
scen.params$initP <- c(0.2,0,0,1)
scenDeaths <- plSIR(8000,N.week,LOOPN=1,Y=Ydeaths,trueX=sampScen$X,verbose="CI",model.params=scen.params)

########################
# To run SimsSIRD.RData examples
sims.params <- list(    initP=c(0.05, 0, 0, 0),
    initX=c(16000,2, 0,0),    hyperPrior = c(15,10,5,10,0.1,10,0.5,10),
    trueTheta = array(c(0.8, 0.5, 0, 0.002),dim=c(4,1)) )
N.week <- 52

for (j in 4:4)
{
  sims.params$initP <- probs[j,]
  sims.params$initP[3:4] <- 0
  sims.params$trueTheta <- gammas[j,]
  simY <- as.matrix(sims[[j]]$y)
  simY[,3:4] <- 0
  simX <- as.matrix(sims[[j]]$x)
  simPL <- plSIR(8000,N.week-1,LOOPN=1,Y=simY,trueX=simX,verbose="CI",model.params=sims.params)
}



