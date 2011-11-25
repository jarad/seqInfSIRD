#Use particles to filter for SIR model with discrete-time sampled 
#transitions 
particleSampledSIR <- function(dt, N, T, LOOPN=1,aLW=0.99,verbose="NONE",col="blue",trueX=NULL,Y=NULL)
{
    initP <- c(0.25, 0.25, 0.5, 0.5)
    initX <- c(19980,20, 0,0)
    hyperPrior <- c(22.5,30,15,30,0,1000,20,1000)
    trueTheta <- array(c(0.75, 0.5, 0, 0.02),dim=c(4,1))
    saved.stats <- array(0, dim=c(LOOPN,ceil(T/dt)+1,24))

    if (is.null(trueX)) # generate a scenario
    {
       trueX <- array(0, dim=c(ceil(T/dt)+1,4))
       Y <- array(0,dim=c(ceil(T/dt)+1,4))
       
       trueX[1,] <- initX

       # Fix the true state trajectory for replicability
       #set.seed(21)
       # Construct the true process and the observations
       for (i in 1:(T/dt)) {
           out <- tauLeap(t(t(trueX[i,])),t(t(trueTheta)),NULL,t(initP))
           Y[i+1,] <- out$Y
           trueX[i+1,] <- out$X
       }
       #browser()
    }

    
    for (loop in 1:LOOPN) {
    
       # Initialize particles
       X <- t(array(rep(initX, N),dim=c(N.RXNS,N)))
       lambda <- t(array(rep(initP, N), dim=c(N.RXNS,N)))
       Suff <- t(array(rep(hyperPrior,N), dim=c(2*N.RXNS,N))) # sufficient conjugate statistics
       fixTheta <- array(0,dim=c(N,N.RXNS))

       for (jj in 1:N.RXNS)
          fixTheta[,jj] <- rgamma(N, hyperPrior[2*jj-1],hyperPrior[2*jj])
       
       theta <- fixTheta
       p.weights <- array(1, N)/N

       curt <- 0 
       i <- 1
       totalWeight <- 1
       saved.stats[loop,i,] <- saveStats(X, fixTheta)

       while (curt < T-dt & trueX[i,2] >0) {
         if (aLW > 1)  # Storvik filter: sample from the posterior Gamma mixture
           for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
           
         else # use the Liu-West thetas
            theta <- fixTheta
        
        out <- tauLeap(t(X), t(theta),Suff,lambda,Y[i+1,])
        X <- t(out$X); dX <- t(out$dX); Suff <- out$hyper
        
        p.weights <- mutate(dX,Y[i+1,],lambda,p.weights)
        curt <- (i-1)*dt

         if ( i %%10  == 0 & verbose =="HIST") # plot posterior of infectiousness parameter
         {
            
            gridx <- seq(0.5,1,by=0.005)
            gridy <- array(0, length(gridx))
            for (jj in 1:length(gridx))
               gridy[jj] <- sum(dgamma(gridx[jj],Suff[,1],Suff[,2]))
         
            if (aLW < 1)
              hist(theta[,1],25,freq=F,main=sprintf("t=%d",curt),xlab='S->I Rate') 
            else
               plot(gridx,gridy/N,type="l",col=col,main=sprintf("t=%d",curt),xlab='S->I Rate')
         
            abline(v=trueTheta[1], col="red")
         }
         totalWeight <- totalWeight*sum(p.weights)

         # resample
         ESS <- 1/sum(N^2*p.weights^2) 
         CoV <- sqrt(sum( (N*p.weights-1)^2) ) 
         Ent <- sum( p.weights*log(N*p.weights))
         newIndex <- multinomial.resample(p.weights) # branchMinVar
         #browser()
         X <- X[newIndex,] 
         Suff <- Suff[newIndex,]
         lambda <- lambda[newIndex,]
         fixTheta <- fixTheta[newIndex,]
         
         # Liu and West move
         meanLam <- varLam <- array(0,N.RXNS)
         if (aLW < 1) # & curt > 1
         {
            for (jj in 1:4) {
              meanLam[jj] <- mean(fixTheta[,jj])
              varLam[jj] <- var(fixTheta[,jj])
              fixTheta[,jj] <- pmax(0, aLW*fixTheta[,jj] + (1-aLW)*meanLam[jj] +
                  sqrt( (1-aLW^2)*varLam[jj])*rnorm(N))
            }
          }
          p.weights <- rep(1/N, N)
        

        i <- i+1
        saved.stats[loop,i,] <- saveStats(X,theta)
        # end of main loop
      }
   }
   
   if (verbose=='CI')
      plot.ci(saved.stats[1,1:i,],trueX[1:i,],trueTheta,1,col)
 
   return( list(stat=saved.stats,trueX=trueX,Y=Y))
}

################### branching function
#resample using Crisan min-variance method
# branching from Crisan (2006) p .10
branchMinVar <- function(p.weights)
{
        N <- length(p.weights)
        branch <- p.weights/mean(p.weights)
        fracBranch <- branch - floor(branch)
        ub <- runif(N-1)
        newNdx <- array(N, dim=c(N,1))
        curNdx <- 1
        gb <- N
        hb <- N
        for (j in 1:(N-1))
        {
          if (fracBranch[j] + gb - branch[j] - floor(gb - branch[j]) < 1)
            ob <- floor(branch[j]) + (ub[j] > 1-fracBranch[j]/(gb-floor(gb)+1e-8))*(hb - floor(gb))
          else
            ob <- floor(branch[j]) + (hb - floor(gb)) + (1-hb+floor(gb))*(ub[j] < 1 - (1-fracBranch[j])/(1-gb+floor(gb)))

          gb <- gb - branch[j]
          if (ob == 0)
            next
          
          hb <- hb - ob
          newNdx[curNdx:(curNdx+ob-1)] <- j
          curNdx <- curNdx + ob
        }
        #browser()
        # last index saved already as N

        return (newNdx)
}

#####################################################
# move one step of the SIR
# X is a matrix: each column has 4 rows for SIRD states 
# theta is a matrix, each column has 4 rows for SIRD rates 
tauLeap<- function(X, theta, hyper=NULL, prop,curY=NULL)
{
    h <- t(X)
    dX <- X
    hyper2 <- array(1, dim=c(N.RXNS,2,dim(X)[2]))
    
    out <- one.step.C(X, hyper2, theta, t(prop), sample=T)
    
    for (j in 1:dim(X)[2])
    {    
       #out <- simulate.one.step.C(X[,j],theta[,j],sum(X[,j]) )
       #dX[,j] <- out$dX
       #X[,j] <- out$X
       h[j,]  <- hazard.R(X[,j], sum(X[,j]))   
    }
       
    if (is.null(curY))
        Y <- rbinom(N.RXNS, out$newX[,1], prop)
    else
        Y <- curY

    if (!is.null(hyper)) 
      for (i in 1:N.RXNS) {
        hyper[,2*i-1] <- hyper[,2*i-1] + Y[i]
        hyper[,2*i] <- hyper[,2*i] + prop[,i]* h[,i]
      }
    

    return(list(X=out$X,dX=out$newX,hyper=hyper,Y=Y))
}

mutate <- function(dX, curY, prop, weights)
{
    newWeights <- weights
   
    for (jj in 1:N.RXNS)
      newWeights <- newWeights*dbinom(curY[jj], dX[,jj], prop[,jj])
      
    return(newWeights)

}

saveStats <- function(X, theta)
{
    summ.stat <- array(0,dim=c(24,1))
    for (jj in 1:4)
    {
      summ.stat[(jj-1)*6+1] <- mean(X[,jj])
      summ.stat[((jj-1)*6+2):((jj-1)*6+3)] <- quantile(X[,jj],c(0.025,0.975))

      summ.stat[(jj-1)*6+4] <- mean(theta[,jj])
      summ.stat[((jj-1)*6+5):((jj-1)*6+6)] <- quantile(theta[,jj],c(0.025,0.975))
    }
    return(summ.stat)
}

plot.ci <- function(saved.stats,trueX,trueTheta,sir.plotCI=0,col="blue")
{
   par(mfcol=c(2,3),mar=c(4,4,2,1), oma=c(0,0,0,1))
   XLab <- c("S", "I", "R", "D")
   TLab <- c("S->I", "I->R", "S->R","I->D")
   for (jj in c(1,2,4))
   {
       len <- dim(saved.stats)[1]
       plot(1:len,saved.stats[,(6*jj-5)],col=col,lwd=2,type="l",ylab=XLab[jj],
        xlab='',ylim=c(min(saved.stats[,(6*jj-4)]),max(saved.stats[,(6*jj-3)])))
       lines(1:len,saved.stats[,(6*jj-4)],col=col,lty=3)
       lines(1:len,saved.stats[,(6*jj-3)],col=col,lty=3)

       lines(1:len,trueX[,jj],col ='red',xlim=c(0,len),lwd=1.5)
       plot(1:len,saved.stats[,(6*jj-2)],col=col,xlab='time',ylab=TLab[jj],type="l",
       ylim=c(trueTheta[jj]*0.8,trueTheta[jj]*1.2), lwd=2 )
       abline(h=trueTheta[jj],col='red',lwd=1.5)
       if (sir.plotCI == 1) {
         lines(1:len,saved.stats[,(6*jj-1)],col=col,lty=3)
         lines(1:len,saved.stats[,6*jj],col=col,lty=3)
         
       }
   }
}
