####################################################
# Plotting Commands related to the particle learning algorithms
# Plot the 95% CI over time: default is to show beta/gamma only
# which.dim -- what to plot
#
plot.ci <- function(saved.stats,trueX=NULL,trueTheta=NULL,trueProp=NULL,
                    col=NULL,plot.diff=0, which.dim=list(x=c(1,2),theta=NULL,prop=NULL),
                    in.legend=NULL,ltype=1,shade.poly=F)
{
   #par(mfcol=c(2,length(which.dim)),mar=c(4,3,2,1)+0.1, oma=c(0,0,0,1))
   XLab <- c("S", "I", "R", "D")
   TLab <- c("S->I", "I->R", "S->R","I->D")
   
   n.comps <- dim(saved.stats)[1]
   if (length(col) == 1 && is.numeric(col))
      col.list <- c(col:(col+n.comps-1))
   else if (length(col) == 0)
      col.list <- c("#bbbbbbcc", "#6666bbcc", "#bb6666cc", "#66bb66cc", "#666666cc") #
   else
      col.list <- col
   
   if (length(ltype)==1)
     ltype = rep(ltype, n.comps)
   
   for (jj in which.dim$x)
   {
      len <- dim(saved.stats)[2]
      if (plot.diff==1) {  # absolute difference
         for (k in 1:n.comps) {
            if (k==1)  
               plot(1:len,(saved.stats[k,,(3*jj-2)]-trueX[,jj]),col=col.list[k],lwd=2,type="l",
                   ylab=XLab[jj], xlab='',lty=ltype[k])
            else 
               lines(1:len,(saved.stats[k,,(3*jj-2)]-trueX[,jj]),col=col.list[k],lwd=2,lty=ltype[k])
            
            lines(1:len,(saved.stats[k,,(3*jj-1)]-trueX[,jj]),col=col.list[k],lty=3)
            lines(1:len,(saved.stats[k,,(3*jj)]-trueX[,jj]),col=col.list[k],lty=3)
          
         }
      }
      else if (plot.diff==2) {  # percent-difference
         for (k in 1:n.comps) {
            if (k==1) 
               plot(1:len,log(saved.stats[k,,(3*jj-2)]/trueX[,jj]),col=col.list[k],lwd=2,type="l",
                   ylab=XLab[jj], xlab='',lty=ltype[k])
            else
               lines(1:len,log(saved.stats[k,,(3*jj-2)]/trueX[,jj]),col=col.list[k],lwd=2,lty=ltype[k])
      
            lines(1:len,log(saved.stats[k,,(3*jj-1)]/trueX[,jj]),col=col.list[k],lty=3)
            lines(1:len,log(saved.stats[k,,(3*jj)]/trueX[,jj]),col=col.list[k],lty=3)
          
         }
      }     
      else { # plot absolute values
         for (k in 1:n.comps) {
            if (k==1)
               plot(1:len,saved.stats[k,,(3*jj-2)],col=col.list[k],type="l",ylab='',main=XLab[jj],lty=3,
                  xlab='',ylim=c(min(saved.stats[,,(3*jj-1)]),max(max(saved.stats[,,(3*jj-2)]),0.97*max(saved.stats[,,(3*jj)]))))
            else
               lines(1:len,saved.stats[k,,(3*jj-2)],col=col.list[k],lty=3)
            
            lines(1:len,saved.stats[k,,(3*jj-1)],col=col.list[k],lty=ltype[k],lwd=2.5)
            lines(1:len,saved.stats[k,,(3*jj)],col=col.list[k],lty=ltype[k],lwd=2.5)

            if (!is.null(trueX))
               lines(1:len,trueX[1:len,jj],col ='red',xlim=c(0,len),lwd=1.5)
            if (n.comps == 1)
              legend("topright", c("Truth", "Median", "95% CI"), col=c("red",col.list[1],col.list[1]), lty=c(ltype[1],3,ltype[1]), lwd=c(1.5,1,2.5))
          }
      }
      if (n.comps > 1 & !is.null(in.legend))
         legend("topright", in.legend, lty=ltype, col=col.list, lwd=2)
   } 
   # loop to plot X-states
   
   
   for (jj in which.dim$theta)
   {
      for (k in 1:n.comps) {
         if (k==1) {
            plot(1:len,saved.stats[k,,(3*(jj+N.STATES)-2)],col=col.list[k],xlab='Period',ylab='',main=TLab[jj],type="l",lty=3,
               ylim=c(min(saved.stats[,,(3*(jj+N.STATES)-1)])*1.02, max(saved.stats[,,3*(jj+N.STATES)]))*0.98)
      #c(trueTheta[jj]*0.75,trueTheta[jj]*1.25)
            if (!is.null(trueTheta))
               abline(h=trueTheta[jj],col='red',lwd=1.5)  
         }
         else 
            lines(1:len,saved.stats[k,,(3*(jj+N.STATES)-2)],col=col.list[k], lty=3)  
          
        
         len <- dim(saved.stats)[2]
         if (shade.poly == T)
             polygon(c(1:len,len:1),c(saved.stats[k,,(3*(jj+N.STATES)-1)],saved.stats[k,len:1,3*(jj+N.STATES)]),col=col.list[k],border=NA)
         else {
            lines(1:len,saved.stats[k,,(3*(jj+N.STATES)-1)],col=col.list[k],lty=ltype[k],lwd=2.5)
            lines(1:len,saved.stats[k,,3*(jj+N.STATES)],col=col.list[k],lty=ltype[k],lwd=2.5)
         }
        
      }   

   } 
   # loop to plot thetas
   
   for (jj in which.dim$prop)
   {
      for (k in 1:n.comps) {
         if (k==1) {
            plot(1:len,saved.stats[k,,(3*(jj+N.STATES+N.RXNS)-2)],col=col.list[k],
                xlab='Period',ylab='',main=paste(TLab[jj]," Sampling"), type="l",lty=3,
               ylim=c(min(saved.stats[,,(3*(jj+N.STATES+N.RXNS)-1)])*1.02, max(saved.stats[,,3*(jj+N.STATES+N.RXNS)]))*0.9)
            if (!is.null(trueProp))
               abline(h=trueProp[jj],col='red',lwd=1.5)  
         }
         else 
            lines(1:len,saved.stats[k,,(3*(jj+N.STATES+N.RXNS)-2)],col=col.list[k],lty=3)  
          
        
         len <- dim(saved.stats)[2]
         if (shade.poly == T & max(saved.stats[k,len:1,3*(jj+N.STATES+N.RXNS)] > 0) )
             polygon(c(1:len,len:1),c(saved.stats[k,,(3*(jj+N.STATES+N.RXNS)-1)],saved.stats[k,len:1,3*(jj+N.STATES+N.RXNS)]),col=col.list[k],border=NA)
         else {
            lines(1:len,saved.stats[k,,(3*(jj+N.STATES+N.RXNS)-1)],col=col.list[k],lty=ltype[k],lwd=2.5)
            lines(1:len,saved.stats[k,,3*(jj+N.STATES+N.RXNS)],col=col.list[k],lty=ltype[k],lwd=2.5)
         }
        
      }
      if (n.comps > 1 & !is.null(in.legend))
         legend("topright", in.legend, lty=ltype, col=col.list, lwd=2)   

   }
}

plot.ci.mcError <- function(in.stats, which.var,col1="#ddddddaa", col2="blue")
#  plots the Monte Carlo variability around the mean/95-quantiles of the posterior for the
#  specified parameter which.var. The shaded region is the min/max across the MC runs in in.stats
#  col1: color of the shaded region; col2: color of the mean line
#
#  First, need to initialize via: plot(c(1,len),y.range,type="n",xaxs="i")
{
  ii <-1
  #for (ii in 1:3) {
    ndx <- which.var*6-3+ii

    len <- dim(in.stats)[2]
    plm <- rbind(apply(in.stats[,,ndx],2, min),apply(in.stats[,,ndx],2,max),apply(in.stats[,,ndx],2,mean))

  #}
    ndx1 <- which.var*6-3+2
    ndx2 <- which.var*6-3+3

    plm2 <- rbind(apply(in.stats[,,ndx1],2,min),apply(in.stats[,,ndx2],2,max),apply(in.stats[,,ndx1],2,mean),apply(in.stats[,,ndx2],2,mean))

    #polygon(c(1:len,len:1),c(plm2[1,1:len],plm2[2,len:1]),col=col1,border=NA)
    lines(1:len,plm2[1,],col=col2,lwd=0.5,lty=2)
    lines(1:len,plm2[2,],col=col2,lwd=0.5,lty=2)

    lines(1:len,plm2[3,],col=col2,lwd=2,lty=2)
    lines(1:len,plm2[4,],col=col2,lwd=2,lty=2)
    polygon(c(1:len,len:1),c(plm[1,1:len],plm[2,len:1]),col=col1,border=NA)
    lines(1:len,plm[3,],col=col2,lwd=3)

}

#######################################
# Compute absolute difference between two cdfs specified as (x,y) lists
cdf.diff <- function( density1, density2)
{
   cdf1 <- cumsum(density1$y)/sum(density1$y)
   step1 <- stepfun(density1$x, c(cdf1,1))
   
   cdf2 <- cumsum(density2$y)/sum(density2$y)
   step2 <- stepfun(density2$x, c(cdf2,1))
   
   ss <- seq( min(min(density1$x),min(density2$x)), max(max(density1$x),max(density2$x)), len=1000)
   
   df <- sum(abs(step1(ss)-step2(ss)))*(ss[2]-ss[1])

   return(list(diff=df,cdf1=step1,cdf2=step2))
}

