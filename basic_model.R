# - - - - - - - - - - - - - - - - - 
# Analysis of tracing
# - - - - - - - - - - - - - - - - - 

library(tidyverse)

# setwd("~/Documents/GitHub/2020-ncov-forward-back-tracing/")


# Define functions --------------------------------------------------------

tertiary_cases <- function(R=1.2, # reproduction number
                            k=0.3, # dispersion
                            q=1, # probability cluster identified
                            p=0.5, # probability contact traced
                            c1=0.5, # proportional reduction from tracing index case
                            c2=0.5, # proportional reduction from tracing contact of cluster 
                            d=0.5 # probability detection of cases
                           ){
  tertiary_avert <- function(R, k, q, p, c1, c2, d){
      R^2*p*c1 + (1-(1-q*p)*(1-d))*p*c2*R^3*(1+1/k)
  }
  
  # No tracing
  baseline_tertiary <- tertiary_avert(R,k,1,1,1,1,1)
  
  # Cases averted with forward tracing
  forward_tertiary_avert <- tertiary_avert(R,k,0,p,c1,c2,d)

  # Cases averted with forward + backward tracing
  forward_backward_tertiary_avert <- tertiary_avert(R,k,q,p,c1,c2,d)

  list(base=baseline_tertiary,
    forward_avert=forward_tertiary_avert,
    fwd_back_avert=forward_backward_tertiary_avert,
    backward_gain=forward_backward_tertiary_avert-forward_tertiary_avert)
}


#secondary_cases()

# Plot functions ----------------------------------------------------------

#kk <- seq(0.1,0.5,0.1)

plotbypc <- function(R=1.2, k=0.3, q=1, d = 0.5, ps=seq(0.0,1,0.1), cs=seq(0.2,1,0.2), noplot=F){
    psmat <- matrix(ps,length(ps),length(cs))
    csmat <- matrix(cs,length(ps),length(cs),byrow=T)
    tertiary <- tertiary_cases(R,k,q,psmat,csmat,csmat,d)
    if(noplot){return(tertiary)} # skip plotting
    
    cols1=terrain.colors(length(cs)+2)[1:length(cs)]
    cols2=terrain.colors(length(cs)+2)[1:length(cs)]
    cols3=terrain.colors(length(cs)+2)[1:length(cs)]
    
    matplot(tertiary$forward_avert/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols1,lwd=2,ylab="proportion cases averted (forward tracing)",xlab="proportion of contacts traced",ylim=c(0,1.02))
    legend("topleft",legend=rev(paste("reduction if traced = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$fwd_back_avert/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols2,lwd=2,ylab="proportion cases averted (forward + backward tracing)",xlab="proportion of contacts traced",ylim=c(0,1.02))
    legend("topleft",legend=rev(paste("reduction if traced = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$backward_gain/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols3,lwd=2,ylab="additional averted with backward tracing",xlab="proportion of contacts traced",ylim=c(0,1.02))
    legend("topleft",legend=rev(paste("reduction if traced = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    return(tertiary)
}

par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(2,0.6,0))
plotbypc()

dev.copy(png,paste0("plot1.png"),units="cm",width=20,height=12,res=150)
#dev.copy(pdf,paste0(dir_pick,"Figure_1.pdf"),width=6,height=8)
dev.off()

