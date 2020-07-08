# - - - - - - - - - - - - - - - - - 
# Analysis of tracing
# - - - - - - - - - - - - - - - - - 

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


# Plot functions ----------------------------------------------------------


plotbypc <- function(R=1.2, k=0.3, q=1, d = 0.5, ps=seq(0.0,1,0.1), cs=seq(0.2,1,0.2), noplot=F, panel=rep("",3)){
    psmat <- matrix(ps,length(ps),length(cs))
    csmat <- matrix(cs,length(ps),length(cs),byrow=T)
    tertiary <- tertiary_cases(R,k,q,psmat,csmat,csmat,d)
    if(noplot){return(tertiary)} # skip plotting
    
    cols1=terrain.colors(length(cs)+2)[1:length(cs)]
    cols2=terrain.colors(length(cs)+2)[1:length(cs)]
    cols3=terrain.colors(length(cs)+2)[1:length(cs)]
    
    pars = paste0(" R = ",R,", k = ",k,", b = ",format(q,nsmall=1),", d = ",d,paste0(rep(" ",10),collapse=""))
     matplot(tertiary$forward_avert/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols1,lwd=2,ylab="effectiveness (forward)",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c(pars,paste(c(panel[1],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$fwd_back_avert/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols2,lwd=2,ylab="effectiveness (forward + backward)",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c("",paste(c(panel[2],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$backward_gain/tertiary$base,x=ps,type="o",pch=19,lty=1,col=cols3,lwd=2,ylab="increment with backward tracing",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c("", paste(c(panel[3],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    return(tertiary)
}

