# - - - - - - - - - - - - - - - - - 
# Analysis of tracing
# - - - - - - - - - - - - - - - - - 

library(tidyverse)

# setwd("~/Documents/GitHub/2020-ncov-forward-back-tracing/")


# Define functions --------------------------------------------------------

secondary_cases <- function(r=1.2, # reproduction number
                            k=0.3, # dispersion
                            q=1, # probability cluster identified
                            p=0.5, # probability contact traced
                            c1=0.5, # scaled reduction from tracing index case
                            c2=0.5, # scaled proportional reduction from tracing contact of cluster case
                            d=0.5 # probability detection of cases
                           ){
  
  # No tracing
  baseline_secondary <- r^3*(1+1/k)
  
  # Cases with forward tracing
  forward_secondary <- (r*d)*(r)*(r*p*c1) + (r^3*(1-d))*(1+1/k)

  # Forward + backward tracing
  backward_secondary <- (r*d)*(r)*(r*p*c1) + (r*(1-d)*q)*(r)*(r*p*c2)*(1+1/k)

  c(base=baseline_secondary,forward=forward_secondary,forward_back=backward_secondary)

}


secondary_cases()

# Plot functions ----------------------------------------------------------

kk <- seq(0.1,0.5,0.1)

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
plot(data_infections,col="black",lwd=2,yaxs="i",ylab="number",xlab="days",ylim=c(0,1200))
points(data_deaths,col="red",lwd=2)
#lines(estimated_cases,col="blue",lwd=2)
#lines(estimated_cases2,col="purple",lwd=2)
