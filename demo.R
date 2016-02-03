library(ks) ## for nonparametric estimation of the symmetric densities
library(copula) ## for handling copulas

source("simulation.R")
source("estimation.R")

## Simualte 300 observations of a 2-dimensional SPLS copula-based mixture model
## with Student and Gaussian symmetric densities and Gaussian copulas (3 groups)

dataGroup <- simul(300,
                   Q=list(function(x){qt(x,df=4)},
                          function(x){qnorm(x,sd=2)}),
                   mu=matrix(nrow=2,ncol=3,
                             c( 0,  3,
                               3,   0,
                               -3, 0)),
                   pz=rep(1/3,3),
                   theta=c(.5,0,-.5),
                   copulaFamilies=rep("gaussian",3))

head(dataGroup)

plot(dataGroup[,1:2]) ## plot the raw data

plot(dataGroup[dataGroup[,3]==1,1:2]) ## plot Group 1 only

## Estimation

data <- dataGroup[,1:2]

em <- EMalgo(data, copulaFamilies=rep("gaussian",3),
                  psi=function(t,j){t},
                  psiInv=function(t,j){t},
                  psiInvPrime=matrix(rep(1,nrow(data)*ncol(data)),
                                     nrow=nrow(data),ncol=ncol(data)),
                  nbit=50, bandwidth=NULL)

em

plot(1:length(em$Ovalue),em$Ovalue,type="b",
     ylab="conditional log-likelihood",
     xlab="steps in the EM algorithm") ## objective function to increase.
## If the plot does not show an upward trend,
## we have doubts about our estimates...

## estimated parameters:

em$mu ## estimated location parameters

dim(em$xtilde)
head(em$xtilde) ## sample of the symmetric distributions g_1 and g_2

em$theta ## estimated copula parameters

em$pz ## estimated cluster weighs
