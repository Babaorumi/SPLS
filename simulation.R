## NOTATIONS 

## n: sample size
## d: dimension (number of variables)
## K: number of clusters

## ARGUMENTS

## n: the number of observations to generate
## Q: either a list of the quantile functions G_j^{-1}, j=1,...,d,
## either a list of (n*1)-vectors: (\tilde x_1,...,\tilde x_n) where the
## \tilde x_j are a sample of G_j.
## pz: numeric vector of size K giving the cluster weights
## theta: numeric vector of size K giving the copula parameters.
## WARNING: only one real parameter per copula!
## copulaFamilies: character vector giving the names of the copula families.
## Can be "gaussian", "frank", "clayton" or "indep".
## psi: a function. MUST be the identity function function(x){x}
## otherwise the code is not guaranteed to work.

## MORE DETAILS

## For more details, see the paper at
## https://hal.archives-ouvertes.fr/hal-01263382
## See also the file demo.R

simul <- function(n,Q,mu,pz,theta,copulaFamilies,psi=NULL){

    d <- nrow(mu)
    K <- ncol(mu)
    x <- matrix(nrow=n,ncol=d)
    z <- numeric(K)
    cop <- list()
    intervalTheta <- matrix(nrow=2,ncol=K)

    if(is.null(psi)){
        psi <- function(t,j){t}
    }
    
    for(z in 1:K){
        if(copulaFamilies[z]=="gaussian"){
            cop[[z]] <- normalCopula(dim=d)
            intervalTheta[,z] <- c(-.99,.99)
        }else if(copulaFamilies[z]=="clayton"){
            cop[[z]] <- claytonCopula(dim=d)
            intervalTheta[,z] <- c(0.001,50)
        }else if(copulaFamilies[z]=="frank"){
            cop[[z]] <- frankCopula(dim=d)
            intervalTheta[,z] <- c(-50,50)
        }else if(copulaFamilies[z]=="indep"){
            cop[[z]] <- indepCopula(dim=d)
            intervalTheta[,z] <- c(0,0)
        }else{
            stop("'copulaFamily' is wrong.")
        }
    }

    for(i in 1:n){
        z[i] <- sample(1:K,1,prob=pz)
        cop[[z[i]]]@parameters <- theta[z[i]]
        u <- rCopula(1,cop[[z[i]]])
        for(j in 1:d){
            if(class(Q[[j]])=="function"){
                x[i,j] <- psi(Q[[j]](u[j])+mu[j,z[i]],j)
            }else{
                x[i,j] <- psi(quantile(Q[[j]], u[j])+mu[j,z[i]],j)
            }
        }
    }
    
    return( cbind(x,z) )
}
