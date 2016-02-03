## It is assumed that the support of G is (-\infty,+\infty) and that
## G is centered at 0. The copulas must have only one parameter.
## effInterval is the effective support of G (from a computer
## point of view).
## Q is either a list of the quantile functions G^{-1} either
## a list of (n*1)-vectors : (\tilde x_1,...,\tilde x_n) where the
## \tilde x_j are a sample of G_j.
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


## G: array of size n*d*K representing G_j(y_{ij}-mu_{jz}) for
## j=1,...,d, i=1,...,n and z=1,...,K (same for g);
## y_{ij}=\psi_j^{-1}(x_{ij})
## mu (d*K)-matrix
## psiInvPrime is a (n*d)-matrix containing the (\psi_j^{-1})'(x_{ij})

EMalgo <- function(data, copulaFamilies, psi, psiInv, psiInvPrime, nbit, bandwidth){

    x <- data
    n <- nrow(x)
    d <- ncol(x)
    K <- length(copulaFamilies)
    km <- kmeans(x,centers=K)
    mu <- t(km$centers)
    pz <- km$size/nrow(x)
    theta <- numeric(K) ## NOT SURE THIS IS NECESSARY
    G <- array(dim=c(n,d,K))
    g <- array(dim=c(n,d,K))
    intervalTheta <- matrix(nrow=2,ncol=K)
    cop <- list()
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
            stop("'copulaFamilies' is wrong.")
        }
    }
    sampleOfG <- NA
    for(j in 1:d){
        for(z in 1:K){
            sampleOfG <- c( sampleOfG, psi(x[km$cluster==z,j],j)-mu[j,z] )
            sampleOfG <- sampleOfG[-1]
            G[,j,z] <- ks::kcde(sampleOfG,
                                eval.points=psiInv(x[,j],j)-mu[j,z])$estimate
            g[,j,z] <- ks::kde(sampleOfG,
                               eval.points=psiInv(x[,j],j)-mu[j,z])$estimate

            
        }
    }

    for(z in 1:K){
        if(copulaFamilies[z]=="indep"){
            cop[[z]]@parameters <- theta[z] <- 0
        }else{
            theta[z] <- fitCopula(cop[[z]], pobs(x[km$cluster==z,]),
                                  method="itau")@estimate
            cop[[z]]@parameters <- theta[z]
        }
    }
    
    objectiveValue <- 0
    
    ## loop
    for(t in 1:nbit){
        
        ## 0. Calcul de h(z|x^i), x^i multivariate (OK)
        num <- hzx <- matrix(nrow=n, ncol=K)
        denom <- numeric(n)
        for(z in 1:K){
            for(i in 1:n){ 
                num[i,z] <- dCopula(G[i,,z],copula=cop[[z]]
                                    )*prod(g[i,,z])*prod(psiInvPrime[i,])*pz[z]
            }
        }
        for(i in 1:n){
            denom[i] <- sum(num[i,],na.rm=TRUE)
            for(z in 1:K){
                hzx[i,z] <- if(denom[i]==0) 1/K else num[i,z]/denom[i]
            }
        }

        ## 1. Update pi_z for each z (OK)
        newpz <- numeric(K)
        for(z in 1:K){
            newpz[z] <- mean(hzx[,z],na.rm=TRUE)
        }

        ## 2. Update mu_{jz} for each j=1,...,d, z=1,...,K (OK)
        newmu <- matrix(nrow=d,ncol=K)
        for(z in 1:K){
            for(j in 1:d){
                newmu[j,z] <- sum(psiInv(x[,j],j)*hzx[,z])/sum(hzx[,z])
            }
        }
        
        ## Update G, ie create a sample of it
        zsim <- numeric(n)
        xtilde <- matrix(nrow=n,ncol=d)
        for(i in 1:n){
            zsim[i] <- sample(1:K,1,prob=hzx[i,]) 
            for(j in 1:d){
                xtilde[i,j] <- psiInv(x[i,j],j)-mu[j,zsim[i]] 
            }
        }

        ## Compute the "G" (OK)
        newG1 <- newG2 <- newg1 <- newg2 <- array(dim=c(n,d,K))
        newG <- newg <- array(dim=c(n,d,K))
        for(z in 1:K){
            for(j in 1:d){
                newG1[,j,z] <- ks::kcde(xtilde[,j],
                                        eval.points=psiInv(x[,j],j)-mu[j,z])$estimate
                newG2[,j,z] <- ks::kcde(xtilde[,j],
                                        eval.points=-psiInv(x[,j],j)+mu[j,z])$estimate
                newg1[,j,z] <- ks::kde(xtilde[,j],
                                       eval.points=psiInv(x[,j],j)-mu[j,z])$estimate
                newg2[,j,z] <- ks::kde(xtilde[,j],
                                       eval.points=-psiInv(x[,j],j)+mu[j,z])$estimate
                ## symmetrization
                newg[,j,z] <- (newg1[,j,z]+newg2[,j,z])/2
                newG[,j,z] <- (newG1[,j,z]+1-newG2[,j,z])/2
            }
        }

        ## Update copulas parameters
        newTheta <- numeric(K)
        for(z in 1:K){
            if(copulaFamilies[z]=="indep"){
                cop[[z]]@parameters <- newTheta[z] <- 0
            }else{
                optimizeFoo <- function(par,z){
                    cop[[z]]@parameters <- par
                    sum(log(dCopula(newG[,,z],copula=cop[[z]]))*hzx[,z])
                }
                newTheta[z] <- optimize(
                    optimizeFoo,intervalTheta[,z],z=z,maximum=TRUE)$maximum
                cop[[z]]@parameters <- newTheta[z]
            }
        }

        ## Compute the new objective value
        term1 <- term2 <- term3 <- matrix(nrow=n,ncol=K)
        for(z in 1:K){
            for(i in 1:n){
                term1[i,z] <- log(dCopula(newG[i,,z],copula=cop[[z]]))*hzx[i,z]
                term2[i,z] <- sum(log(newg[i,,z]))*hzx[i,z]
                term3[i,z] <- log(newpz[z])*hzx[i,z]
            }
        }

        objectiveValue <- c(objectiveValue,sum(term1+term2+term3))
        print(objectiveValue[length(objectiveValue)])
        
        ## Updates the programme variables
        mu <- newmu
        g <- newg
        G <- newG
        pz <- newpz
        theta <- newTheta
    }
    
    return( list(G=G,g=g,theta=theta,mu=mu,pz=pz,xtilde=xtilde,
                 hzx=hzx, copulaFamilies=copulaFamilies,
                 Ovalue=objectiveValue[-1]) )
}


