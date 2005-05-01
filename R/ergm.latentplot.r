ergm.latentplot <- function(gY, dimSpace=2, theta0=NULL,
                        MCMCSampleSize=1000, burnin=0, interval=1,
                        z.delta=0.1, z.prior.mu=0, z.prior.sd=10,
                        b.delta=0.5, b.prior.mu=0, b.prior.sd=10,
                        penalty.sigma=c(10,0.5), maxit=200, trace=0,
                        nsubsample=100, verbose=TRUE)
{
   Y <- sociomatrix(gY)
   directed <- is.directed(gY)
#  Ydesign <- get.network.attribute(gY,"design")
#  if(!is.null(Ydesign)){
#    Ydesign <- sociomatrix(Ydesign)==0
#  }
   D <- ergm.geodesicmatrix(gY)
#
#  Make the maximum distance close to the rest
#
   reach <- D!=Inf
   D[D==Inf] <- max(D[D!=Inf])+1
#
#  Make the maximum distance close to the rest
#
   D[D==gY$gal$n] <- max(D[D!=gY$gal$n])+1
#
   g <- gY$gal$n
   if(dimSpace>0){
    cat("Calling MDS\n")
    Z <- cmdscale(D,dimSpace)
    vZ <- as.vector(Z)
    Z <- matrix(vZ,nrow=g)
    ## Center the positions around the origin
    Z <- sweep(Z,2,apply(Z,2,mean))
#
#   Find the solution with space from MDS positions
#
    abz.list <- list(Y=Y,Z=Z,nnodes=g,dimSpace=dimSpace,
                     reach=reach,directed=directed)
    abvZ.mds <- optim(par=0,
                   fn=mlpYmdsZ.plot, gr=mlpYmdsZ.grad.plot,
                   method="BFGS", 
                   control=list(fnscale=-1, maxit=maxit,trace=trace),
                   abz.list=abz.list)
   }else{
    abvZ.mds <- list(par=rep(0,length=1))
   }
   if(verbose){trace <- 4}else{trace <- 0}
   if(length(theta0)==1){
    beta <- theta0
   }else{
    ## start the coeffs for beta to 0
    beta <- abvZ.mds$par
   }

   cat(paste("The latent space has dimension",dimSpace,"\n"))
   cat(paste("b.prior.mu",b.prior.mu,""))
   cat(paste("b.prior.sd",b.prior.sd,"\n"))
   cat(paste("z.prior.mu",z.prior.mu,""))
   cat(paste("z.prior.sd",z.prior.sd,"\n"))
   if(penalty.sigma[1]>0){
     cat(paste("A distance penalty sigma of",penalty.sigma[1],"was used\n"))
   }else{
     cat(paste("No distance penalty was used\n"))
   }

   ## now find the mle
   abvZ <- c(beta,Z)

   #BFGS  use previous found MDS fit to plug into quasi newton raphson
   if(penalty.sigma[1]>0){
     penalty.factor <- c(1/(penalty.sigma[1]*penalty.sigma[1]),
                         penalty.sigma[2])
   }else{
     penalty.factor <- c(0,penalty.sigma[2])
   }
   abz.list <- list(Y=Y,nnodes=g,dimSpace=dimSpace,
                    penalty.factor=penalty.factor,reach=reach,
                    directed=directed)
   cat("Calling MLE fit\n")

   MLE.fit <- try(
               optim(par=abvZ,fn=mlpY.plot, gr=mlpY.grad.plot,
                method="BFGS",
                control=list(fnscale=-1, maxit=maxit,trace=trace),
                abz.list=abz.list)
                )
   if(inherits(MLE.fit,"try-error")){
     warning("MLE could not be found.")
     MLE.like <- abvZ.mds$value
   }else{
     abvZ <- MLE.fit$par
     MLE.like <- MLE.fit$value
   }
   Z.mle <- matrix(abvZ[-1],nrow=g,ncol=dimSpace)
#
#  Now set up the Bayesian MCMC
#
   # assign the mle values to proper spot 
   beta.mle <- as.vector(abvZ[1],mode="numeric")
   storage.mode(beta.mle) <- "double"
   # get dimensions for the vectors
   lz <- gY$gal$n * dimSpace
#  n.sample <- ceiling( (MCMCSampleSize-burnin)/interval ) 
   n.sample <- MCMCSampleSize
   MCMCSampleSize <- ceiling( MCMCSampleSize*interval + burnin) 
   if(n.sample<1){
     stop (paste("Invalid MCMC sample size and/or burnin."))
   }

   # make the vectors which will return all of the MCMC iteration values for 
   # the parameters
   Beta <-  vector( mode="numeric",length=n.sample)
   Llike <- vector( mode="numeric", length=n.sample )
   vZ.post <- vector(mode="numeric",length=lz*n.sample)
   storage.mode(Beta) <- "double"
   storage.mode(Llike) <- "double"
   storage.mode(vZ.post) <- "double"

   # make the vectors which will return the rates of acceptance 
   Z.rate <-  vector( mode="numeric", length=n.sample )
   B.rate <-  vector( mode="numeric", length=n.sample )
   storage.mode(Z.rate) <- "double"
   storage.mode(B.rate) <- "double"

   # make sure that there is a vector of length p filled with the deltas for each 
   # beta.  should it be vector or the same for each beta?  i say diff
   if(length(b.delta)==1){
     b.delta <-as.vector(matrix(b.delta, nrow=1, ncol=dimSpace), mode="numeric")
   }
   storage.mode(b.delta) <- "double"
#
   n.edges <- gY$gal$mnext-1 
   edgelist <- as.matrix.network(gY, matrix.type="edgelist") - 1
   heads<-edgelist[,1]  # C indexes starting at 0
   tails<-edgelist[,2]
##
## Krista: crate the edgelist for the missing (like this above)
## 
#   edgelistmis <- as.matrix.network(Ydesign, matrix.type="edgelist") - 1
#   headsmis<-edgelist[,1]  # C indexes starting at 0
#   tailsmis<-edgelist[,2]
##
   if(directed){
    directed <- 1
   }else{
    directed <- 0
   }
  
   # get the Z.mle into vector form
   vZ <- abvZ[-1] 
   storage.mode(vZ) <- "double"

   # call the C code with the proper arguements
   cat("Calling MCMC_latent\n")
   ret <- .C("MCMC_latent_wrapperplot",  
          as.integer(heads),              as.integer(tails),    
          as.integer(n.edges),            as.integer(gY$gal$n),
          as.integer(MCMCSampleSize),     as.integer(burnin),
          as.integer(interval),           as.integer(dimSpace),
          as.double(z.prior.mu),          as.double(z.prior.sd),
          as.double(b.prior.mu),          as.double(b.prior.sd),
          as.double(vZ),                  as.double(z.delta),
          as.double(vZ.post),             as.double(Z.rate),              
          as.double(beta.mle),            as.double(b.delta),
          as.double(Beta),                as.double(B.rate),
          as.integer(directed),
          as.double(Llike),               as.integer(nsubsample),
          PACKAGE="latentnet"
          )

   # Take the positions out of vector form and put into a list of 
   # matrices
#
#  MSH: keep it simple
#
   Zp <- array(ret[[15]],dim=c(g,dimSpace,n.sample))

  # Return all the results of the chain
   Results<-list(Z.mle=Z.mle,beta.mle=beta.mle,
              Llik=ret[[22]],
              Beta=ret[[19]],
              Beta.rate=ret[[20]],
              Z=Zp,Z.rate=ret[[16]], newnetwork=gY, mle.like = MLE.like
              )
   Results
}

##############################################################################

##############################################################################
# MINUS log prob of network, with Z in vector form, 
# With Z[1]=alpha and the rest as Z[1,1], Z[1,2], Z[2,1]..........
# to be used in by optim or nlm
mlpY.plot<-function(abvZ, abz.list)
{   
  penalty.factor <- abz.list$penalty.factor
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Ydesign <- abz.list$Ydesign #K#
  # extract parameters
  beta <- abvZ[1]

  # put Z in matrix form
  Z<-matrix( abvZ[-1],nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # see above for formula
  eta <- beta + lpZ
  # Return llk
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
##
## Select only the non-missing values
## #here
  if(!is.null(Ydesign)){ #K#
   diagY <- diagY & abz.list$Ydesign #K#
  }  #K#
  diagYnr <- diagY & !abz.list$reach
  diagY <- diagY & abz.list$reach
  Y <- abz.list$Y[diagY]
  logexpeta <- log( 1+exp(eta) )
  logexpeta[is.na(logexpeta)] <- eta
  eta <- eta[diagY]
# print(c(sum(Y*eta), sum(- logexpeta[diagY]), - 0.5*sum(logexpeta[diagYnr])))
#
  loglik <- sum( Y*eta - logexpeta[diagY] )
  loglik <- loglik - penalty.factor[2]*sum(logexpeta[diagYnr])
#
# Note this uses a penalized distance that
# makes the mean distance small
#
# loglik - 1*(-mean(lpZ[diagY]))
#
# This should be right, but this works better
# diamZ <- dimSpace*mean(Z*Z)
# but this works better
  diamZ <- mean(Z%*%t(Z))
#
# diamZ.nr <- -mean(lpZ[diagYnr])
# diamZ    <- -mean(lpZ[diagY])
# diamZ <- mean(Z%*%t(Z))
# loglik - (diamZ - 50*diamZ.nr) * penalty.factor
# loglik - (diamZ - diamZ.nr) * penalty.factor
# print(c(loglik, penalty.factor, diamZ))
# print(loglik)
  loglik - (diamZ) * penalty.factor[1]
}

mlpY.grad.plot<-function(abvZ, abz.list)
{   

  penalty.factor <- abz.list$penalty.factor
  nnodes <- abz.list$nnodes
  dimSpace <- abz.list$dimSpace
  Ydesign<- abz.list$Ydesign
  # extract parameters
  beta <- abvZ[1]

  # put Z in matrix form
  Z<-matrix( abvZ[-1],nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # see above for formula
  eta <- beta + lpZ

  probeta <- exp(eta)/(1+exp(eta))
  probeta[is.na(probeta)] <- 1
  Dlliknuij <-  abz.list$Y - probeta
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
##
## Select only the non-missing values
##
  if(!is.null(Ydesign)){ #K#
   diagY <- diagY & abz.list$Ydesign #K#
  } #K#
  diagYnr <- diagY & !abz.list$reach
  diagY <- diagY & abz.list$reach

  Dllik <- vector(length=1+nnodes*dimSpace)

  cd <- 1
  Dllik[cd] <- sum(Dlliknuij[diagY])
  for(j in 1:dimSpace){
   for(i in 1:nnodes){
    v <- (Z[i,j]-Z[,j])/lpZ[i,]
    v[is.infinite(v)|is.na(v)] <- 0
    M <- matrix(0,ncol=nnodes,nrow=nnodes)
    M[i,] <- v
    M[,i] <- v
    cd <- cd + 1
    Dllik[cd] <- sum(Dlliknuij[diagY] * M[diagY])
#
#   Add penalty
#
#   print(c(penalty.factor, sum(Dlliknuij[diagYnr]*M[diagYnr]),Z[i,j]/nnodes))

    Dllik[cd] <- Dllik[cd] + penalty.factor[2]*sum(Dlliknuij[diagYnr]*M[diagYnr])
    Dllik[cd] <- Dllik[cd] - 2*penalty.factor[1]*Z[i,j]/nnodes
   }
  }

  # Return the gradient of the llk
  Dllik
}

mlpYmdsZ.plot<-function(abvZ,abz.list)
{   
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  beta <- abvZ[1]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # see above for formula
  eta <- beta + lpZ
  # Return llk
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta) < col(eta)
  }
##
## Krista:  Select only the non-missing values
##
#  if(!is.null(Ydesign)){
#   diagY <- diagY & abz.list$Ydesign
#  }
  diagY <- diagY & abz.list$reach
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  logexpeta <- log( 1+exp(eta) )
  logexpeta[is.na(logexpeta)] <- eta
  loglik <- sum( Y*eta - logexpeta )
  loglik
}

mlpYmdsZ.grad.plot<-function(abvZ,abz.list)
{   
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  beta <- abvZ[1]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # see above for formula
  eta <- beta + lpZ

  probeta <- exp(eta)/(1+exp(eta))
  probeta[is.na(probeta)] <- 1
  Dlliknuij <-  abz.list$Y - probeta
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
##
## Select only the non-missing values
##
  if(!is.null(Ydesign)){
   diagY <- diagY & abz.list$Ydesign
  }
  diagY <- diagY & abz.list$reach

  Dllik <- vector(length=length(beta))

  cd <- 1
  Dllik[cd] <- sum(Dlliknuij[diagY])

  # Return the gradient of the llk
  Dllik
}
