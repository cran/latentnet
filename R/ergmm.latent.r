#############################################################################
###  ergmm.latent: R function prepares the R data to be passed into the C  ###
###  function, calls the C function to estimate the latent space model    ###
###  and then puts the C data back into readable R storage.               ###
#############################################################################
###                Description of variables                               ###
#############################################################################
### Arguments:
###
###      gY: a network object representing a social network
###
### dimSpace: The dimension of the latent space
###
###       p: The number of covariates in the model
###
###       X: A list of matrices with each matrix containing the covariate
###          differences between each of the actors.
###
### MCMCSampleSize: number of MCMC iterations to preform
###
###  burnin: number of iterations to allow the chain to converge
###
### interval: interval between saving the MCMC interation
###
### z.delta: The proposal distribution for all of the positions is  a
###          normal distribution with standard deviation equal to z.delta
###          centered around the precvious position.
###
### z.prior.mu: The prior distriution for each of the positions is a normal
###          distribution centered this prior mu value
###
### z.prior.sd: The prior distribution for each of the positions is  a
###          normal distribtuion with this standard deviation
###
### b.delta: This is a vector of length p (if not of length p   the values
###          will be recycled to fill a vector of length p). Each value
###          represents the standard deviation on the  of the propsal
###          value for that coefficent.  The proposal  distribution for
###          the betas is a normal centered around the  previous
###          distribution.
###
### b.prior.mu: This is a vector of length p with the prior distribution
###          for the betas is normal with this mean.  If this is not a
###          vector of length p, values will be recylced or truncated to
###          fill a p length vector 
###
### b.prior.mu: This is a vector of length p with the prior distribution
###          for the betas is normal with this standard deviation.  If
###        this is not a vector of length p,values will ne recycled of
###          truncated to fill a p length vector.
#############################################################################
###                    Description of Return Value                        ###
#############################################################################

ergmm.latent <- function(gY, dimSpace=2, p=0, X=NULL, theta0=NULL,
                        MCMCSampleSize=1000, burnin=0, interval=1,
                        z.delta=0.1, z.prior.mu=0, z.prior.sd=10,
                        b.delta=0.5, b.prior.mu=0, b.prior.sd=10,
                        penalty.sigma=c(10,0.5), maxit=200, trace=0,
                        MLEonly=FALSE,
                        verbose=FALSE, ...)
{
   Y <- as.sociomatrix(gY)
#  Ydesign <- get.network.attribute(gY,"design")
#  if(!is.null(Ydesign)){
#    Ydesign <- sociomatrix(Ydesign)
#  }
   directed <- is.directed(gY)
   efit <- network.layout.fruchtermanreingold(Y, NULL)
   D <- as.matrix(dist(efit))
#
#  Make the maximum distance close to the rest
#
   reach <- D!=Inf
   D[D==Inf] <- max(D[D!=Inf])+1
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
    abz.list <- list(Y=Y,dp=p,X=X,Z=Z,nnodes=g,dimSpace=dimSpace,
                     reach=reach,directed=directed)
    abvZ.mds <- optim(par=vector(length=p,mode="numeric"),
                   fn=mlpYmdsZ, gr=mlpYmdsZ.grad,
                   method="BFGS", 
                   control=list(fnscale=-1, maxit=maxit,trace=trace),
                   abz.list=abz.list)
   }else{
    abvZ.mds <- list(par=rep(0,length=p))
   }
   if(verbose){trace <- 4}else{trace <- 0}
   if(!is.null(theta0) && length(theta0)==p){
    beta <- theta0
   }else{
    ## start the coeffs for beta to 0
    beta <- abvZ.mds$par
   }

   cat(paste("The latent space has dimension",dimSpace,"\n"))
   if(!MLEonly){
    cat(paste("b.prior.mu",b.prior.mu,""))
    cat(paste("b.prior.sd",b.prior.sd,"\n"))
    cat(paste("z.prior.mu",z.prior.mu,""))
    cat(paste("z.prior.sd",z.prior.sd,"\n"))
    if(penalty.sigma[1]>0){
     cat(paste("A distance penalty sigma of",penalty.sigma[1],"was used\n"))
    }else{
     cat(paste("No distance penalty was used\n"))
    }
   }

   if(dimSpace>0){
    # this can be speeded up with c and also what to do when get
    # infinite gradient?????????
    ## there is a problem that needs to be addressed sometimes the R optim
    ##   function breaks on finding the mle of certain networks but if called
    ##   a few more times it will work - very mysterious didnt spend much time
   
    ## now find the mle
    abvZ <- c(beta,Z)

    #simulated annealing to find rough mle
#   abvZ <- optim(par=abvZ,fn=mlpY,Y=Y,dp=p,X=X,nnodes=g,k=dimSpace,
#                method="SANN",
#                control=list(fnscale=-1, maxit=maxit))$par
    #BFGS  use previous found mle to plug into quasi newton raphson
#   cat("Calling MLE fit\n")
#   MLE.fit0 <- optim(par=abvZ,fn=mlpY,Y=Y,dp=p,X=X,nnodes=g,k=dimSpace,
#                method="BFGS",
#                control=list(fnscale=-1, maxit=maxit,trace=trace))
#   cat("Calling MLE fit\n")
#   abvZ <- MLE.fit0$par
#
#   msh: ADDED gradient information
#
    #BFGS  use previous found MDS fit to plug into quasi newton raphson
    if(penalty.sigma[1]>0){
      penalty.factor <- c(1/(penalty.sigma[1]*penalty.sigma[1]),
                          penalty.sigma[2])
    }else{
      penalty.factor <- c(0,penalty.sigma[2])
    }
    abz.list <- list(Y=Y,dp=p,X=X,nnodes=g,dimSpace=dimSpace,
                     penalty.factor=penalty.factor,reach=reach,
                     directed=directed)
    cat("Calling MLE fit\n")

    MLE.fit <- try(
                optim(par=abvZ,fn=mlpY, gr=mlpY.grad,
                 method="BFGS",
                 control=list(fnscale=-1, maxit=maxit,trace=trace),
                 abz.list=abz.list)
                 )
    if(inherits(MLE.fit,"try-error")){
     warning("MLE could not be found.")
     MLE.like <- abvZ.mds$value
     MLE.iterations <- abvZ.mds$iter
    }else{
     abvZ <- MLE.fit$par
     MLE.like <- MLE.fit$value
     MLE.iterations <- MLE.fit$iter
    }
    Z.mle <- matrix(abvZ[-(1:p)],nrow=g,ncol=dimSpace)
    if(MLEonly){
     beta.mle <- as.vector(abvZ[1:p],mode="numeric")
     Results<-list(Z.mle=Z.mle,beta.mle=beta.mle,
                Llik=MLE.like, mle.like=MLE.like,
                Hessian=diag(beta.mle-beta.mle+1),Beta=matrix(beta.mle,nrow=1),
                Beta.rate=MLE.fit$counts[1],
                Z=array(Z.mle,dim=c(dim(Z.mle),1)),
                Z.rate=MLE.fit$counts[2], newnetwork=gY
                )
     return(Results)
    }
   }else{
    abvZ <- beta
    abz.list <- list(Y=Y,dp=p,X=X,Z=NULL,nnodes=g,dimSpace=dimSpace,
                     reach=reach,directed=directed)
    abvZ.0 <- optim(par=abvZ,fn=mlpY0,gr=mlpY0.grad,
                    hessian=TRUE, method="BFGS",
                    control=list(fnscale=-1, maxit=maxit,trace=trace),
                    abz.list=abz.list)
    abvZ <- c(abvZ.0$par, abvZ.0$value)
   }
   # assign the mle values to proper spot 
   beta.mle <- as.vector(abvZ[1:p],mode="numeric")
   if(dimSpace==0){
    # Return all the results of the chain
     Results<-list(Z.mle=NULL,beta.mle=beta.mle,
                Llik=abvZ[p+1],
                mle.like=MLE.like,
                Hessian=solve(-abvZ.0$hessian),Beta=beta.mle,
                Beta.rate=abvZ.0$counts[1],
                Z=NULL,Z.rate=abvZ.0$counts[2], newnetwork=gY
                )
     return(Results)
   }
   # else fit the model with MCMC
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
   Beta <-  vector( mode="numeric",length=n.sample*p)
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
     b.delta <-as.vector(matrix(b.delta, nrow=p, ncol=dimSpace), mode="numeric")
   }
   storage.mode(b.delta) <- "double"
#
   edgelist <- as.matrix.network(gY, matrix.type="edgelist") - 1
#
   n.edges <- gY$gal$mnext-1 
   heads<-edgelist[,1]  # C indexes starting at 0
   tails<-edgelist[,2]
   if(directed){
    directed <- 1
   }else{
    directed <- 0
   }
  
   # get the Z.mle into vector form
   vZ <- abvZ[-(1:p)] 
   storage.mode(vZ) <- "double"

   # put the covariates into a vector
   vX <- as.numeric(unlist(X))
   storage.mode(vX) <- "double"

    # call the C code with the proper arguements
   cat("Calling MCMC_latent\n")
   ret <- .C("MCMC_latent_wrapper",  
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
          as.integer(p),                  as.integer(directed),
          as.double(vX),                  as.double(Llike),
          PACKAGE="latentnet"
          )

   # Take the positions out of vector form and put into a list of 
   # matrices
#
#  MSH: keep it simple
#
   Zp <- array(ret[[15]],dim=c(g,dimSpace,n.sample))
#
#  New Procrustes of the fit including reflections
#
   Zp<-array(apply(Zp,3,ergmm.procAdj,a=NULL,fa=ergmm.procAdj.fcnt(Z.mle)),
             dim=dim(Zp))

  # Return all the results of the chain
   Results<-list(Z.mle=Z.mle,beta.mle=beta.mle,
              Llik=ret[[24]],
              Beta=ret[[19]],
              Beta.rate=ret[[20]],
              Z=Zp,Z.rate=ret[[18]], newnetwork=gY, mle.like = MLE.like,
              mle.iterations=MLE.iterations
              )
   Results
}

##############################################################################

##############################################################################
# MINUS log prob of network, with Z in vector form, 
# With Z[1]=alpha and the rest as Z[1,1], Z[1,2], Z[2,1]..........
# to be used in by optim or nlm
mlpY<-function(abvZ, abz.list)
{   
  dp <- abz.list$dp
# dp <- 1
  penalty.factor <- abz.list$penalty.factor
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  beta <- abvZ[1:dp]

  # put Z in matrix form
  Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
     for(i in 2:dp){
       tX <- tX + (abz.list$X[[i]] * beta[i])
     }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
# diag(eta)<-0
  # Return llk
# sum( abz.list$Y*eta - log( 1+exp(eta) ) ) + nnodes*log(2)
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
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
  loglik - (diamZ) * penalty.factor[1]
}

mlpY.grad<-function(abvZ, abz.list)
{   

  penalty.factor <- abz.list$penalty.factor
  dp <- abz.list$dp
# dp <- 1
  nnodes <- abz.list$nnodes
  dimSpace <- abz.list$dimSpace
  # extract parameters
  beta <- abvZ[1:dp]

  # put Z in matrix form
  Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
     for(i in 2:dp){
       tX <- tX + (abz.list$X[[i]] * beta[i])
     }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0

  probeta <- exp(eta)/(1+exp(eta))
  probeta[is.na(probeta)] <- 1
  Dlliknuij <-  abz.list$Y - probeta
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
  diagYnr <- diagY & !abz.list$reach
  diagY <- diagY & abz.list$reach

  Dllik <- vector(length=dp+nnodes*dimSpace)

  cd <- 0
  for(i in 1:dp){
   cd <- cd + 1
   M <- abz.list$X[[i]]
   Dllik[cd] <- sum(Dlliknuij[diagY] * M[diagY])
  }
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

mlpYmdsZ<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
# dp <- 1
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  beta <- abvZ[1:dp]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
     for(i in 2:dp){
       tX <- tX + abz.list$X[[i]] * beta[i]
     }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
# diag(eta)<-0
  # Return llk
# sum( abz.list$Y*eta - log( 1+exp(eta) ) ) + nnodes*log(2)
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta) < col(eta)
  }
  diagY <- diagY & abz.list$reach
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  logexpeta <- log( 1+exp(eta) )
  logexpeta[is.na(logexpeta)] <- eta
  loglik <- sum( Y*eta - logexpeta )
  loglik
}

mlpYmdsZ.grad<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
# dp <- 1
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  beta <- abvZ[1:dp]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
     for(i in 2:dp){
       tX <- tX + (abz.list$X[[i]] * beta[i])
     }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0

  probeta <- exp(eta)/(1+exp(eta))
  probeta[is.na(probeta)] <- 1
  Dlliknuij <-  abz.list$Y - probeta
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
  diagY <- diagY & abz.list$reach

  Dllik <- vector(length=length(beta))

  cd <- 0
  for(i in 1:dp){
   cd <- cd + 1
   M <- abz.list$X[[i]]
   diag(M)<-0
   Dllik[cd] <- sum(Dlliknuij[diagY] * M[diagY])
  }

  # Return the gradient of the llk
  Dllik
}

mlpYnull<-function(abvZ,Y,dp,X,nnodes,k=2)
{   
  # extract parameters
  beta <- abvZ[1:dp]

  # put Z in matrix form
  Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=k)

  # find the covariate contribution
  tX <- X[[1]] * beta[1]
  if(dp>1){
    for(i in 2:dp){
      tX <- tX + (X[[i]] * beta[i])
    }
  }

  # see above for formula
  eta <- tX
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0
  # Return llk
  logexpeta <- log( 1+exp(eta) )
  logexpeta[is.na(logexpeta)] <- eta
  -(  sum( Y*eta - logexpeta ) + dim(Y)[1]*log(2)  )
}

mlpY0<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
# dp <- 1
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  # extract parameters
  beta <- abvZ[1:dp]

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
    for(i in 2:dp){
      tX <- tX + abz.list$X[[i]] * beta[i]
    }
  }

  # see above for formula
  eta <- tX
# diag(eta)<-0
  # Return llk
# -(  sum( Y*eta - log( 1+exp(eta) ) ) + dim(Y)[1]*log(2)  )
  if(abz.list$directed){
   diagY <- row(eta)!=col(eta)
  }else{
   diagY <- row(eta)< col(eta)
  }
  diagY <- diagY & abz.list$reach
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  logexpeta <- log( 1+exp(eta) )
  logexpeta[is.na(logexpeta)] <- eta
  sum( Y*eta - logexpeta )
}

mlpY0.grad<-function(abvZ, abz.list)
{   
  dp <- abz.list$dp
# dp <- 1
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  # extract parameters
  beta <- abvZ[1:dp]

  Deta<-matrix(0, nrow=nnodes*nnodes, ncol=dp)
  cd <- 0
  for(i in 1:dp){
   cd <- cd + 1
   M <- abz.list$X[[i]]
   diag(M)<-0
   Deta[,cd] <- M
  }

  # find the covariate contribution
  tX <- abz.list$X[[1]] * beta[1]
  if(dp>1){
     for(i in 2:dp){
       tX <- tX + abz.list$X[[i]] * beta[i]
     }
  }

  # see above for formula
  eta <- tX
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0
  # Return the gradient of the llk
  probeta <- exp(eta)/(1+exp(eta))
  probeta[is.na(probeta)] <- 1
  aaa <-  matrix(abz.list$Y - probeta, nrow=1)
  aaa <- as.vector(aaa %*% Deta)
  aaa
}

##############################################################################
##############################################################################
#  Gives the negative distance between nodes
## This uses the distance formula and fancy tricks
## for dim k: d[i,j] <- sqrt( (x1i-x2j)^2 + ..... + (xki-xkj)^2 ) 
## ==> 
##  d[i,j]^2 = x1i^2 + x2j^2 +...+ x1k^2 + x2k^2 - 2*x1i*x1j -...- 2*x1k*x2k
## that is what the matrix manipulations below do.
lpz.dist<-function(Z)
{
  if(!is.matrix(Z)){
   Z <- matrix(Z,ncol=1)
  }
  # the diagonals of this will be the sum of square of all the cordinates
  # for each node.  We want to add all these together for each point.
  ZtZ<-Z%*%t(Z)
  mg<-as.matrix(diag(ZtZ))%*%rep(1,dim(Z)[1]) 

  # adding together all of the squared distances of each pair
  mg<-mg+t(mg)
   
  # This takes the square root of the the squares minus 2 times coords
  d<-sqrt(abs(mg-2*ZtZ))
  #returns the negative distance
  -d             
}

###############################################################################
   "ergmm.procAdj.fcnt" <- function(a){sweep(a,2,apply(a,2,mean))}
   "ergmm.procAdj" <- function(b,a,fa=ergmm.procAdj.fcnt(a)){
#
#   Ordinary Procustes analysis of b onto a: 
#    using translation, rotation and reflections by least squares.
#
    fb   <- ergmm.procAdj.fcnt(b)
    Rsvd <- svd(t(fa) %*% fb)
    R <- Rsvd$v %*% t(Rsvd$u)
#   returns the estimated orthogonal rotation matrix
    fb %*% R
   }
