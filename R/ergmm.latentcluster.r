#############################################################################
###  ergmm.latentCluster: R function prepares the R data to be passed into ###
###  the C function, calls the C function to estimate the latent space    ###
###  model and then puts the C data back into readable R storage.         ###
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
### b.prior.sd: This is a vector of length p with the prior distribution
###          for the betas is normal with this standard deviation.  If
###        this is not a vector of length p,values will ne recycled of
###          truncated to fill a p length vector.
###
### epsilon:  - to be specified within?
###
### mu:  - to be specified within?
###
### Sigma:  - to be specified within?
###
### Ki:  - to be specified within?
###
#############################################################################
###                    Description of Return Value                        ###
#############################################################################


ergmm.latentcluster <- function(gY, dimSpace=2, p=0, X=NULL, theta0=NULL, ng = 1,
                        MCMCSampleSize=1000, burnin=0, interval=1,
                        z.delta=0.1, z.prior.mu=0, z.prior.sd=10,
                        b.delta=0.5, b.prior.mu=0, b.prior.sd=10,
                        maxit=200,Sigprior = 10, muSigprior = 75, dirprior = 1,
                        alphaprior = 1, chisqprop = 2, thetaprop = 0.1,
                        verbose=FALSE, User.mle=NULL,Z=NULL)
{
  if(verbose){trace <- 4}else{trace <- 0}
  Y <- as.matrix.network(gY,matrix.type="adjacency")
  D <- DistMatrix(Y)
  g <- gY$gal$n
  if(is.null(Z))
    if(dimSpace>0){
      cat("Calling MDS\n")
      Z <- cmdscale(D,dimSpace)
      vZ <- as.vector(Z)
      Z <- matrix(vZ,nrow=g)
      ## Center the positions around the origin
      Z <- sweep(Z,2,apply(Z,2,mean))
    }
#
#  Find the solution with space from MDS positions
#
  abz.list <- list(Y=Y,dp=p,X=X,Z=Z,nnodes=g,dimSpace=dimSpace)
#   abvZ.mds <- optim(par=vector(length=p,mode="numeric"),
#                     fn=mlpY.clustermdsZ, gr=mlpY.clustermdsZ.grad,
#                     method="BFGS", 
#                     control=list(fnscale=-1, maxit=maxit),
#                     abz.list=abz.list)
  if(p>0)
  {
    abvZ.mds <- try(optim(par=vector(length=p,mode="numeric"),
                      fn=mlpY.clustermdsZ, gr=mlpY.clustermdsZ.grad,
                      method="BFGS", 
                      control=list(fnscale=-1, maxit=maxit),
                      abz.list=abz.list))
    if(inherits(abvZ.mds,"try-error"))
    {
      abz.list <- list(Y=Y,dp=p,X=X,Z=Z+rnorm(dim(Z)[1]),nnodes=g,dimSpace=dimSpace)
      abvZ.mds <- optim(par=vector(length=p,mode="numeric"),
                        fn=mlpY.clustermdsZ, gr=mlpY.clustermdsZ.grad,
                        method="BFGS", 
                        control=list(fnscale=-1, maxit=maxit),
                        abz.list=abz.list)
    }
  }
  else
    abvZ.mds <- optim(par=vector(length=1,mode="numeric"),
                      fn=mlpY.clustermdsZ, gr=mlpY.clustermdsZ.grad,
                      method="BFGS", 
                      control=list(fnscale=-1, maxit=maxit),
                      abz.list=abz.list)    
  if((length(theta0)==p) && (p>0)){
    beta <- theta0
  }else{
    ## start the coeffs for beta to 0
    beta <- abvZ.mds$par
  }
  
  cat(paste("The latent space has dimension",dimSpace,"\n"))
  cat(paste("b.prior.sd",b.prior.sd,"\n"))
  cat(paste("z.prior.sd",z.prior.sd,"\n"))

  if(dimSpace>0){
    # this can be speeded up with c and also what to do when get
    # infinite gradient?????????
    ## there is a problem that needs to be addressed sometimes the R optim
    ##   function breaks on finding the mle of certain networks but if called
    ##   a few more times it will work - very mysterious didnt spend much time
   
    ## now find the mle
    if(p>0)
      abvZ <- c(beta,Z)
    else abvZ <- as.vector(Z)

    #simulated annealing to find rough mle
#  abvZ <- optim(par=abvZ,fn=mlpY.cluster,Y=Y,dp=p,X=X,nnodes=g,k=dimSpace,
#                method="SANN",
#                control=list(fnscale=-1, maxit=maxit))$par
   #BFGS  use previous found mle to plug into quasi newton raphson
#  cat("Calling MLE fit\n")
#  MLE.fit0 <- optim(par=abvZ,fn=mlpY.cluster,Y=Y,dp=p,X=X,nnodes=g,k=dimSpace,
#                method="BFGS",
#                control=list(fnscale=-1, maxit=maxit,trace=trace))
#  cat("Calling MLE fit\n")
#  abvZ <- MLE.fit0$par
#
#  msh: ADDED gradient information
#
    #BFGS  use previous found MDS fit to plug into quasi newton raphson
    abz.list <- list(Y=Y,dp=p,X=X,nnodes=g,dimSpace=dimSpace)
    cat("Calling MLE fit\n")
#     MLE.fit <- optim(par=abvZ,fn=mlpY.cluster,gr=mlpY.cluster.grad,
#                      method="BFGS",
#                      control=list(fnscale=-1, maxit=maxit,trace=trace),
#                      abz.list=abz.list)
    MLE.fit <- try(optim(par=abvZ,fn=mlpY.cluster,gr=mlpY.cluster.grad,
                         method="BFGS",
                         control=list(fnscale=-1, maxit=maxit,trace=trace),
                         abz.list=abz.list))
    if(inherits(MLE.fit,"try-error"))
    {
      if(p>0)
        abvZ <- c(beta,Z+rnorm(dim(Z)[1]))
      else abvZ <- as.vector(Z+rnorm(dim(Z)[1]))

      MLE.fit <- optim(par=abvZ,fn=mlpY.cluster,gr=mlpY.cluster.grad,
                       method="BFGS",
                       control=list(fnscale=-1, maxit=maxit,trace=trace),
                       abz.list=abz.list)
    }

    abvZ <- MLE.fit$par
    MLE.like <- MLE.fit$value

   # scale back down usually too big 
   # MSH: do not do this
#  abvZ <- (abvZ) * 2 /abvZ[1] 
#  abvZ <- (abvZ) * 2 /abvZ.0[1]
    if(p>0)
      Z.mle <- matrix(abvZ[-(1:p)],nrow=g,ncol=dimSpace)
    else
      Z.mle <- matrix(abvZ,nrow=g,ncol=dimSpace)
  }else{
    abvZ <- beta
    abz.list <- list(Y=Y,dp=p,X=X,Z=Z,nnodes=g,dimSpace=dimSpace)
    abvZ.0 <- optim(par=abvZ,fn=mlpY.cluster0,gr=mlpY.cluster0.grad,
                    hessian=T, method="BFGS",
                    control=list(fnscale=-1, maxit=maxit,trace=trace),
                    abz.list=abz.list)
    abvZ <- c(abvZ.0$par, abvZ.0$value)
  }
   # assign the mle values to proper spot
  if(p>0)
    beta.mle <- as.vector(abvZ[1:p],mode="numeric")
  else
    beta.mle <- 0
  if(dimSpace==0){
    # Return all the results of the chain
    Results<-list(Z.mle=NULL,beta.mle=beta.mle,
                  Llik=abvZ[p+1],
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
  n.sample <- MCMCSampleSize
  MCMCSampleSize <- ceiling( MCMCSampleSize*interval + burnin) 
  if(n.sample<1){
    stop (paste("Invalid MCMC sample size and/or burnin."))
  }

  # make the vectors which will return all of the MCMC iteration values for 
  # the parameters
  KiList <- vector( mode="numeric",length=n.sample*g)
  muList <- vector( mode="numeric",length=n.sample*dimSpace*ng)
  SigmaList <- vector( mode="numeric",length=n.sample*ng)
  Llike <- vector( mode="numeric", length=n.sample )
  vZ.post <- vector(mode="numeric",length=lz*n.sample)
  storage.mode(muList) <- "double"
  storage.mode(SigmaList) <- "double"
  storage.mode(Llike) <- "double"
  storage.mode(vZ.post) <- "double"

   # make the vectors which will return the rates of acceptance 
  Z.rate <-  vector( mode="numeric", length=n.sample )
  B.rate <-  vector( mode="numeric", length=n.sample )
  storage.mode(Z.rate) <- "double"
  storage.mode(B.rate) <- "double"

   # make sure that there is a vector of length p filled with the deltas for each 
   # beta.  should it be vector or the same for each beta?  i say diff

   # get the edge list ready to pass in.  if the network is undirected included the
   # ij edge as well as the ji edge
  edgelist <- as.matrix.network(gY, matrix.type="edgelist") - 1
  if(gY$gal$directed)
    dir <- 1
  else
    dir <- 0
  n.edges <- gY$gal$mnext-1 
  heads<-edgelist[,1]  # C indexes starting at 0
  tails<-edgelist[,2]
  
   # get the Z.mle into vector form
  if(p>0)
    vZ <- abvZ[-(1:p)]
  else 
    vZ <- abvZ
  if(!is.null(User.mle))
    vZ <- as.vector(User.mle)
  rms <- sqrt(mean(vZ^2))
  vZ <- vZ/rms
  if(p==0) p <- 1
  Beta <- vector( mode="numeric",length=n.sample*(p+1))
##  Beta <- vector( mode="numeric",length=n.sample*p)
  storage.mode(Beta) <- "double"
  storage.mode(vZ) <- "double"
  if(length(b.delta)==1){
    b.delta <-as.vector(matrix(b.delta, nrow=p, ncol=dimSpace), mode="numeric")
  }
  storage.mode(b.delta) <- "double"

  if(!require(mclust,quietly=TRUE, warn.conflicts = FALSE)){
   stop("You need the 'mclust' package to fit latent cluster models.")
  }
  if(dimSpace > 1){
   el.hc <- hc(modelName="VII",data=matrix(vZ,ncol=dimSpace))
  }else{
   el.hc <- hc(modelName="V",data=vZ)
  }
  cl <- hclass(el.hc,ng)
  if(dimSpace > 1){
    el.me <- me(modelName="VII",data=matrix(vZ,ncol=dimSpace),z=unmap(cl))
  }else{
    el.me <- me(modelName="V",data=vZ,z=unmap(cl))
  }
  if(any(is.na(el.me$parameters$mean)))
    {
     if(dimSpace > 1){
      el.me <- mstep(modelName="VII",data=matrix(vZ,ncol=dimSpace),z=unmap(cl))
     }else{
      el.me <- mstep(modelName="V",data=vZ,z=unmap(cl))
     }
     el.me$z <- unmap(cl)
    }
  mu <- t(el.me$parameters$mean)
  if(dimSpace > 1){
   Sigma <- el.me$parameters$variance$sigma[1,1,]
  }else{
   Sigma <- el.me$parameters$variance$sigma
  }
  Ki <- map(el.me$z)
#  cat("mu:\n")
#  print(mu)
#  cat("Sigma:\n")
#  print(Sigma)
#  cat("Ki:\n")
#  print(Ki)
  epsilon <- table(map(el.me$z))/sum(table(map(el.me$z)))
    
   # put the covariates into a vector
  vX <- as.numeric(unlist(X))
  storage.mode(vX) <- "double"

    # call the C code with the proper arguements
  cat("Calling MCMC_latent\n")
#
#       The standard non-parallel version
#
    ret <- .C("MCMC_latent2_wrapper",  
              as.integer(heads),              as.integer(tails),    
              as.integer(n.edges),            as.integer(gY$gal$n),
              as.integer(MCMCSampleSize),     as.integer(burnin),
              as.integer(interval),           as.integer(dimSpace),
              as.double(z.prior.mu),          as.double(z.prior.sd),
              as.double(b.prior.mu),          as.double(b.prior.sd),
              as.double(vZ),                  as.double(z.delta),
              as.double(vZ.post),             as.double(Z.rate),              
              ##            as.double(beta.mle),            as.double(b.delta),
              as.double(c(beta.mle,1)),            as.double(b.delta),
              as.double(Beta),                as.double(B.rate),
              as.integer(p),                  as.integer(dir),
              as.double(vX), 
              as.double(Llike),#after this, everything needs to be added...
              as.integer(ng),                 as.double(epsilon),
              as.double(mu),                  as.double(Sigma),
              as.integer(Ki),                 KiList =as.integer(KiList),
              muList = as.double(muList),     SigmaList = as.double(SigmaList),
              as.double(Sigprior),            as.double(muSigprior),
              as.double(dirprior),            as.double(alphaprior),
              as.integer(chisqprop),          as.double(thetaprop),
              PACKAGE="latentnet")
  
   # Take the positions out of vector form and put into a list of 
   # matrices
#
#  MSH: keep it simple
#

    Zp <- array(ret[[15]],dim=c(g,dimSpace,n.sample))

  # Return all the results of the chain
    Results<-list(Z.mle=Z.mle,beta.mle=beta.mle,
                  Llik=ret[[24]],
                  Beta=ret[[19]],
                  Beta.rate=ret[[20]],ngroups = ng,
                  Ki = matrix(ret$KiList,ncol=g),
                  mu = matrix(ret$muList,ncol=ng * dimSpace),
                  Sigma = matrix(ret$SigmaList,ncol=ng),
                  Z=Zp,Z.rate=ret[[18]], newnetwork=gY, mle.like = MLE.like,
                  mu.mle = mu, Sigma.mle = Sigma, Ki.mle = Ki
                  )
  Results
}

############################################################################
# Distance Matrix Function for networks calculations using sociomatrix
# for both directed and undirected i think that this one is fast
# it doesnt allow of infinite distances though.
DistMatrix<-function(X)
{
  # D is the end distance matrix
  # Tp is the adjacency martix raised to the pth power
  # X is the adjacency matrix
  D<-Tp<-X
  g<-dim(X)[1]

  # Calculate the distances.
  #  The adjacency matrix to the pth power contains in the i,j box
  #  the number of paths of length p conecting i to j.
  #  Thus we raise X to a maximum of g-1 power and record the lowest
  #  value of p that gives a number in each of the positions.
  for( p in 2:(g-1) )
  {
     # raise X one more power than before and store in Tp
     Tp<-Tp%*%X
     # if there was a 0 in the position and there is a nonzero
     #  entry in Tp store p in that position
     D<- D + p*( (D==0) & (Tp!=0) )
  }
  # make the distance between two nodes which are not connected g
  D<-D*(D>0) + 5*g*(D==0)
  # make sure diagnols are undefined.  (cycles undefined.)
  diag(D)<-0
  D 
}
##############################################################################

##############################################################################
# MINUS log prob of network, with Z in vector form, 
# With Z[1]=alpha and the rest as Z[1,1], Z[1,2], Z[2,1]..........
# to be used in by optim or nlm
mlpY.cluster<-function(abvZ, abz.list)
{
  dp <- abz.list$dp
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # put Z in matrix form
  if(dp>0)
    Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)
  else
    Z<-matrix(abvZ,nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + (abz.list$X[[i]] * beta[i])
      }
    }
  }#end ifelse BETAZERO

  # see above for formula
  eta <- tX + lpZ
#  eta <- tX*(1 + lpZ) ##JMT 22/3/04
  # Make sure that the distance between node[i] and node[i] is zero
# diag(eta)<-0
  # Return llk
# sum( abz.list$Y*eta - log( 1+exp(eta) ) ) + nnodes*log(2)
  diagY <- row(eta)!=col(eta)
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  sum( Y*eta - log( 1+exp(eta) ) )
}

mlpY.cluster.grad<-function(abvZ, abz.list)
{   

  dp <- abz.list$dp
  nnodes <- abz.list$nnodes
  dimSpace <- abz.list$dimSpace
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # put Z in matrix form
  if(dp>0)
    Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)
  else
    Z<-matrix( abvZ,nrow=nnodes,ncol=dimSpace)

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + (abz.list$X[[i]] * beta[i])
      }
    }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0

  Dlliknuij <-  abz.list$Y - exp(eta)/(1+exp(eta))

  Dllik <- vector(length=dp+nnodes*dimSpace)

  cd <- 0
  if(dp>0)
    for(i in 1:dp){
      cd <- cd + 1
      M <- abz.list$X[[i]]
      diag(M)<-0
      Dllik[cd] <- sum(Dlliknuij * M)
    }
  for(j in 1:dimSpace){
   for(i in 1:nnodes){
    v <- (Z[i,j]-Z[,j])/lpZ[i,]
    v[is.infinite(v)|is.na(v)] <- 0
    M <- matrix(0,ncol=nnodes,nrow=nnodes)
    M[i,] <- v
    M[,i] <- v
    diag(M)<-0
    cd <- cd + 1
    Dllik[cd] <- sum(Dlliknuij * M)
   }
  }

  # Return the gradient of the llk
  Dllik
}

mlpY.clustermdsZ<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + abz.list$X[[i]] * beta[i]
      }
    }
  }

  # see above for formula
  eta <- tX + lpZ
#  eta <- tX*(1 + lpZ)  ##JMT 22/3/04
  # Make sure that the distance between node[i] and node[i] is zero
# diag(eta)<-0
  # Return llk
# sum( abz.list$Y*eta - log( 1+exp(eta) ) ) + nnodes*log(2)
  diagY <- row(eta)!=col(eta)
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  sum( Y*eta - log( 1+exp(eta) ) )
}

mlpY.clustermdsZ.grad<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  Z <- abz.list$Z
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # This returns the negative distance between all the nodes
  lpZ<-lpz.dist(Z)

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + (abz.list$X[[i]] * beta[i])
      }
    }
  }

  # see above for formula
  eta <- tX + lpZ
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0

  Dlliknuij <-  abz.list$Y - exp(eta)/(1+exp(eta))

  Dllik <- vector(length=length(beta))

  cd <- 0
  if(dp>0)
    for(i in 1:dp){
      cd <- cd + 1
      M <- abz.list$X[[i]]
      diag(M)<-0
      Dllik[cd] <- sum(Dlliknuij * M)
    }

  # Return the gradient of the llk
  Dllik
}

mlpY.clusternull<-function(abvZ,Y,dp,X,nnodes,k=2)
{   
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # put Z in matrix form
  if(dp>0)
    Z<-matrix( abvZ[-(1:dp)],nrow=nnodes,ncol=k)
  else
    Z<-matrix( abvZ,nrow=nnodes,ncol=k)

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + (X[[i]] * beta[i])
      }
    }
  }

  # see above for formula
  eta <- tX
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0
  # Return llk
  -(  sum( Y*eta - log( 1+exp(eta) ) ) + dim(Y)[1]*log(2)  )
}

mlpY.cluster0<-function(abvZ,abz.list)
{   
  dp <- abz.list$dp
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + abz.list$X[[i]] * beta[i]
      }
    }
  }

  # see above for formula
  eta <- tX
# diag(eta)<-0
  # Return llk
# -(  sum( Y*eta - log( 1+exp(eta) ) ) + dim(Y)[1]*log(2)  )
  diagY <- row(eta)!=col(eta)
  Y <- abz.list$Y[diagY]
  eta <- eta[diagY]
  sum( Y*eta - log( 1+exp(eta) ) )
}

mlpY.cluster0.grad<-function(abvZ, abz.list)
{   
  dp <- abz.list$dp
  dimSpace <- abz.list$dimSpace
  nnodes <- abz.list$nnodes
  # extract parameters
  # extract parameters
  if(dp>0)
    beta <- abvZ[1:dp]

  Deta<-matrix(0, nrow=nnodes*nnodes, ncol=dp)
  cd <- 0
  if(dp>0)
    for(i in 1:dp){
      cd <- cd + 1
      M <- abz.list$X[[i]]
      diag(M)<-0
      Deta[,cd] <- M
    }

  # find the covariate contribution
  if(dp==0)
    tX <- 0
  else{
    tX <- abz.list$X[[1]] * beta[1]
    if(dp>1){
      for(i in 2:dp){
        tX <- tX + abz.list$X[[i]] * beta[i]
      }
    }
  }

  # see above for formula
  eta <- tX
  # Make sure that the distance between node[i] and node[i] is zero
  diag(eta)<-0
  # Return the gradient of the llk
  aaa <-  matrix(abz.list$Y - exp(eta)/(1+exp(eta)), nrow=1)
  aaa <- as.vector(aaa %*% Deta)
  aaa
}
