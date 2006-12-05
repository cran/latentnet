ergmm.statseval.latent <- function (z, Clist, m, MCMCsamplesize, burnin, 
            interval, formula, 
            X,dimSpace, maxit, pmodes=TRUE, penalty.sigma=c(10,0.5),
            redo.mle=FALSE, verbose=TRUE) 
{
  vnames <- m$coef.names
  if(!is.null(z$Z)){
#
#  So MCMC values
#
   l <- list(sample=NA, iterations=z$mle.iterations,
#            MCMCtheta = z$beta.mle, 
             loglikelihood=max(c(z$Llik,z$mle.lik)),
             mcmc.loglikelihood=z$Llik)
#            gradient = NA)
   Nnodes <- dim(z$Z)[1]
   ndim <- dim(z$Z)[2]
   samplesize <- dim(z$Z)[3]
   l$Beta <- matrix(z$Beta,nrow=samplesize)
   if(verbose){trace <- 4}else{trace <- 0}
#
#  Posterior means
#
   Z.pm <- apply(z$Z,c(1,2),mean)
#
   if(pmodes){
#
#   MSH: Add posterior modes
#   Warning: Only does two-dimensional stuff
#
#   require(KernSmooth,quietly=TRUE)
#   pmode <- function(x){
#    est <- bkde2D(x=x[,1:2], bandwidth=c(2,2))
#    dmax <- order(-est$fhat)[1]
#    c(est$x1[row(est$fhat)[dmax]],est$x2[col(est$fhat)[dmax]])
#   }
##  Z.pmode <- array(0, dim=dim(z$Z)[1:2])
#   Z.pmode <- array(0, dim=c(Nnodes, ndim))
#
    pmode <- function(x){
      mvimode(x)$theta
    }
    Z.pmode <- array(0, dim=c(Nnodes,ndim))
    for(k in 1:Nnodes){
     aaa <- t(z$Z[k,,])
     if(dim(z$Z)[2]==1){aaa <- t(aaa)}
     Z.pmode[k,] <- pmode(aaa)
    }
   }else{
    Z.pmode <- Z.pm
   }
#
#  Next reports posterior modes as primary
#
#  l$coef <- apply(l$Beta,2,mean)
#  vcov <- var(l$Beta)
#  vcov <- -solve(vcov)
#  dimnames(vcov) <- list(vnames,vnames)
#
#  Calculate the correct MLE
#
   Y <- as.sociomatrix(z$newnetwork)
   reach <- ergmm.geodesicmatrix(z$newnetwork)!=Inf
   nnodes <- dim(z$Z)[1]
   dp <- length(X)
#  
#  Z minimizing the mean posterior KL distance
#
   Xm <- matrix(unlist(X),ncol=length(X))
   eta <- Xm %*% t(l$Beta)
   for(i in (1:samplesize)){
    eta[,i] <- eta[,i] + as.vector(lpz.dist(z$Z[,,i]))
   }
   Z.pp <- exp(eta)/(1+exp(eta))
   Z.pp[eta > 700] <- 1
   Z.pp <- apply(Z.pp,1,mean)
   pY <- matrix(Z.pp, ncol=nnodes)
   diag(pY) <- 0
   abvZ <- c(apply(l$Beta,2,mean),Z.pmode)
   if(penalty.sigma[1]>0){
     penalty.factor <- c(1/(penalty.sigma[1]*penalty.sigma[1]),penalty.sigma[2])
   }else{
     penalty.factor <- c(0,penalty.sigma[2])
   }
   abz.list <- list(Y=pY,dp=dp,X=X,nnodes=nnodes,dimSpace=dimSpace,
                    penalty.factor=penalty.factor,
                    reach=reach,directed=is.directed(z$newnetwork))
   cat("Calling min KL fit to",dimSpace,"dimensions.\n")
   MKL.fit <- try(optim(par=abvZ,fn=mlpY, gr=mlpY.grad,
                 method="BFGS", hessian=TRUE,
                 control=list(fnscale=-1, maxit=maxit,trace=trace),
                 abz.list=abz.list))
   if(inherits(MKL.fit,"try-error")){
    warning("MKL could not be found.")
    MKL.like <- NA
    vcov <- diag(abz.list$dp)
   }else{
    abvZ <- MKL.fit$par
    MKL.like <- MKL.fit$value
    vcov <- MKL.fit$hessian[1:abz.list$dp, 1:abz.list$dp]
    vcov <- matrix(vcov,ncol=abz.list$dp)
    dimnames(vcov) <- list(vnames,vnames)
   }
#
   l$Z.mkl <- matrix( abvZ[-(1:abz.list$dp)],nrow=nnodes,ncol=dimSpace)
   l$beta.mkl <- abvZ[1:abz.list$dp]
   l$coef <- l$beta.mkl
   if(redo.mle){
#
#   Use KL to seed MLE
#
    abvZ <- c(l$beta.mkl,l$Z.mkl)
    abz.list <- list(Y=Y,dp=dp,X=X,nnodes=nnodes,dimSpace=dimSpace,
                     penalty.factor=penalty.factor,
                     reach=reach,directed=is.directed(z$newnetwork))
    cat("Calling true MLE fit\n")
    MLE.fit <- try(optim(par=abvZ,fn=mlpY,gr=mlpY.grad,
                  method="BFGS",
                  control=list(fnscale=-1, maxit=maxit,trace=trace),
                  abz.list=abz.list))
    if(inherits(MLE.fit,"try-error")){
     warning("MLE could not be found.")
     MLE.like <- NA
     l$mle.lik <- MKL.like
     vcov <- diag(abz.list$dp)
    }else{
     abvZ <- MLE.fit$par
     MLE.like <- MLE.fit$value
     l$mle.lik <- MLE.fit$value
    }
#   l$MCMCtheta <- abvZ[(1:dp)]
    Z.mle <- matrix(abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)
   }else{
    l$mle.lik <- z$mle.like
    Z.mle <- z$Z.mle
   }
   l$Z <- z$Z
  }else{
#
#  So no MCMC and use MLE fits
#
   l <- list(sample=NA, iterations=z$mle.iterations,
#            MCMCtheta = z$beta.mle, 
             loglikelihood=max(c(z$Llik,z$mle.lik)),
             mcmc.loglikelihood=z$Llik)
#            gradient = z$Z.rate)
   Z.pm <- NULL
   Z.pmode <- NULL
   Z.mle <- NULL
   samplesize <- 0
   l$Beta <- z$Beta.mle
   l$coef <- z$Beta
   vcov <- z$Hessian
   dimnames(vcov) <- list(vnames,vnames)
  }
#
# Next for PM version (cut later!)
# l$MCMCtheta <- apply(l$Beta,2,mean)
#
  l$hessian <- vcov
  names(l$coef) <- vnames
# Z.mle <- z$Z.mle
# Z.mle <- (mean(z$Alpha)/z$alpha.mle)*(z$Z.mle)
#
  dlogit <- dp
  bicLR <- 2 * l$loglikelihood - dlogit * log(sum(Y))
  aicLR <- 2 * l$loglikelihood - dlogit * 2
#
  l$aic <- aicLR  #AIC
  l$BIC <- bicLR  #BIC
  l$bic <- l$BIC
#
  l$samplesize <- samplesize
  l$Z.mle <- Z.mle
  l$Z.pmean <- Z.pm
  l$Z.pmode <- Z.pmode
  l$latent <- TRUE
  l$newnetwork <- z$newnetwork
  l$formula <- formula
  structure(l, class = "ergmm")
}
