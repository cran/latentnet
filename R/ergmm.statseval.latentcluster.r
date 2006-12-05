ergmm.statseval.latentcluster <- function (z, Clist, m, MCMCsamplesize, burnin, 
            interval, formula, 
            X,dimSpace, maxit, pmodes=TRUE, penalty.sigma=c(10,0.5),
            redo.mle=FALSE, verbose=FALSE,iter.max=10) 
{
  ngroups <- z$ng
  samplesize <- dim(z$Z)[3]
  ndim <- dim(z$Z)[2]
  Nnodes <- dim(z$Z)[1]

# For procOPA
  if(!require(shapes,quietly=TRUE)){
   stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
  }

  Z.mean <- matrix(0,nrow=samplesize,ncol=ndim)
  for(i in 1:samplesize)
    Z.mean[i,] <- apply(z$Z[,,i],2,mean)

  Z.proc <- z$Z
  for(i in 1:samplesize)
    Z.proc[,,i] <- Z.proc[,,i] - matrix(Z.mean[i,],Nnodes,ndim,byrow=TRUE)

  Mu.proc <- array(0,c(ngroups,ndim,samplesize))
  for(i in 1:samplesize)
  {
    temp <- procOPA(z$Z.mle,z$Z[,,i],FALSE,TRUE)
    Z.proc[,,i] <- fcnt(z$Z[,,i]) %*% temp$R
    Mu.proc[,,i] <- fcnt(matrix(z$mu[i,],ngroups,ndim,byrow=TRUE)) %*% temp$R
  }

  
  vnames <- m$coef.names
  if(!is.null(z$Z)){
   l <- list(sample=NA, iterations=NA,
             MCMCtheta = z$beta.mle, 
             loglikelihood=max(c(z$Llik,z$mle.lik)),
             mcmc.loglikelihood=z$Llik,
             gradient = NA)
   l$Beta <- matrix(z$Beta,nrow=samplesize)
   if(verbose){trace <- 4}else{trace <- 0}
#
#  Posterior means
#
   Z.pm <- apply(Z.proc,c(1,2),mean)
#
   if(pmodes){
#
#   MSH: Add posterior modes
#   Warning: Only does two-dimensional stuff
#
#    require(KernSmooth,quietly=TRUE)
#    pmode <- function(x){
#     est <- bkde2D(x=x[,1:2], bandwidth=c(2,2))
#     dmax <- order(-est$fhat)[1]
#     c(est$x1[row(est$fhat)[dmax]],est$x2[col(est$fhat)[dmax]])
#    }
##   Z.pmode <- array(0, dim=dim(z$Z)[1:2])
#    Z.pmode <- array(0, dim=c(dim(z$Z)[1],2))
#
    pmode <- function(x){
      mvimode(x)$theta
    }
    Z.pmode <- array(0, dim=c(Nnodes,ndim))
    for(k in 1:Nnodes){
     Z.pmode[k,] <- pmode(t(Z.proc[k,,]))
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

   eta <- Xm %*% t(l$Beta[,-(dp+1)])
   for(i in (1:samplesize)){
    eta[,i] <- eta[,i] + as.vector(lpz.dist(Z.proc[,,i]))
#    eta[,i] <- eta[,i] + l$Beta[i,(dp+1)] *as.vector(lpz.dist(Z.proc[,,i]))
   }
   Z.pp <- exp(eta)/(1+exp(eta))
   Z.pp[eta > 700] <- 1
   Z.pp <- apply(Z.pp,1,mean)
   pY <- matrix(Z.pp, ncol=nnodes)
   diag(pY) <- 0
   if(dp>1)
     abvZ <- c(apply(l$Beta[,-(dp+1)],2,mean),Z.pmode)
   else
     abvZ <- c(mean(l$Beta[,-(dp+1)]),Z.pmode)
   if(penalty.sigma[1]>0){
     penalty.factor <- c(1/(penalty.sigma[1]*penalty.sigma[1]),penalty.sigma[2])
   }else{
     penalty.factor <- c(0,penalty.sigma[2])
   }
   abz.list <- list(Y=pY,dp=dp,X=X,nnodes=nnodes,dimSpace=dimSpace,
                    penalty.factor=penalty.factor,
                    reach=reach,directed=is.directed(z$newnetwork))
   cat("Calling min KL fit to",dimSpace,"dimensions.\n")
   MKL.fit <- try(optim(par=abvZ,fn=mlpY.cluster, gr=mlpY.cluster.grad,
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
#  Use KL to seed MLE
#
     abvZ <- c(l$beta.mkl,l$Z.mkl)
     abz.list <- list(Y=Y,dp=dp,X=X,nnodes=nnodes,dimSpace=dimSpace,
                      penalty.factor=penalty.factor,
                      reach=reach,directed=is.directed(z$newnetwork))
     cat("Calling true MLE fit\n")
     MLE.fit <- try(optim(par=abvZ,fn=mlpY.cluster,gr=mlpY.cluster.grad,
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
     l$MCMCtheta <- abvZ[(1:dp)]
     Z.mle <- matrix(abvZ[-(1:dp)],nrow=nnodes,ncol=dimSpace)
   }else{
     l$mle.lik <- z$mle.like
     Z.mle <- z$Z.mle
   }
   l$Z <- Z.proc
 }else{
   l <- list(sample=NA, iterations=z$Beta.rate,
             MCMCtheta = z$beta.mle, 
             loglikelihood=max(c(z$Llik,z$mle.lik)),
             mcmc.loglikelihood=z$Llik,
             gradient = z$Z.rate)
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
  l$samplesize <- samplesize
  l$Z.mle <- Z.mle
  l$Z.pmean <- Z.pm
  l$Z.pmode <- Z.pmode
  l$latent <- TRUE
  l$cluster <- TRUE
  l$newnetwork <- z$newnetwork
  l$formula <- formula

  ####################################################################
  if(ngroups>1)
  {
    mu.0 <- apply(z$Z.mle,2,function(x,ki)tapply(x,ki,mean),ki = z$Ki.mle)

    permute <- ergmm.permutation(ngroups)

    vt.c <- matrix(0,samplesize,ngroups)
    Z.mu <- array(0,c(ngroups,ndim,samplesize))
    Z.Ki <- z$Ki
    Z.sigma <- matrix(NA,ngroups,samplesize)
    for(loop in 1:samplesize)
    {
      mu.1 <- apply(Z.proc[,,loop],2,function(x,ki)tapply(x,ki,mean),ki = Z.Ki[loop,])
      mutab <- table(Z.Ki[loop,])
      mu.names <- as.numeric(names(mutab))
      n1 <- length(mutab)
      d1 <- as.matrix(dist(rbind(mu.1,mu.0)))[1:n1,(n1+1):(n1+ngroups)]
      if(length(mutab)==1)
      {
        d1.min <- order(d1)[1]
        Z.Ki[loop,] <- d1.min
        Z.mu[d1.min,,loop] <- mu.1
      }
      else
      {
        d1.use <- rep(TRUE,ngroups)
        if(length(mutab)<ngroups)
        {
          d1.use <- rep(FALSE,ngroups)
          d1.new <- matrix(0,ngroups,ngroups)
          mu.new <- matrix(0,ngroups,ndim)
          j <- 1
          for(i in 1:ngroups)
            if(any(mu.names == i))
            {
              d1.new[i,] <- d1[j,]
              mu.new[i,] <- mu.1[j,]
              j <- j + 1
              d1.use[i] <- TRUE
            }
          d1 <- d1.new
          mu.1 <- mu.new
        }
    
        d1.vec <- rep(0,nrow(permute))
        for(j in 1:nrow(permute))
          for(i in 1:ngroups)
            d1.vec[j] <- d1.vec[j] + d1.use[i] * d1[i,permute[j,i]]
    
        d1.min <- order(d1.vec)[1]
        per.to <- order(permute[d1.min,])
        vt.c[loop,] <- per.to
        Z.Ki[loop,] <-ergmm.labelswitch(z$Ki[loop,],per.to)
        Z.mu[,,loop] <- mu.1[per.to,]
        Z.sigma[,loop] <- z$Sigma[loop,per.to]
      }
    }

    vt <- vt.c

    sum.diff <- 1e5
    iter <- 0
    minat.old <- minat <- rep(0,samplesize)

    while((sum.diff>10) & (iter<iter.max))
    {
      qig <- matrix(0,Nnodes,ngroups)
      for(i in 1:Nnodes)
        for(g in 1:ngroups)
        {
          qig[i,g] <- 0
          for(k in 1:samplesize)
            qig[i,g] <- qig[i,g] + 1/samplesize * prod(dnorm(Z.proc[i,,k],Mu.proc[vt[k,g],,k],z$Sigma[k,vt[k,g]]))
        }
      qig <- qig / apply(qig,1,sum)

      minat <- klswitch.c(qig,permute,Z.proc,Mu.proc,z$Sigma)$minat + 1
      vt <- permute[minat,]
      sum.diff <- sum(minat!=minat.old)
      minat.old <- minat
      iter <- iter + 1
    }
    vt <- permute[minat,]
    qig <- matrix(0,Nnodes,ngroups)
    pnames <- network.vertex.names(z$newnetwork)
    dimnames(qig) <- list(pnames, 1:ngroups)
    for(i in 1:Nnodes)
      for(g in 1:ngroups)
        {
          qig[i,g] <- 0
          for(k in 1:samplesize)
            qig[i,g] <- qig[i,g] + 1/samplesize * prod(dnorm(Z.proc[i,,k],Mu.proc[vt[k,g],,k],z$Sigma[k,vt[k,g]]))
        }
    qig <- qig / apply(qig,1,sum)
  
    labs <- apply(qig,1,function(x)order(x)[ngroups])

    Z.Ki <- matrix(0,dim(Z.proc)[3],dim(Z.proc)[1])
    for(i in 1:nrow(Z.Ki))
      Z.Ki[i,] <- ergmm.labelswitch(z$Ki[i,],vt[i,])
  } else {
    qig <- matrix(1,Nnodes,ngroups)
    vt <- matrix(1,samplesize,ngroups)
    labs <- rep(1,Nnodes)
    Z.Ki <- z$Ki
  }
  
  logit.negloglike <- function(beta,Y,Y.dist)
    {
      eta <- beta - Y.dist
      return(-sum(log(exp(Y * eta)/(1+ exp(eta)))))
    }

  distmat <- as.matrix(dist(l$Z.mkl))
  Y <- as.sociomatrix(z$newnetwork)[!diag(Nnodes)]
  Y.dist <- distmat[!diag(Nnodes)]

  abvZ <- c(0.36)
  logit.fit <- optim(par=abvZ,fn=logit.negloglike,method="BFGS", hessian=TRUE,
                     control=list(maxit=200,trace=trace),
                     Y=Y,Y.dist=Y.dist)

# Next original
# bicLR <- -2 * logit.fit$value - 1 * log(Nnodes*(Nnodes-1))
  dlogit <- dp
  bicLR <- -2 * logit.fit$value - dlogit * log(sum(Y))
  aicLR <- -2 * logit.fit$value - dlogit * 2
  labs.use <- labs
  Z.mkl.use <- l$Z.mkl
  
  temp <-   as.numeric(names(table(labs.use))[table(labs.use) < 3])
  for(i in temp)
    {
      Z.mkl.use <- Z.mkl.use[!labs.use==i,]
      labs.use <- labs.use[!labs.use==i]
    }

  mbc.fit <- me(modelName="VII",Z.mkl.use,unmap(labs.use))
  bicMBC <- bic("VII",mbc.fit$loglik,mbc.fit$n,mbc.fit$d,ngroups)

  l$d.mbc <- (ngroups) * 4 - 1  #d
  l$ngroups <- ngroups  #G
  l$logl.lr <- -logit.fit$value  #lllr
  l$logl.mbc <- mbc.fit$loglik  #llmbc
  aicMBC <- 2*l$logl.mbc - l$d.mbc*2

  l$aic <- aicMBC + aicLR  #AIC
  l$BIC <- bicMBC + bicLR  #BIC
  l$bic <- l$BIC
  l$logl <- -logit.fit$value + 2*mbc.fit$loglik  #llik

#
# Two-stage MLE
#
# mbc.fit.to.mle <- me(modelName="VII",Z.mkl.use,unmap(labs.use))
# bicMBC.to.mle <- bic("VII",mbc.fit.to.mle$loglik,mbc.fit.to.mle$n,mbc.fit.to.mle$d,ngroups)
# l$logl.mbc <- mbc.fit.to.mle$loglik  #llmbc
# aicMBC.to.mle <- 2*l$logl.mbc.to.mle - l$d.mbc*2
# l$mle2.aic <- aicMBC.to.mle + aicLR  #AIC
# l$mle2.BIC <- bicMBC.to.mle + bicLR  #BIC

  l$class <- labs
  l$Ki <- Z.Ki
  l$Ki.mle <- z$Ki.mle
  l$qig <- qig
  Mu.new <- Mu.proc
  for(g in 1:ngroups)
    for(k in 1:samplesize)
      Mu.new[g,,k] <- Mu.proc[vt[k,g],,k]
  Sigma.new <- z$Sigma
  for(k in 1:samplesize)
    Sigma.new[k,] <- z$Sigma[k,vt[k,]]
  
  
  l$mu <- Mu.new
  l$mu.mle <- z$mu.mle
  l$Sigma.mle <- z$Sigma.mle
  l$Sigma <- Sigma.new
  ####################################################################
  
  structure(l, class = "ergmm")
}
