latent.wrapper <- function(theta0, trms, g, m, Clist, mClist,
                           MCMCsamplesize=1000, burnin=0, interval=1,
                           formula,
                           latent.control, verbose, ...
                          ){
#
      P <- 0
      for(i in 1:(length(attr(trms,"variables"))-2)){
       if(m$options[[i]]$name =="latentcov") {
	 P <- P + 1
       }
       if(m$options[[i]]$name =="latent"){
	 uo <- m$options[[i]]$inputs
       }
       if(m$options[[i]]$name =="latentcluster"){
         uo <- m$options[[i]]$inputs
       }
      }
      plm <- m
      killlatentcov <- rep(FALSE,length(plm$options))
      addintercept <- TRUE
      for(i in 1:length(killlatentcov)){
       if(plm$options[[i]]$name =="latentcov") {
	 killlatentcov[i] <- TRUE
       }
       if(plm$options[[i]]$name =="latent") {
	 killlatentcov[i] <- TRUE
       }
       if(plm$options[[i]]$name =="latentcluster") {
         killlatentcov[i] <- TRUE
       }
       if(plm$options[[i]]$name =="-1") {
	 addintercept <- FALSE
       }
      }
      plm$options[killlatentcov] <- NULL
      killlatentcov <- rep(FALSE,length(plm$coef.names))
      for(i in seq(along=plm$coef.names)){
       killlatentcov[i] <- pmatch("latentcov.",plm$coef.names[i],nomatch=0)==1
      }
      plm$coef.names <- plm$coef.names[!killlatentcov]
      if(!is.null(plm$coef.names)){
        pl.info <- ergmm.plinfo.latent(ergmm.Cprepare.latent(g, plm), mClist, plm)$xmat
        P <- P + ncol(pl.info) + addintercept
      }else{
        pl.info <- matrix(0,ncol=0,nrow=1)
        P <- P + addintercept
      }
      X <- vector(mode="list",length=P)
      p <- 0
      nnodes <- network.size(g)
      temp <- matrix(0,ncol=nnodes,nrow=nnodes)
      base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
      base <- base[base[, 2] > base[, 1], ]
      if (is.directed(g)) {
	base <- rbind(base, cbind(base[, 2], base[, 1]))
      }
      if(addintercept){
	  p <- p + 1
	  temp[row(temp)!=col(temp)] <- 1
	  X[[p]] <- temp
	  plm$coef.names <- c("density",plm$coef.names)
	  m$coef.names <- c("density",m$coef.names)
      }
      if(ncol(pl.info)>0){
       for(i in 1:ncol(pl.info)){
	  p <- p + 1
	  temp[base] <- pl.info[,i]
	  X[[p]] <- temp
       }
      }
      for(i in 1:(length(attr(trms,"variables"))-2))
	{
	  if(m$options[[i]]$name =="latentcov")
	  {
	    p <- p + 1
            xmat <- m$options[[i]]$inputs
	    X[[p]] <-  matrix(xmat[-c(1:4)],ncol=xmat[4])
	  }
	}
     
      if(is.null(latent.control$MLEonly)){latent.control$MLEonly <- FALSE}
      if(is.null(latent.control$maxit)){latent.control$maxit <- 40}
      if(is.null(latent.control$penalty.sigma)){
         latent.control$penalty.sigma <- c(10,0.5)}

      if(!is.latent.cluster(m))
      {
        z <- ergmm.latent(gY=g,dimSpace=uo[4],p=p,X=X,
                         theta0=theta0,
                         MCMCSampleSize=MCMCsamplesize,
                         burnin=burnin, interval=interval,
                         z.delta=uo[5], z.prior.mu=uo[6], z.prior.sd=uo[7],
                         b.delta=uo[8], b.prior.mu=uo[9], b.prior.sd=uo[10],
                         maxit=latent.control$maxit,
                         penalty.sigma=latent.control$penalty.sigma,
                         MLEonly=latent.control$MLEonly,
                         verbose=verbose, ...)
        ###NEED to organise output of the function.
        v <- ergmm.statseval.latent(z, Clist, m, MCMCsamplesize, 
	    burnin, interval, formula, 
	    X,dimSpace=uo[4],
	    maxit=latent.control$maxit,
	    penalty.sigma=latent.control$penalty.sigma)
      }else{
        #do latentcluster
        z <- ergmm.latentcluster(gY=g,dimSpace=uo[4],ng=uo[5],p=p,X=X,
                         theta0=theta0,
                         MCMCSampleSize=MCMCsamplesize,
                         burnin=burnin, interval=interval,
                         z.delta=uo[6], z.prior.mu=uo[7], z.prior.sd=uo[8],
                         b.delta=uo[9], b.prior.mu=uo[10], b.prior.sd=uo[11],
                         Sigprior=uo[12], muSigprior = uo[13],dirprior =uo[14],
                         alphaprior = uo[15], chisqprop = uo[16],
                         thetaprop = uo[17], maxit=latent.control$maxit, ...)
#                         penalty.sigma=latent.control$penalty.sigma,
#                         verbose=verbose, ...)
        ###NEED to organise output of the function.
        if(uo[4] > 1){
         v <- ergmm.statseval.latentcluster(z, Clist, m, MCMCsamplesize, 
	    burnin, interval, formula, 
	    X,dimSpace=uo[4],
	    maxit=latent.control$maxit,
	    penalty.sigma=latent.control$penalty.sigma)
        }else{
         v <- ergmm.statseval.latent1cluster(z, Clist, m, MCMCsamplesize, 
	    burnin, interval, formula, 
	    X,dimSpace=uo[4],
	    maxit=latent.control$maxit,
	    penalty.sigma=latent.control$penalty.sigma)
        }
        v$logl <- NULL
      }#end if cluster
      yij <- as.matrix.network(g, matrix.type="adjacency")
      X.l <- X
      form <- "yij ~ -1 + "
      for(i in 1:p)
        {
          X.l[[i]] <- X[[i]][row(yij) != col(yij)]
          if(i > 1)
            form <- paste(form," + X.l[[",i,"]]",sep="")
          else form <- paste(form,"X.l[[",i,"]]",sep="")
        }
      form <- formula(form)
      yij <- yij[row(yij) != col(yij)]
      glm.out <- glm(form,family="binomial")
      aaa <- t(v$Z[1,,])
      if(dim(v$Z)[2]==1){aaa <- matrix(aaa,ncol=1)}
      if(is.latent.cluster(m)){
        statsmatrix <- cbind(v$mcmc.loglikelihood[1:nrow(v$Beta)],
                             v$Beta[,-ncol(v$Beta)], aaa)
        colnames(statsmatrix) <- c("mcmc.loglikelihood", 
                                   m$coef.names,
#                                 "beta1",
                                   paste("Z",1:dim(v$Z)[2])
                                  )
      }else{
        statsmatrix <- cbind(v$mcmc.loglikelihood[1:nrow(v$Beta)],
                             v$Beta, aaa)
        colnames(statsmatrix) <- c("mcmc.loglikelihood", 
                                   m$coef.names,
                                   paste("Z",1:dim(v$Z)[2])
                                  )
      }
      endrun <- burnin+interval*(MCMCsamplesize-1)
      attr(statsmatrix, "mcpar") <- c(burnin+1, endrun, interval)
      attr(statsmatrix, "class") <- "mcmc"
      v$sample <- statsmatrix
#     v$mc.se <- NA 
      v$mc.se <- NULL 
      v$gradient <- NULL 
      v$MCMCtheta <- NULL 
      v$network <- g
#     v$newnetwork <- z$newnetwork
      v$network <- z$newnetwork
      v$newnetwork <- NULL
#     v$glm <- glm.out
      v$glm <- NULL
      v$glm.names <- NULL
#     v$samplesize <- NULL
      v$coef.names <- m$coef.names
#     v$null.deviance <- v$glm$null.deviance
      v$null.deviance <- glm.out$null.deviance
#     v$mle.lik <- -0.5*v$glm$deviance
#     Nullify some clustering stuff
      v$d.mbc <- NULL
      if(!is.latent.cluster(m)){
        v$cluster <- FALSE
#       v$bic <- NULL
#       v$aic <- NULL
      }else{
        v$cluster <- TRUE
      }
      v
}
