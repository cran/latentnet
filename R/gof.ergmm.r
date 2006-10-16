gof <- function(object, ...){
 UseMethod("gof")
}

gof.default <- function(object,...)
{
  stop("Either a ergmm object or a formula argument must be given")
}

gof.ergmm <- function (object, ..., nsim=100,
                      GOF=~degree+espartners+distance, 
		      verbose=FALSE) {

# g <- as.network(object$graph)
# if(!is.network(g)){
#   stop("A graph in the network object must be given")
# }
# n <- network.size(g)
# if(is.null(seed)){seed <- sample(10000000, size=1)}

   
  formula <- object$formula

  trms <- ergmm.getterms.latent(formula)
  if(length(trms)>2){
    nw <- eval(trms[[2]], sys.parent())
  }else{
    stop("A network object on the RHS of the formula argument must be given")
  }

  nsim <- max(nsim, dim(object$Z)[3])

  all.gof.vars <- all.vars(GOF)

# match variables

  for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree','triadcensus','model'
     )
                               )
  }
  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

  if(!is.network(nw)){
    stop("A network object on the RHS of the formula argument must be given")
  }

  pval.model<-pval.triadcensus<-pval.dist<-pval.deg<-pval.espart<-pval.espart<-NULL
#
  obs.model<-pobs.model<-sim.model<-psim.model<-pval.model<-bds.model<-NULL
  obs.triadcensus<-pobs.triadcensus<-sim.triadcensus<-psim.triadcensus<-pval.triadcensus<-bds.triadcensus<-NULL
  obs.dist<-pobs.dist<-sim.dist<-psim.dist<-pval.dist<-bds.dist<-NULL
  obs.deg<-pobs.deg<-sim.deg<-psim.deg<-pval.deg<-bds.deg<-NULL
  obs.espart<-pobs.espart<-sim.espart<-psim.espart<-pval.espart<-bds.espart<-NULL
  obs.dspart<-pobs.dspart<-sim.dspart<-psim.dspart<-pval.dspart<-bds.dspart<-NULL

  obs.ideg<-pobs.ideg<-sim.ideg<-psim.ideg<-pval.ideg<-bds.ideg<-pval.ideg<-NULL
  obs.odeg<-pobs.odeg<-sim.odeg<-psim.odeg<-pval.odeg<-bds.odeg<-pval.odeg<-NULL

  n <- network.size(nw)

  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables

  if ('model' %in% all.gof.vars) {
   obs.model <- summary(formula)
   sim.model <- array(0,dim=c(nsim,length(obs.model)))
   dimnames(sim.model) <- list(paste(c(1:nsim)),names(obs.model))
  }

  if ('distance' %in% all.gof.vars) {
   obs.dist <- ergm.geodistdist(nw)
   obs.dist[obs.dist==Inf] <- n
   sim.dist <-array(0,dim=c(nsim,n))
   dimnames(sim.dist)  <- list(paste(c(1:nsim)),paste(1:n))
  }

  if ('odegree' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.odeg <- summary(as.formula(paste('nw ~ odegree(',mesp,')',sep="")),drop=FALSE)
   sim.odeg <- array(0,dim=c(nsim,n))
   dimnames(sim.odeg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.odeg) <- dimnames(sim.odeg)[[2]]
  }

  if ('idegree' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.ideg <- summary(as.formula(paste('nw ~ idegree(',mesp,')',sep="")),drop=FALSE)
   sim.ideg <- array(0,dim=c(nsim,n))
   dimnames(sim.ideg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.ideg) <- dimnames(sim.ideg)[[2]]
  }

  if ('degree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     obs.deg <- summary(as.formula(paste('nw ~ degree(',mesp,')',sep="")),drop=FALSE)
   sim.deg <- array(0,dim=c(nsim,n))
   dimnames(sim.deg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.deg) <- dimnames(sim.deg)[[2]]
  }
 
  if ('espartners' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.espart <- summary(as.formula(paste('nw ~ esp(',mesp,')',sep="")), drop=FALSE)
   sim.espart <- array(0,dim=c(nsim,n-1))
   dimnames(sim.espart) <- list(paste(c(1:nsim)),paste(0:(n-2)))
  }
 
  if ('dspartners' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.dspart <- summary(as.formula(paste('nw ~ dsp(',mesp,')',sep="")), drop=FALSE)
   sim.dspart <- array(0,dim=c(nsim,n-1))
   dimnames(sim.dspart) <- list(paste(c(1:nsim)),paste(0:(n-2)))
  }

  if ('triadcensus' %in% all.gof.vars) {
    obs.triadcensus <- SimCond$summary.triadcensus[,"mean"]
   sim.triadcensus <- array(0,dim=c(nsim,16))
   dimnames(sim.triadcensus) <- list(paste(c(1:nsim)),
    c("003","012", "102", "021D", "021U", "021C", "111D", "111U", "030T",
      "030C", "201", "120D", "120U", "120C", "210", "300"))
  }
 
  # Simulate an exponential family random graph model

  SimNetworkSeriesObj <- rergm(object,mkl=FALSE,n=nsim)

  if(verbose){cat("\nCollating simulations\n")}

  for (i in 1:nsim)
  { 
    if(verbose){
     cat("\nCalculating statistics for simulation",i,"\n")
    }

    if ('model' %in% all.gof.vars) {
     sim.model[i,] <- summary(update(formula,SimNetworkSeriesObj$networks[[i]] ~ .))
    }

    if ('distance' %in% all.gof.vars) {
     sim.dist[i,] <- ergm.geodistdist(SimNetworkSeriesObj$networks[[i]])
    }
    if ('idegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- SimNetworkSeriesObj$networks[[i]]
     sim.ideg[i,] <- summary(as.formula(paste('gi ~ idegree(',mesp,')',sep="")),drop=FALSE)
#    temp <- table(degreedist(SimNetworkSeriesObj$networks[[i]], print=verbose)[1,])
#    sim.ideg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('odegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- SimNetworkSeriesObj$networks[[i]]
     sim.odeg[i,] <- summary(as.formula(paste('gi ~ odegree(',mesp,')',sep="")),drop=FALSE)
    }
    if ('degree' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj$networks[[i]]
     if(is.bipartite(gi)){
      temp <- degreedist(gi, print=FALSE)$event
      sim.deg[i,] <- c(temp,rep(0,n-length(temp)))
     }else{
      mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
      sim.deg[i,] <- summary(as.formula(paste('gi ~ degree(',mesp,')',sep="")),drop=FALSE)
     }
    }
    if ('espartners' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj$networks[[i]]
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.espart[i,] <- summary(as.formula(paste('gi ~ esp(',mesp,')',sep="")), drop=FALSE)
    }
    if ('dspartners' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj$networks[[i]]
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.dspart[i,] <- summary(as.formula(paste('gi ~ dsp(',mesp,')',sep="")), drop=FALSE)
    }
    if ('triadcensus' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj$networks[[i]]
     sim.triadcensus[i,] <- summary(as.formula('gi ~ triadcensus(1:16)'), drop=FALSE)
    }
  }

  # calculate p-values

 if ('model' %in% all.gof.vars) {
  pval.model <- apply(sim.model <= obs.model[col(sim.model)],2,mean)
  pval.model.top <- apply(sim.model >= obs.model[col(sim.model)],2,mean)
  pval.model <- cbind(obs.model,apply(sim.model, 2,min), apply(sim.model, 2,mean),
                apply(sim.model, 2,max), pmin(1,2*pmin(pval.model,pval.model.top)))
  dimnames(pval.model)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.model <- pval.model.top
  psim.model <- apply(sim.model,2,rank)/nrow(sim.model)
  bds.model <- apply(psim.model,2,quantile,probs=c(0.025,0.975))
 }

 if ('distance' %in% all.gof.vars) {
  pval.dist <- apply(sim.dist <= obs.dist[col(sim.dist)],2,mean)
  pval.dist.top <- apply(sim.dist >= obs.dist[col(sim.dist)],2,mean)
  pval.dist <- cbind(obs.dist,apply(sim.dist, 2,min), apply(sim.dist, 2,mean),
                apply(sim.dist, 2,max), pmin(1,2*pmin(pval.dist,pval.dist.top)))
  dimnames(pval.dist)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.dist <- obs.dist/sum(obs.dist)
  psim.dist <- sweep(sim.dist,2,apply(sim.dist,1,sum),"/")
  bds.dist <- apply(psim.dist,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for minimum geodesic distance\n\n")
# print(pval.dist)

 if ('idegree' %in% all.gof.vars) {
  pval.ideg <- apply(sim.ideg <= obs.ideg[col(sim.ideg)],2,mean)
  pval.ideg.top <- apply(sim.ideg >= obs.ideg[col(sim.ideg)],2,mean)
  pval.ideg <- cbind(obs.ideg,apply(sim.ideg, 2,min), apply(sim.ideg, 2,mean),
                apply(sim.ideg, 2,max), pmin(1,2*pmin(pval.ideg,pval.ideg.top)))
  dimnames(pval.ideg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.ideg <- obs.ideg/sum(obs.ideg)
  psim.ideg <- sweep(sim.ideg,2,apply(sim.ideg,1,sum),"/")
  bds.ideg <- apply(psim.ideg,2,quantile,probs=c(0.025,0.975))
 }

 if ('odegree' %in% all.gof.vars) {
  pval.odeg <- apply(sim.odeg <= obs.odeg[col(sim.odeg)],2,mean)
  pval.odeg.top <- apply(sim.odeg >= obs.odeg[col(sim.odeg)],2,mean)
  pval.odeg <- cbind(obs.odeg,apply(sim.odeg, 2,min), apply(sim.odeg, 2,mean),
                apply(sim.odeg, 2,max), pmin(1,2*pmin(pval.odeg,pval.odeg.top)))
  dimnames(pval.odeg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.odeg <- obs.odeg/sum(obs.odeg)
  psim.odeg <- sweep(sim.odeg,2,apply(sim.odeg,1,sum),"/")
  bds.odeg <- apply(psim.odeg,2,quantile,probs=c(0.025,0.975))
 }

 if ('degree' %in% all.gof.vars) {
  pval.deg <- apply(sim.deg <= obs.deg[col(sim.deg)],2,mean)
  pval.deg.top <- apply(sim.deg >= obs.deg[col(sim.deg)],2,mean)
  pval.deg <- cbind(obs.deg,apply(sim.deg, 2,min), apply(sim.deg, 2,mean),
                apply(sim.deg, 2,max), pmin(1,2*pmin(pval.deg,pval.deg.top)))
  dimnames(pval.deg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.deg <- obs.deg/sum(obs.deg)
  psim.deg <- sweep(sim.deg,2,apply(sim.deg,1,sum),"/")
  bds.deg <- apply(psim.deg,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for degree\n\n")
# print(pval.deg)

 if ('espartners' %in% all.gof.vars) {
  pval.espart <- apply(sim.espart <= obs.espart[col(sim.espart)],2,mean)
  pval.espart.top <- apply(sim.espart >= obs.espart[col(sim.espart)],2,mean)
  pval.espart <- cbind(obs.espart,apply(sim.espart, 2,min), apply(sim.espart, 2,mean),
                apply(sim.espart, 2,max), pmin(1,2*pmin(pval.espart,pval.espart.top)))
  dimnames(pval.espart)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.espart <- obs.espart/sum(obs.espart)
  psim.espart <- sweep(sim.espart,2,apply(sim.espart,1,sum),"/")
  bds.espart <- apply(psim.espart,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for edgewise shared partner\n\n")
# print(pval.espart)

 if ('dspartners' %in% all.gof.vars) {
  pval.dspart <- apply(sim.dspart <= obs.dspart[col(sim.dspart)],2,mean)
  pval.dspart.top <- apply(sim.dspart >= obs.dspart[col(sim.dspart)],2,mean)
  pval.dspart <- cbind(obs.dspart,apply(sim.dspart, 2,min), apply(sim.dspart, 2,mean),
                apply(sim.dspart, 2,max), pmin(1,2*pmin(pval.dspart,pval.dspart.top)))
  dimnames(pval.dspart)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.dspart <- obs.dspart/sum(obs.dspart)
  psim.dspart <- sweep(sim.dspart,2,apply(sim.dspart,1,sum),"/")
  bds.dspart <- apply(psim.dspart,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for dyadwise shared partner\n\n")
# print(pval.dspart)

 if ('triadcensus' %in% all.gof.vars) {
  pval.triadcensus <- apply(sim.triadcensus <= obs.triadcensus[col(sim.triadcensus)],2,mean)
  pval.triadcensus.top <- apply(sim.triadcensus >= obs.triadcensus[col(sim.triadcensus)],2,mean)
  pval.triadcensus <- cbind(obs.triadcensus,apply(sim.triadcensus, 2,min), apply(sim.triadcensus, 2,mean),
                apply(sim.triadcensus, 2,max), pmin(1,2*pmin(pval.triadcensus,pval.triadcensus.top)))
  dimnames(pval.triadcensus)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.triadcensus <- obs.triadcensus/sum(obs.triadcensus)
  psim.triadcensus <- sweep(sim.triadcensus,2,apply(sim.triadcensus,1,sum),"/")
  bds.triadcensus <- apply(psim.triadcensus,2,quantile,probs=c(0.025,0.975))
 }

# Return

  returnlist <- list(n,
   pval.model, pval.triadcensus, pval.dist, pval.ideg, pval.odeg, pval.deg, pval.espart, pval.dspart,
   obs.model, pobs.model, sim.model, psim.model, pval.model, bds.model,
   obs.triadcensus, pobs.triadcensus, sim.triadcensus, psim.triadcensus, pval.triadcensus, bds.triadcensus,
   obs.dist, pobs.dist, sim.dist, psim.dist, pval.dist, bds.dist,
   obs.ideg, pobs.ideg, sim.ideg, psim.ideg, pval.ideg, bds.ideg,
   obs.odeg, pobs.odeg, sim.odeg, psim.odeg, pval.odeg, bds.odeg,
   obs.deg, pobs.deg, sim.deg, psim.deg, pval.deg, bds.deg,
   obs.espart, pobs.espart, sim.espart, psim.espart, pval.espart, bds.espart,
   obs.dspart, pobs.dspart, sim.dspart, psim.dspart, pval.dspart, bds.dspart,
   GOF
                   )

  names(returnlist) <- c(
  "network.size",
  "summary.model",
  "summary.triadcensus",
  "summary.dist",
  "summary.ideg",
  "summary.odeg",
  "summary.deg",
  "summary.espart",
  "summary.dspart",
  "obs.model", "pobs.model", "sim.model", "psim.model", "pval.model", "bds.model",
  "obs.triadcensus", "pobs.triadcensus", "sim.triadcensus", "psim.triadcensus", "pval.triadcensus", "bds.triadcensus",
  "obs.dist", "pobs.dist", "sim.dist", "psim.dist", "pval.dist", "bds.dist",
  "obs.ideg", "pobs.ideg", "sim.ideg", "psim.ideg", "pval.ideg", "bds.ideg",
  "obs.odeg", "pobs.odeg", "sim.odeg", "psim.odeg", "pval.odeg", "bds.odeg",
  "obs.deg", "pobs.deg", "sim.deg", "psim.deg", "pval.deg", "bds.deg",
  "obs.espart", "pobs.espart", "sim.espart", "psim.espart", "pval.espart", "bds.espart",
  "obs.dspart", "pobs.dspart", "sim.dspart", "psim.dspart", "pval.dspart", "bds.dspart",
  "GOF"
                        )
  class(returnlist) <- "gofobject"
  returnlist
  }

print.gofobject <- function(x, ...){

 all.gof.vars <- all.vars(x$GOF)

# match variables

 for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'triadcensus', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree','model'
     )
                               )
 }
 GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

 for(statname in all.gof.vars){
   if ('model' == statname) {
    cat("\nGoodness-of-fit for model statistics\n\n")
    print(x$summary.model)
   }

   if ('distance' == statname) {
    cat("\nGoodness-of-fit for minimum geodesic distance\n\n")
    print(x$summary.dist)
   }
  
   if ('idegree' == statname) {
    cat("\nGoodness-of-fit for in degree\n\n")
    print(x$summary.ideg)
   }
  
   if ('odegree' == statname) {
    cat("\nGoodness-of-fit for out degree\n\n")
    print(x$summary.odeg)
   }
  
   if ('degree' == statname) {
    cat("\nGoodness-of-fit for degree\n\n")
    print(x$summary.deg)
   }
  
   if ('espartners' == statname) {
    cat("\nGoodness-of-fit for edgewise shared partner\n\n")
    print(x$summary.espart)
   }
  
   if ('triadcensus' == statname) {
    cat("\nGoodness-of-fit for triad census\n\n")
    print(x$summary.triadcensus)
   }
  
   if ('dspartners' == statname) {
    cat("\nGoodness-of-fit for dyadwise shared partner\n\n")
    print(x$summary.dspart)
   }
  }

  invisible()
}

summary.gofobject <- function(object, ...){

 all.gof.vars <- all.vars(object$GOF)

# match variables

 for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'triadcensus', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree','model'
     )
                               )
 }
 GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

 for(statname in all.gof.vars){
   if ('model' == statname) {
    cat("\nGoodness-of-fit for model statistics\n\n")
    print(object$summary.model)
   }

   if ('distance' == statname) {
    cat("\nGoodness-of-fit for minimum geodesic distance\n\n")
    print(object$summary.dist)
   }
  
   if ('idegree' == statname) {
    cat("\nGoodness-of-fit for in degree\n\n")
    print(object$summary.ideg)
   }
  
   if ('odegree' == statname) {
    cat("\nGoodness-of-fit for out degree\n\n")
    print(object$summary.odeg)
   }
  
   if ('degree' == statname) {
    cat("\nGoodness-of-fit for degree\n\n")
    print(object$summary.deg)
   }
  
   if ('espartners' == statname) {
    cat("\nGoodness-of-fit for edgewise shared partner\n\n")
    print(object$summary.espart)
   }
  
   if ('triadcensus' == statname) {
    cat("\nGoodness-of-fit for triad census\n\n")
    print(object$summary.triadcensus)
   }
  
   if ('dspartners' == statname) {
    cat("\nGoodness-of-fit for dyadwise shared partner\n\n")
    print(object$summary.dspart)
   }
  }

  invisible()
}

plot.gofobject <- function(x, ..., 
         cex.axis=0.7, plotodds=FALSE,
         main="Goodness-of-fit diagnostics", 
         normalize.reachability=FALSE,
         verbose=FALSE) {

 color <- "gray75"
#par(oma=c(0.5,2,1,0.5))

#statsno <- (sum(stats=='deg')>0) + (sum(stats=='espart')>0) + (sum(stats=='d
 all.gof.vars <- all.vars(x$GOF)
 statsno <- length(all.gof.vars)

# match variables

 for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'triadcensus', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree', 'model'
     )
                               )
 }
 GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

 if(statsno==0){
  stop("The gof object does not contain any statistics!\n")
 }
 n <- x$network.size

#attach(x)
  
 ###model####

 for(statname in all.gof.vars){
  if ('model' == statname) {

   nstats <- length(x$obs.model)
   if( min(x$pval.model[,"MC p-value"]) <0) {
    pval.max <- max((1:nstats)[x$pval.model[1:nstats, "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:nstats)[x$obs.model[1:nstats] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < nstats) {
        model <- c(1:pval.max)
    }
    else {
        model <- c(1:nstats)
    }
    if (plotodds) {
        odds <- x$psim.model
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.model
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.model
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for the statistic"
    }
    else {
        out <- x$psim.model
        out.obs <- x$pobs.model
        out.bds <- x$bds.model
        ylab <- "statistic"
    }
    pnames <- names(out.obs)
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, model]), xlab = "model statistics", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	ylim=c(ymin,ymax)
	   )

    points(seq(along = model), out.bds[1,model], pch = 1,cex=0.75)
    points(seq(along = model), out.bds[2,model], pch = 1,cex=0.75)
    lines(seq(along = model), out.bds[1, model], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = model), out.bds[2, model], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = model), out.obs[model], pch = 16,cex=0.75)
    lines(seq(along = model), out.obs[model], lty = 1,lwd=3)
  }

 ###degree####

  if ('degree' == statname) {

   if( min(x$pval.deg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.deg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.deg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        deg <- c(1:pval.max)
    }
    else {
        deg <- c(1:n)
    }
    if (plotodds) {
        odds <- x$psim.deg
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.deg
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.deg
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.deg
        out.obs <- x$pobs.deg
        out.bds <- x$bds.deg
        ylab <- "proportion of nodes"
    }
    pnames <- c(deg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, deg]), xlab = "degree", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = deg), out.bds[1,deg], pch = 1,cex=0.75)
    points(seq(along = deg), out.bds[2,deg], pch = 1,cex=0.75)
    lines(seq(along = deg), out.bds[1, deg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = deg), out.bds[2, deg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = deg), out.obs[deg], pch = 16,cex=0.75)
    lines(seq(along = deg), out.obs[deg], lty = 1,lwd=3)
  }

  ###odegree####

  if ('odegree' == statname) {

   if( min(x$pval.odeg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.odeg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.odeg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        odeg <- c(1:pval.max)
    }
    else {
        odeg <- c(1:n)
    }
    if (plotodds) {
        odds <- x$psim.odeg
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.odeg
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.odeg
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.odeg
        out.obs <- x$pobs.odeg
        out.bds <- x$bds.odeg
        ylab <- "proportion of nodes"
    }
    pnames <- c(odeg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, odeg]), xlab = "out degree", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = odeg), out.bds[1,odeg], pch = 1,cex=0.75)
    points(seq(along = odeg), out.bds[2,odeg], pch = 1,cex=0.75)
    lines(seq(along = odeg), out.bds[1, odeg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = odeg), out.bds[2, odeg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = odeg), out.obs[odeg], pch = 16,cex=0.75)
    lines(seq(along = odeg), out.obs[odeg], lty = 1,lwd=3)
  }

  ###idegree####

  if ('idegree' == statname) {

   if( min(x$pval.ideg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.ideg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.ideg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        ideg <- c(1:pval.max)
    }
    else {
        ideg <- c(1:n)
    }
    if (plotodds) {
        odds <- x$psim.ideg
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.ideg
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.ideg
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.ideg
        out.obs <- x$pobs.ideg
        out.bds <- x$bds.ideg
        ylab <- "proportion of nodes"
    }
    pnames <- c(ideg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, ideg]), xlab = "in degree", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = ideg), out.bds[1,ideg], pch = 1,cex=0.75)
    points(seq(along = ideg), out.bds[2,ideg], pch = 1,cex=0.75)
    lines(seq(along = ideg), out.bds[1, ideg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = ideg), out.bds[2, ideg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = ideg), out.obs[ideg], pch = 16,cex=0.75)
    lines(seq(along = ideg), out.obs[ideg], lty = 1,lwd=3)
  }

  ###espart####

  if ('espartners' == statname) {

   pval.max <- max((1:(n - 1))[x$pval.espart[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        espart <- c(1:pval.max)
    }
    else {
        espart <- c(1:(n-1))
    }
    if (plotodds) {
        odds <- x$psim.espart
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.espart
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.espart
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for an edge"
    }
    else {
        out <- x$psim.espart
        out.obs <- x$pobs.espart
        out.bds <- x$bds.espart
        ylab <- "proportion of edges"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    pnames <- c(espart)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, espart]), xlab = "edge-wise shared partners", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = espart), out.bds[1,espart], pch = 1,cex=0.75)
    points(seq(along = espart), out.bds[2,espart], pch = 1,cex=0.75)
    lines(seq(along = espart), out.bds[1, espart], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = espart), out.bds[2, espart], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = espart), out.obs[espart], pch = 16, cex=0.75)
    lines(seq(along = espart), out.obs[espart], lty = 1,lwd=3)

  }

  ###dspart####

  if ('dspartners' == statname) {
   pval.max <- max((1:(n - 1))[x$pval.dspart[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        dspart <- c(1:pval.max)
    }
    else {
        dspart <- c(1:n)
    }
    if (plotodds) {
        odds <- x$psim.dspart
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.dspart
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.dspart
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for an edge"
    }
    else {
        out <- x$psim.dspart
        out.obs <- x$pobs.dspart
        out.bds <- x$bds.dspart
        ylab <- "proportion of dyads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    pnames <- c(dspart)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, dspart]), xlab = "dyad-wise shared partners", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = dspart), out.bds[1,dspart], pch = 1,cex=0.75)
    points(seq(along = dspart), out.bds[2,dspart], pch = 1,cex=0.75)
    lines(seq(along = dspart), out.bds[1, dspart], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = dspart), out.bds[2, dspart], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = dspart), out.obs[dspart], pch = 16,cex=0.75)
    lines(seq(along = dspart), out.obs[dspart], lty = 1,lwd=3)
  }

  ###triadcensus####

  if ('triadcensus' == statname) {

    triadcensus <- c(1:16)
    if (plotodds) {
        odds <- x$psim.triadcensus
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.triadcensus
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.triadcensus
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a triad"
    }
    else {
        out <- x$psim.triadcensus
        out.obs <- x$pobs.triadcensus
        out.bds <- x$bds.triadcensus
        ylab <- "proportion of triads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    pnames <- c("003","012", "102", "021D", "021U", "021C", "111D",
                "111U", "030T",
                "030C", "201", "120D", "120U", "120C", "210", "300")
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, triadcensus]), xlab = "triad census", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
	)

    points(seq(along = triadcensus), out.bds[1,triadcensus], pch = 1,cex=0.75)
    points(seq(along = triadcensus), out.bds[2,triadcensus], pch = 1,cex=0.75)
    lines(seq(along = triadcensus), out.bds[1, triadcensus], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = triadcensus), out.bds[2, triadcensus], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = triadcensus), out.obs[triadcensus], pch = 16, cex=0.75)
    lines(seq(along = triadcensus), out.obs[triadcensus], lty = 1,lwd=3)

  }

  ###distance####

  if ('distance' == statname) {

    pval.max <- max((1:(n - 1))[x$pval.dist[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        dist <- c(1:pval.max, n)
    }
    else {
        dist <- c(1:n)
    }
    pnames <- paste(dist)
    pnames[length(dist)] <- "NR"
    if (plotodds) {
        odds <- x$psim.dist
        odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.dist
        odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.dist
        odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        oodds <- is.infinite(odds) | is.na(odds)
        oodds.obs <- is.infinite(odds.obs) | is.na(odds.obs)
        oodds.bds <- is.infinite(odds.bds) | is.na(odds.bds)
        mininf <- min(min(odds[!oodds]),min(odds.obs[!oodds.obs]),min(odds.bds[!oodds.bds]))
        maxinf <- max(max(odds[!oodds]),max(odds.obs[!oodds.obs]),max(odds.bds[!oodds.bds]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        odds.bds[1,][is.na(odds.bds[1,])] <- mininf
        odds.bds[2,][is.na(odds.bds[2,])] <- maxinf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a dyad"
    }
    else {
        out <- x$psim.dist
        out.obs <- x$pobs.dist
        out.bds <- x$bds.dist
        ylab <- "proportion of dyads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }

    if(normalize.reachability){
      mdist <- max(dist,na.rm=TRUE)
      totrange <- range(out.bds[1,][out.bds[1,] > out.bds[1,mdist]],
                        out.bds[2,][out.bds[2,] < out.bds[2,mdist]])
      out[,mdist] <- (out[,mdist]-out.bds[1,mdist]) * 
        diff(totrange) / diff(out.bds[,mdist]) + totrange[1]
      out.obs[mdist] <- (out.obs[mdist]- out.bds[1,mdist]) *
        diff(totrange) / diff(out.bds[,mdist]) + totrange[1]
      out.bds[,mdist] <- totrange
    }

    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))
    if(!plotodds){
     ymin <- max(0,ymin)
     ymax <- min(1,ymax)
    }

    boxplot(data.frame(out[, dist]), xlab = "minimum geodesic distance", 
        ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
	  ylim=c(ymin,ymax)
    )

    points(seq(along = dist), out.bds[1,dist], pch = 1,cex=0.75)
    points(seq(along = dist), out.bds[2,dist], pch = 1,cex=0.75)
    lines(seq(along = dist)[-length(dist)], out.bds[1, dist][-length(dist)], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = dist)[-length(dist)], out.bds[2, dist][-length(dist)], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = dist), out.obs[dist], pch = 16,cex=0.75)
    lines(seq(along = dist)[-length(dist)], out.obs[dist][-length(dist)],
		 lty = 1,lwd=3)
    }
   }

   mtext(main,side=3,outer=TRUE,cex=1.5,padj=2)
   invisible()
  }
