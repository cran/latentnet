rergm <- function(object, ...)
UseMethod("rergm")

rergm.default <- function(object,...,prob,theta0,n=1,directed=TRUE,
                          numedges=NULL)
{
  #should take in nnodes, and one of prob or theta0
  #returns a bernouli network.
  if(is.latent(object))
    return(rergm.latent(object,n=n,...))
  nnodes <- object
  if(directed){
   ndyads <- nnodes*(nnodes-1)
  }else{
   ndyads <- nnodes*(nnodes-1)/2
  }

  if(missing(prob)){
    if(missing(theta0)){
#     prob <- 0.5
#
#     So the expected number of ties is the same as
#     the number of nodes
#
      prob <- nnodes/ndyads
    }else{
      prob <- exp(theta0)/(1+exp(theta0))
    }
  }
# cat(paste("prob =",prob,"\n"))

  networks <- list()
  for(k in 1:n)
  {
    networks <- vector("list",length=n)
    g.mat <- matrix(0,nnodes,nnodes)
    dimnames(g.mat) <- list(1:nnodes,1:nnodes)
    if(is.null(numedges)){
     gij <- runif(ndyads)<prob
    }else{
     gij <- rep(0,ndyads)
     gij[sample(1:ndyads,size=numedges,replace=FALSE)] <- 1
    }
    if(directed){
     g.mat[row(g.mat)!=col(g.mat)] <- gij
    }else{
     g.mat[row(g.mat) < col(g.mat)] <- gij
     g.mat <- g.mat + t(g.mat)
#    g.mat[col(g.mat) < row(g.mat)] <- gij
    }
    networks[[k]] <- network(g.mat,directed=directed)
  }
  if(n > 1){
   out.list <- list(formula = ~1, networks = networks,
                    stats = numeric(0),coef=prob)
   class(out.list) <- "network.series"
  }else{
   out.list <- networks[[1]]
  }
  return(out.list)
}
rergm.ergm <- function(object, ..., theta0=NULL, n=1,
                       burnin=1000, interval=1000, 
                       randseed=NULL,
		       sequential=TRUE, summarizestats=FALSE,
		       verbose=FALSE)
{
  if(is.latent(object)){
    return(rergm.latent(object,n=n,...))
  }
}    
