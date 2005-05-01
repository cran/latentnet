rergm.latent <- function(object,mkl=TRUE,n=1,...)
{
  out.list <- list()
  theta0 <- object$beta.mkl
  if(mkl)
    {
      for(i in 1:n)
        out.list[[i]] <- network(rergm.latent.sociomatrix(object,mkl=mkl))
    }
  else
    {
      if(n > object$samplesize)
        n <- object$samplesize
      for(i in 1:n)
        out.list[[i]] <- network(rergm.latent.sociomatrix(object,mkl=mkl,which=i))
    }

  if(n>1)
  {
    out.list <- list(formula = object$formula, networks = out.list, 
                     stats = NULL, coef=theta0)
    class(out.list) <- "network.series"
  } else {
    out.list <- out.list[[1]]
  }
  return(out.list)
}

rergm.latent.sociomatrix <- function(object, ..., mkl=TRUE,which)
{
  if(is.null(object$newnetwork)){
    newnetwork <- object$network
  }else{
    newnetwork <- object$newnetwork
  }
  nnodes <- network.size(newnetwork)
  if(missing(which)) mkl <- TRUE
  
  if(!mkl)
    {
      z.dist <- as.matrix(dist(object$Z[,,which]))
      beta <- object$Beta[which,]
    }
  else
    {
    z.dist <- as.matrix(dist(object$Z.mkl))
    beta <- object$beta.mkl
  }
  ##currently only renames the match'ed variables
  varnames <- sub("nodematch.","",object$glm.names)[-1]
  vars <- list()
  eta <- matrix(beta[1],nnodes,nnodes)
  if(length(varnames)>0)
    for(i in 1:length(varnames))
      {
        vars[[i]] <- unlist(get.vertex.attribute(newnetwork,varnames[i]))
        eta <- eta + beta[i+1] * outer(vars[[i]],vars[[i]],"==")
      }
  if(is.latent.cluster(object))
    eta <- eta - z.dist * mean(object$Beta[,ncol(object$Beta)])
  else
    eta <- eta - z.dist
  eta <- exp(eta)/(1+exp(eta))
  Yij <- matrix(runif(nnodes*nnodes),nnodes,nnodes)
  Yij <- 1*(Yij < eta)
}
