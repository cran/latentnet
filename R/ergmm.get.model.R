ergmm.get.model <- function(formula,response,family,fam.par,prior){
  
  terms<-terms(formula)

  if(!attr(terms,"response") || terms[[1]]!="~") stop("Formula must be of form 'network ~ model'.")

  Yg <- try(as.network(eval(terms[[2]],attr(terms,".Environment"))))
  if(inherits(Yg,"try-error")){
    stop("Invalid network. Is the left-hand-side of the formula correct?")
  }

  model<-list(formula=formula,
              Yg=Yg,
              Ym=getYm(Yg,response),
              response=response,
              family=family,
              familyID=family.IDs[[family]],
              fam.par=fam.par,
              coef.names=character(0),
              X=list(),
              p=0,
              d=0,
              G=0,
              intercept=as.logical(attr(terms,"intercept")),
              prior=list() ## Only here for convenience.
              )

  model<-fam.par.check(model)
  
  if(model$intercept){
    model<-InitErgmm.latentcov(model,matrix(1,network.size(Yg),network.size(Yg)),"density")
  }
              
  for (term in as.list(attr(terms,"variables"))[-(1:2)]){
    if (is.call(term)){
      init.call<-list()
      init.call<-list(as.name(paste("InitErgmm.", term[[1]], sep = "")),
                      model=model)
      
      init.call<-c(init.call,as.list(term)[-1])
    }else{
      init.call <- list(as.name(paste("InitErgmm.", term, sep = "")),model=model)
    }
    model <- eval(as.call(init.call), attr(terms,".Environment"))
  }
  
  if(!("Z.var" %in% names(model$prior))) model$prior$Z.var<-model$prior$Z.var.mul*(network.size(model$Yg)/max(1,model$G))^(2/model$d)
  if(!("Z.mean.var" %in% names(model$prior))) model$prior$Z.mean.var<-model$prior$Z.mean.var.mul*model$prior$Z.var*max(1,model$G)^(2/model$d)
  if(!("Z.var.df" %in% names(model$prior))) model$prior$Z.var.df<-model$prior$Z.var.df.mul*sqrt(network.size(model$Yg)/max(1,model$G))
  if(!("Z.pK" %in% names(model$prior))) model$prior$Z.pK<-model$prior$Z.pK.mul*sqrt(network.size(model$Yg)/max(1,model$G))
  
  if(prior$adjust.beta.var) model$prior$beta.var<-model$prior$beta.var/sapply(1:model$p,function(i) mean((model$X[[i]][observed.dyads(model$Yg)])^2))
  
  for(name in names(prior)){
    model$prior[[name]]<-prior[[name]]
  }

  prior<-model$prior
  model$prior<-NULL
  
  class(model)<-"ergmm.model"  
  list(model=model,prior=prior)
}