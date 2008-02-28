print.ergmm.model<-function(x,...){
  cat("Exponential Random Graph Mixed Model definition\n")
  cat("Formula: ")
  print(x$formula)
  if(!is.null(x$response) && mode(x$response)=="character")
    cat("Attribute of interest:",x$response,"\n")
  cat("Family:",x$family,"\n")

  cat("Terms:\n")
  if(x$intercept) cat("- intercept\n")
  if(x$d>0){
    cat("- latent space of",x$d,"dimensions")
    if(x$G>0){
      cat(" in",x$G,"clusters\n")
    }else cat("\n")
  }
}