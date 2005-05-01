InitErgm.latentcov<-function (g, model, x, attrname=NULL, drop=TRUE, ...) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 4:5))
    stop(paste("latentcov() model term expected 1 or 2 arguments, got ", 
                                   nargs() - 3, sep = ""), call. = FALSE)
  if (nargs()==4){drop <- attrname; attrname <- NULL}
  #Coerce x to an adjacency matrix
  if(is.network(x))
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
  else if(is.character(x))
    xm<-as.matrix.network(g,matrix.type="adjacency",x)
  else
    xm<-as.matrix(x)
  #Update the option number
  optionnumber <- 1 + length(model$options)
  #Update the options list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_latentcov ignores the value of inp->attrib).
  model$options[[optionnumber]] <- list(name = "latentcov", soname="latentnet", 
            inputs = c(1, 1, 1+NROW(xm)*NROW(xm), NROW(xm), as.double(xm)))
  #Update the coefficient name list, adding dyadcov.nameofx
  if(nargs()==5)
    cn<-paste("latentcov", as.character(sys.call(0)[[4]]), 
                         as.character(sys.call(0)[[5]]), sep = ".")
  else
    cn<-paste("latentcov", as.character(sys.call(0)[[4]]), sep = ".")
  model$coef.names <- c(model$coef.names, cn)
  #Return the updated model list
  model
}
