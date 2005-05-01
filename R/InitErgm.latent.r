InitErgm.latent<-function(g, model, k, z.delta=0.1, z.prior.mu=0,
                          z.prior.sd=10, b.delta=0.5, b.prior.mu=0,
                          b.prior.sd=10,...)
{
    if (nargs()<3)
        stop(paste("latent() model term expected at least 1 argument, got ", 
            nargs()-2, sep=""), call.=FALSE)
    optionnumber<-1+length(model$options)
    model$options[[optionnumber]] <-
      list(name="latent", soname="latentnet", inputs=c(0, -1,10,k,
                           z.delta, z.prior.mu, z.prior.sd,
                           b.delta, b.prior.mu, b.prior.sd))
    model$latent <- TRUE
    model
}
