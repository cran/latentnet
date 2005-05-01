InitErgm.latentcluster<-function(g, model, k, ngroups, z.delta=0.1,
                                 z.prior.mu=0, z.prior.sd=10, b.delta=0.5,
                                 b.prior.mu=0, b.prior.sd=10,
                                 Sigprior = qchisq(0.05,3),
                                 muSigprior = 2, dirprior=3,
                                 alphaprior=3,
                                 chisqprop = 6, thetaprop=0.1, ...)
{
#                                Sigprior = 0.33,
#                                muSigprior = 2, dirprior=1, alphaprior=1,
#                                chisqprop = 2, thetaprop=0.1, ...)
    if (nargs()<3)
    stop(paste("latentcluster() model term expected at least 1 argument, got ",
            nargs()-2, sep=""), call.=FALSE)
    optionnumber<-1+length(model$options)
    model$options[[optionnumber]] <-
     list(name="latentcluster", soname="latentnet", inputs=c(0, -1,10,k, ngroups,
                          z.delta, z.prior.mu, z.prior.sd,
                          b.delta, b.prior.mu, b.prior.sd,Sigprior,muSigprior,
                          dirprior,alphaprior,chisqprop,thetaprop))
    model$latent <- TRUE
    model$cluster <- TRUE
    model
}
