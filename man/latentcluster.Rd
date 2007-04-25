\name{latentcluster}
\alias{latentcluster}
\title{Latent cluster models for networks}
\description{
  \code{latentcluster()} is a term to the  %\code{\link{statnet}}
  function \code{\link{ergmm}}
  to fit a latent position cluster model to a 
  given network, \eqn{g}.
  \code{\link{ergmm}} returns a Bayesian model fit based on a Monte Carlo scheme.
  The default prior specifications are diffuse. An approximate MLE fit is
  also returned. 

  The \code{\link{ergmm}} specifies models via: \code{y ~ latentcluster(<options>)}
  where \code{y} is a \code{network} object
  For the list of possible \code{<options>}, see below.
  For the list of other model terms, see
  the manual pages for \code{terms.ergmm}.
}
\usage{
 latentcluster(k=2, ngroups, z.prior.mu=0, z.prior.sd=10, b.delta=0.5,
               b.prior.mu=0, b.prior.sd=10,
               Sigprior = qchisq(0.05,3),
               muSigprior = 2, dirprior=3,
               alphaprior=3,
               chisqprop = 6, thetaprop=0.1, \dots)
}
\arguments{
\item{k}{Dimension of the latent space.}
\item{ngroups}{Number of clusters in the latent space.}
\item{z.prior.mu}{Prior mean for each dimension of the latent positions. If a
constant is passed it is used for each dimension.}
\item{z.prior.sd}{Prior standard deviation for each dimension of the latent positions. If a constant is passed it is used for each dimension.}
\item{b.delta}{Standard deviation of the deviance for covariate parameters.
If a constant is passed it is used for each dimension.}
\item{b.prior.mu}{Prior mean for the covariate parameters. If a constant is
passed it is used for each dimension.}
\item{b.prior.sd}{Prior standard deviation for the covariate parameters. If a
constant is passed it is used for each dimension.}
\item{Sigprior}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{muSigprior}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{dirprior}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{alphaprior}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{chisqprop}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{thetaprop}{Prior standard deviations for the node positions relative to
the cluster mean. If a constant is passed it is used for each dimension.}
\item{\dots}{Specific to the model term.}
}
\value{\code{\link{ergmm}} returns an object of class \code{\link{ergmm}} that
is a list.
}
\seealso{latent, plot.ergmm, sna, network, terms.ergmm}
\references{
 Peter D. Hoff,  Adrian E. Raftery and Mark S. Handcock.
 \emph{Latent space approaches to social network analysis.} 
 Journal of the American Statistical Association, Dec 2002, Vol.97, 
 Iss. 460;  pg. 1090-1098.

 Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum.
 \emph{Model-Based Clustering for Social Networks.} 
 Working Paper Number 46, Center for Statistics and the Social Sciences,
 University of Washington, April 2005.
}
\examples{
%\dontshow{require(statnet)}
\dontrun{
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Get the group labels
samp.labs <- substr(get.vertex.attribute(samplike,"group"),1,1)
#
# Using Sampson's Monk data, lets fit a latent clustering model
#
samp.fit <- ergmm(samplike ~ latentcluster(k=2, ngroups=3), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Lets look at the goodness of fit:
#
plot(samp.fit,label=samp.labs, vertex.col="group")
plot(samp.fit,pie=TRUE,label=samp.labs)
plot(samp.fit,density=c(2,2))
plot(samp.fit,contours=5,contour.color="red")
plot(samp.fit,density=TRUE,drawarrows=TRUE)
ergmm.add.contours(samp.fit,nlevels=8,lwd=2)
points(samp.fit$Z.mkl,pch=19,col=samp.fit$class)
}
}
\keyword{graphs}
