\name{latent}
\alias{latent}
\title{Latent position models for networks}
\description{
  \code{\link{latent}} is a term to the %\code{\link{statnet}}
  function \code{\link{ergmm}}
  which is used to fit a latent position model to a 
  given network, \eqn{g}.
  \code{\link{ergmm}} returns a Bayesian model fit based on a 
  Monte Carlo scheme.
  The default prior specifications are diffuse. An approximate MLE fit is
  also returned. 

  The \code{\link{ergmm}} specifies models via: \code{y ~ latent(<options>)}
  where \code{y} is a \code{network} object.
  For the list of possible \code{<options>}, see below.
  For the list of other model terms, see
  the manual pages for \code{terms.ergmm}.
}
\usage{
 latent(k=2, z.delta=0.1, z.prior.mu=0,
        z.prior.sd=10, b.delta=0.5, b.prior.mu=0,
        b.prior.sd=10)
}
\arguments{
\item{k}{Dimension of the latent space.}
\item{z.delta}{Standard deviation of deviance in the proposal for the latent positions. If a constant is passed it is used for each dimension.}
\item{z.prior.mu}{Prior mean for each dimension of the latent positions. If a
constant is passed it is used for each dimension.}
\item{z.prior.sd}{Prior standard deviation for each dimension of the latent positions. If a constant is passed it is used for each dimension.}
\item{b.delta}{Standard deviation of the deviance for covariate parameters.
If a constant is passed it is used for each dimension.}
\item{b.prior.mu}{Prior mean for the covariate parameters. If a constant is
passed it is used for each dimension.}
\item{b.prior.sd}{Prior standard deviation for the covariate parameters. If a
constant is passed it is used for each dimension.}
%\item{\dots}{Specific to the model term.}
}
\value{\code{\link{ergmm}} returns an object of class \code{\link{ergmm}} that
is a list.
}
\seealso{latentcluster, plot.ergmm, sna, network, terms.ergmm}
\references{
 Peter D. Hoff,  Adrian E. Raftery and Mark S. Handcock.
 \emph{Latent space approaches to social network analysis.} 
 Journal of the American Statistical Association, Dec 2002, Vol.97, 
 Iss. 460;  pg. 1090-1098.
}
\examples{
%\dontshow{require(statnet)}
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Get the group labels
samp.labs <- substr(get.vertex.attribute(samplike,"group"),1,1)
#
samp.fit <- ergmm(samplike ~ latent(k=2), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit
#
plot(samp.fit,label=samp.labs, vertex.col="group")
}
\keyword{graphs}
