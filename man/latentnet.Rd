\name{latentnet}
\alias{latentnet}
\title{Latent position and cluster models for networks}
\description{
  By way of context, the package \code{statnet} is used to fit linear exponential random network
  models, in which the probability of a given network, \eqn{g},
  on a set of nodes is \eqn{\exp(\theta^{T}S(g))/c(\theta)},
  where \eqn{S(g)} is a vector of network statistics, \eqn{\theta}
  is a parameter vector of the same length and \eqn{c(\theta)}
  is the normalizing constant for the distribution.  The \code{ergm}
  function can return either a maximum pseudo-likelihood estimate or
  an approximate MLE based on a Monte Carlo scheme.
  For details, see \url{http://www.csde.washington.edu/statnet/latentnet}.

  The package \code{\link{latentnet}} adds latent variable models 
  to the traditional 
  exponential random network models. These latent models are
  fit using the \code{\link{ergmm}} function.

  The \code{\link{ergmm}} specifies models via: \code{g ~ <model terms>}
  where \code{g} is a \code{network} object.
  For the list of possible \code{<model terms>}, see 
  the manual pages for \code{terms.ergmm}.
}
\usage{
}
\arguments{
\item{\dots}{Specific to the latent variable model terms.}
}
%\usage{
%ergmm(formula, theta0=NULL, 
%     burnin=1000, MCMCsamplesize=1000, interval=10,
%     latent.control=list(maxit=40,penalty.sigma=c(10,0.5),MLEonly=FALSE),
%     returnMCMCstats=TRUE, randseed=NULL, 
%     verbose=FALSE, \dots)
%}
%\arguments{
% \item{formula}{An \R formula object, of the form 
%    \code{g ~ <term 1> + <term 2> ...},
%  where \code{g} is a network object or a matrix that can be coerced to a graph
%  object, and \code{<term 1>}, \code{<term 2>}, etc, are each model terms.
%  For the list of possible \code{<model terms>}, see
%  the manual pages for \code{terms.ergmm}.
%  To create a network object in \R, use the \code{network()} function,
%  then, if necessary, add nodal
%  attributes to it using \code{set.vertex.attributes()}.}
%}
%
\value{\code{\link{ergmm}} returns an object of class \code{\link{ergmm}} that
is a list.
}
\seealso{ergmm, latent, latentcluster, sna, network, terms.ergmm}
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
\keyword{graphs}
