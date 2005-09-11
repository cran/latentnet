\name{ergmm}
\alias{ergmm}
\title{Fit a Latent Space Random Graph Model}
\description{
  \code{\link{ergmm}} is used to fit latent space and latent space cluster random network models, 
    as described in Hoff, Raftery and Handcock (2002) and 
    Handcock, Raftery and tantrum (2005).
  \code{\link{ergmm}} can return either a Bayesian model fit
  or an approximate MLE based on a Monte Carlo scheme.
% \code{\link{ergmm}} can also be used to fit linear exponential random network models, in which
% the probability of a given network, \eqn{g}, on a set of nodes is 
% \eqn{\exp(\theta^{T}S(g))/c(\theta)}, where 
% \eqn{S(g)} is a vector of network statistics,
% \eqn{\theta} is a parameter vector of the same length and \eqn{c(\theta)} is the
% normalizing constant for the distribution.
}
\usage{
ergmm(formula, theta0=NULL, 
     burnin=1000, MCMCsamplesize=1000, interval=10, maxit=5,
     latent.control=list(maxit=40,penalty.sigma=c(10,0.5),MLEonly=FALSE),
     returnMCMCstats=TRUE, randseed=NULL, 
     verbose=FALSE, \dots)
}
\arguments{
 \item{formula}{An \R formula object, of the form
    \code{g ~ <term 1> + <term 2> ...},
  where \code{g} is a network object or a matrix that can be coerced to a
  network object, and \code{<term 1>}, \code{<term 2>}, etc, are each 
  terms chosen from the list given below.
  To create a network object in \R, use the \code{network} function,
  then add nodal attributes to it using \code{set.vertex.attribute} 
  if necessary.}
 \item{theta0}{The parameter value used to generate
  the MCMC sample. By default the MPLE is used (\code{startatMPLE=TRUE}).}
 \item{burnin}{The number of proposals before any MCMC sampling
  is done.  Currently, there is no support for any check of the Markov
  chain mixing, so \code{burnin} should be set to a fairly large number.}
 \item{MCMCsamplesize}{Size of the sample of network statistics,
  randomly drawn from a given distribution on the set of all networks,
  returned by the Metropolis-Hastings algorithm.}
 \item{interval}{The number of proposals between sampled statistics.
  The program prints a warning if too few proposals are being
  accepted in any interval before each sample.
% (if the number of proposals between sampled observations
% ever equals an integral multiple of 100(1+the number of proposals accepted)).
  }
 \item{maxit}{The number of times the parameter for the MCMC
  should be updated by maximizing the MCMC likelihood. At each
  step the parameter is changed to the values that maximizes
  the MCMC likelihood based on the current sample. For each
  step both the MCMCsamplesize and the interval are increased by 10\%.}
 \item{latent.control}{Control
  variables for the latent space algorithm. This are used only if
  a \code{\link{latent}} term is included in the model. \code{maxit}
  sets the maximum number of iterations to use in the
  Quasi-Newton-Raphson algorithm to maximize the MCMC likelihood.
  \code{penalty.sigma} is the penalty on the norm of the
  latent distances to use in the penalized log-likelihood. The multiplier is 
  \code{1/(penalty.sigma[1]^2} so that smaller values offer greater penalties.
  The second component is the multiplier on the component effects. Values
  less than 1 reduce the repulsion between components.
  This is used in the MLE of positions only and not the MCMC
  log-likelihood. It can be interpreted as a surrogate for
  a prior distribution and expresses the belief that the latent distances are
  not too large on the log-odds scale. The default value is
  \code{latent.control=list(maxit=40,penalty.sigma=c(10,0.5))}.
  }
 \item{returnMCMCstats}{If this is \code{TRUE} the matrix of change 
  statistics from the MCMC run is returned as component \code{$sample}.
  This matrix is actually an object of class \code{mcmc} and can be 
  used directly in the \code{CODA} package to assess MCMC convergence.}
%See \url{http://www.mrc-bsu.cam.ac.uk/bugs/classic/coda04/readme.shtml}.}
 \item{randseed}{Random number integer seed.
   The default is \code{sample(10000000, size=1)}.}
 \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program, including (currently) some goodness of fit 
  statistics.}
 \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\section{Model Terms}{
  The \code{\link{latentnet}} package itself allows only three type of terms:
  \code{\link{latent}}, \code{\link{latentcluster}} and \code{\link[latentcov]{latentcov}}.
  The umbrella \code{ergm} package allows the user to explore a large number
  of potential models for their network data in addition to these terms.  The
  terms currently supported by the program,
  and a brief description of each is given in the documentation
  \code{\link[terms.ergm]{terms.ergm}} for the \code{ergm} package.
  In the formula for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.
  The current options are:
 \item{\code{latent(k=2, \ldots)}}{\emph{Latent position model term.}
  where \code{k} is the dimension of the latent space. For information on the
  other arguments look for help on \code{\link{latent}}.
 }
 \item{\code{latentcluster(k=2, ngroups, \dots)}}{\emph{Latent position cluster model term.}
  where \code{k} is the dimension of the latent space
  and \code{ngroups} is the number of clusters in the latent space.
  For information on the other arguments look for help on
  \code{\link{latentcluster}}.
 }
 \item{\code{latentcov(cv, attrname=NULL)}}{\emph{Covariates for the latent model.}  \code{cv} is either a matrix of
    covariates on each pair of vertices, or a network; 
    if the latter, optional argument
    \code{attrname} provides the name of the edge attribute to
    use for edge values.  This option adds one statistic to the
    model, representing the effect of the given covariate on
    the appearance of edges.  \code{edgecov} can be called more
    than once, to model the effects of multiple covariates.
 }
}
\value{\code{\link{ergmm}} returns an object of class \code{\link{ergmm}} that is a list
  consisting of the following elements:
 \item{\$coef}{The Monte Carlo maximum likelihood estimate
  of \eqn{\theta}, the vector of coefficients for the model
  parameters.}
 \item{\$sample}{The \eqn{n\times (p+2+k)} matrix of network statistics, 
  where \eqn{n} is the
  sample size and \eqn{p} is the number of network covariates specified in the
  model via the \code{latentcov} terms (usually 0). The columns are:
  "mcmc.loglikelihood", the log-likelihood value;
  "density", the constant term in the latent model;
  the \code{p} covariates;
  "Z 1", "Z 2", ..., "Z k", the \code{k} dimensional
   positions of the first node. The values are recorded for each sample drawn.}
 \item{\$iterations}{The number of Newton-Raphson iterations required
  before convergence.}
 \item{\$MCMCtheta}{The value of \eqn{\theta} used to produce the Markov chain
  Monte Carlo sample.  As long as the Markov chain mixes sufficiently
  well, \code{$sample} is roughly a random sample from the distribution
  of network statistics specified by the model with the parameter equal
  to \code{$MCMCtheta}.  In the current version, if startatMPLE is 
  \code{TRUE}, then \code{$MCMCtheta} equals the MPLE.}
 \item{\$loglikelihood}{The approximate log-likelihood for the MLE.
  The value is only approximate because it is based on the MCMC random sample.}
 \item{\$gradient}{The value of the gradient vector of the approximated
  loglikelihood function, evaluated at the maximizer.  This vector
  should be very close to zero.}
 \item{\$hessian}{The Hessian matrix of the approximated loglikelihood
  function, evaluated at the maximizer.  This matrix may be inverted to
  give an approximate covariance matrix for the MLE.}
 \item{\$samplesize}{The size of the MCMC sample}
 \item{\$formula}{The original formula entered into the \code{\link{ergmm}} function.}
 \item{\$sample}{If the option \code{$returnMCMCstats=TRUE}, this is 
  the the matrix of change statistics from the MCMC run.}
 \item{\$newnetwork}{The network generated at the end of the MCMC sampling.}

  See the function \code{print.ergmm} for details on how
  an \code{\link{ergmm}} object is printed.  Note that we have written a
  function, \code{\link{summary.ergmm}} that returns a summary of the
  relevant parts of the \code{\link{ergmm}} object in concise summary
  format.
}
\seealso{network, set.vertex.attributes, set.network.attributes,
         summary.ergmm, print.ergmm}
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
#
# See http://www.csde.washington.edu/statnet
# for examples
#
# Use 'data(package = "latentnet")' to list the data sets in a
#
data(package="latentnet")
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Get the group labels
#
group <- get.vertex.attribute(samplike,"group")
samp.labs <- substr(group,1,1)
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
\dontrun{
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
#
# Add contours
#
ergmm.add.contours(samp.fit,nlevels=8,lwd=2)
points(samp.fit$Z.mkl,pch=19,col=samp.fit$class)
#
# Try a covariate on the group
#
samegroup <- outer(group, group, "==")
diag(samegroup) <- 0
samp.fit <- ergmm(samplike ~ latentcov(samegroup) + latent(k=2))
summary(samp.fit)
}
}
\keyword{graphs}
