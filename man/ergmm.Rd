\name{ergmm}
\alias{ergmm}
\title{Fit a Latent Space Random Graph Model}
\description{
  \code{\link{ergmm}} is used to fit latent space and latent space cluster random network models, 
    as described in Hoff, Raftery and Handcock (2002) and 
    Handcock, Raftery and Tantrum (2005).
  \code{\link{ergmm}} produces likelihood-based inference. Approximate maximum likelihood estimators are computed, and Bayesian inference is implemented via a MCMC algorithm.
% \code{\link{ergmm}} can also be used to fit linear exponential random network models, in which
% the probability of a given network, \eqn{g}, on a set of nodes is 
% \eqn{\exp(\theta{\cdot}S(g))/c(\theta)}, where 
% \eqn{S(g)} is a vector of network statistics,
% \eqn{\theta} is a parameter vector of the same length and \eqn{c(\theta)} is the
% normalizing constant for the distribution.
}
\usage{
ergmm(formula, theta0=NULL, 
     burnin=1000, MCMCsamplesize=1000, interval=10,
     latent.control=list(maxit=40,penalty.sigma=c(10,0.5),MLEonly=FALSE),
     returnMCMCstats=TRUE, randseed=NULL, 
     verbose=FALSE, \dots)
}
\arguments{
 \item{formula}{An \R formula object, of the form
    \code{y ~ <term 1> + <term 2> ...},
  where \code{y} is a network object or a matrix that can be coerced to a
  network object, and \code{<term 1>}, \code{<term 2>}, etc, are each 
  terms chosen from the list given below.
  To create a network object in \R, use the \code{network} function.
  For a description of the possible terms see the \code{\link{terms.ergmm}}.
 }
 \item{theta0}{The initial parameter value used to find
  the MLE. The default is based on
  multidimensional scaling fit to the positions.}
 \item{burnin}{The number of proposals before any MCMC sampling
  is done.}
 \item{MCMCsamplesize}{The number of posterior samples to draw.}
 \item{interval}{The number of proposal steps between sampled statistics.}
% The program prints a warning if too few proposals are being
% accepted in any interval before each sample.
% (if the number of proposals between sampled observations
% ever equals an integral multiple of 100(1+the number of proposals accepted)).
% }
%\item{maxit}{The number of times the parameter for the MCMC
% should be updated by maximizing the MCMC likelihood. At each
% step the parameter is changed to the values that maximizes
% the MCMC likelihood based on the current sample. For each
% step both the MCMCsamplesize and the interval are increased by 10\%.}
 \item{latent.control}{Control
  variables for the latent space algorithm. This are used only if
  a \code{\link{latent}} term is included in the model. \code{maxit}
  sets the maximum number of iterations to use in the
  Quasi-Newton-Raphson algorithm to maximize the MCMC likelihood.
  \code{MLEonly} is a logical flag set to compute only the MLE estimates and not
  Bayesian inference based on the MCMC algorithm.
  \code{penalty.sigma} is the penalty on the norm of the
  latent distances to use in the penalized log-likelihood. The multiplier is 
  \code{1/(penalty.sigma[1]^2} so that smaller values offer greater penalties.
  The second component is the multiplier on the component effects. Values
  less than 1 reduce the repulsion between components.
  This is used in the MLE of positions only and not the MCMC
  log-likelihood. It can be interpreted as a surrogate for
  a prior distribution and expresses the belief that the latent distances are
  not too large on the log-odds scale. 
% The default value is
% \code{latent.control=list(maxit=40,penalty.sigma=c(10,0.5))}.
  }
 \item{returnMCMCstats}{If this is \code{TRUE} the matrix of change 
  statistics from the MCMC run is returned as component \code{sample}.
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

\value{\code{\link{ergmm}} returns an object of class \code{\link{ergmm}} that is a list. 
  Fits including a \code{latentcluster} term  will have at least the following
components and fits including a \code{latent} term  will have at least
the components up to and including \code{network}.

 \item{coef}{The maximum likelihood estimate
  the \eqn{p} vector of coefficients for the model
  parameters (excluding the latent positions and cluster parameters). 
  By default this is just the intercept with \eqn{p=1}.}
 \item{coef.names}{A \eqn{p} vector of the coefficient names.}
 \item{Beta}{The \code{MCMCsamplesize}\eqn{\times p} matrix of coefficients for the model
  parameters corresponding to each of the posterior samples. By default this is
  the intercept only.}
 \item{Z}{The \code{MCMCsamplesize}\eqn{\times k} matrix of (Procrustified)
  posterior positions, where \code{MCMCsamplesize} is the
  sample size and \eqn{k} is the number of dimensions of the latent space.}
 \item{Z.mkl}{The \code{network.size(g)}\eqn{\times} \code{k} matrix of 
  minimum Kullback-Leibler positions for each of the nodes.}
 \item{Z.pmean}{The \code{network.size(g)}\eqn{\times} \code{k} matrix of 
  posterior mean positions for each of the nodes.}
 \item{Z.pmode}{The \code{network.size(g)}\eqn{\times} \code{k} matrix of 
  posterior modal positions for each of the nodes.}
 \item{Z.mle}{The \code{network.size(g)}\eqn{\times} \code{k} matrix of 
  MLE positions for each of the nodes.}
 \item{beta.mkl}{The \eqn{p} vector of coefficients for the model
  parameters based on the minimum Kullback-Leibler positions for each of the nodes.}
 \item{samplesize}{The number of MCMC samples drawn from the posterior.}
 \item{sample}{The \code{MCMCsamplesize}\eqn{\times (p+2+k)} matrix of network statistics, 
  where \code{MCMCsamplesize} is the
  sample size and \eqn{p} is the number of network covariates specified in the
  model via the \code{latentcov} terms (usually 0). The columns are:
  ``mcmc.loglikelihood", the log-likelihood value;
  ``density", the constant term in the latent model;
  the \code{p} covariates;
  ``Z 1", ``Z 2", ..., ``Z k", the \code{k} dimensional
   positions of the first node. The values are recorded for each sample drawn.
   This is primarily used for MCMC diagnostics to assess convergence.}
 \item{iterations}{The number of Newton-Raphson iterations required
  before convergence.}
 \item{interval}{The number of proposals between sampled statistics.}
%\item{MCMCtheta}{The value of used to produce the Markov chain
% Monte Carlo sample.  As long as the Markov chain mixes sufficiently
% well, \code{sample} is roughly a random sample from the distribution
% of network statistics specified by the model with the parameter equal
% to \code{MCMCtheta}.  In the current version, if startatMPLE is 
% \code{TRUE}, then \code{MCMCtheta} equals the MPLE.}
 \item{null.deviance}{The deviance for the null model, comparable with
    \code{-2 loglikelihood}. The null model will include the 
    intercept if there is one in the model, but not the latent variables or latent
    clusters.}
 \item{mcmc.loglikelihood}{The log-likelihood values corresponding to each of the posterior samples.}
 \item{loglikelihood}{The log-likelihood for the MLE of positions (and based on
the final fits to the other parameters).}
 \item{mle.lik}{The log-likelihood for the initial MLE fit of positions.}
% The value is only approximate because it is based on the MCMC random sample.}
%\item{gradient}{The value of the gradient vector of the approximated
% loglikelihood function, evaluated at the maximizer.  This vector
% should be very close to zero.}
 \item{hessian}{The Hessian matrix of the approximated loglikelihood
  function, evaluated at the maximizer.  This matrix may be inverted to
  give an approximate covariance matrix for the MLE of the parameters.}
 \item{formula}{The original formula entered into the \code{\link{ergmm}} function.}
 \item{latent}{A flag to indicate that this is a fit of latent variable model.
  This is always \code{TRUE} for \code{ergmm} fits 
  and is included for consistency with the \code{statnet}
  package.} 
 \item{cluster}{A flag to indicate that this is a fit of a latent cluster model.
  This is always \code{TRUE} for \code{ergmm} fits if a \code{latentcluster}
  term is in the model
  and is included for consistency with the \code{statnet}
  package.} 
%\item{newnetwork}{The network generated at the end of the MCMC sampling.}
 \item{network}{The modeled network as an \code{network} object.}
 \item{BIC}{A Bayesian Information Criterion approximation for the model. 
            This is the approximation based on the fully 
            Bayesian estimation method in Section 3.2 of
            Handcock, Raftery and Tantrum (2005).  
            The formula for the approximation is given at the end of Section 4 in that paper.
            See the references for details.}
 \item{class}{The vector of posterior modal classes for each node.}
%\item{d.mbc}{}
 \item{Ki}{The \code{MCMCsamplesize}\eqn{\times}\code{network.size(g)}
  matrix of posterior draws of the classes, where \code{MCMCsamplesize} is the
  sample size and \code{network.size(g)} is the number of nodes in the
  network.}
 \item{Ki.mle}{The \code{network.size(g)} vector of maximum likelihood classes for each node.}
%\item{logl}{The log-likelihood.}
 \item{logl.lr}{The log-likelihood for the latent space component of the model.}
 \item{logl.mbc}{The log-likelihood for the model-based clustering component of
the model.}
 \item{mu}{The \code{ngroups}\eqn{\times}\code{k}\eqn{\times}\code{MCMCsamplesize} array of posterior draws of the mean positions
  of the class, where \code{MCMCsamplesize} is the
  sample size and \code{ngroups} is the number of classes.}
 \item{mu.mle}{The \code{ngroups}\eqn{\times}\code{k} matrix of 
  maximum likelihood mean positions for each class.}
 \item{ngroups}{The number of classes or clusters.}
 \item{qig}{The \code{network.size(g)}\eqn{\times}\code{ngroups}
  matrix of posterior probabilities of class membership for each of the nodes.}
 \item{Sigma}{The \code{MCMCsamplesize}\eqn{\times}\code{ngroups}
  array of posterior draws of the variances of the positions
  of the class, where \code{MCMCsamplesize} is the
  sample size and \code{ngroups} is the number of classes.}
 \item{Sigma.mle}{The maximum likelihood variances of the
  positions for each class.}
% See the function \code{print.ergmm} for details on how
% an \code{\link{ergmm}} object is printed.
  Note that we have written a
  function, \code{\link{summary.ergmm}} that returns a summary of the
  relevant parts of the \code{\link{ergmm}} object in concise summary
  format.
}
\seealso{network, set.vertex.attributes, set.network.attributes,
         summary.ergmm}%, print.ergmm}
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
# See http://www.csde.washington.edu/statnet/latentnet
# for more examples
#
# For an explanation and examples of creating 'network' objects
# see the required 'network' package.
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
