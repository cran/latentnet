\name{rergm}
\alias{rergm}
\alias{rergm.default}
\alias{rergm.ergm}
\title{Draw from the distribution of an Exponential Family Random Graph Model}
\usage{
rergm(object, \dots)

\method{rergm}{default}(object,\dots,prob,theta0,n=1,
                          directed=TRUE,numedges=NULL)
\method{rergm}{ergm}(object, \dots, theta0=NULL, n=1,
                     burnin=1000, interval=1000, 
                     randseed=NULL,
                     sequential=TRUE, summarizestats=FALSE,
                     verbose=FALSE)
}
\description{
  \code{rergm} is used to draw from exponential family random network
  models in their natural parameterizations.
  See \code{\link{ergm}} for more information on these models.
}
\arguments{
 \item{object}{an R object. Either a number of nodes in the network,
   a formula or an ergm object. See documentation for \code{\link{ergm}}.
   If the number of nodes in the network is given then
   Bernoulli networks are drawn.}
 \item{prob}{The probability of a link for Bernoulli networks.  Defaults
    to 0.5 if neither prob nor theta0 are given.}
 \item{theta0}{For Bernoulli networks this is the log-odds of a tie,
    however it is only used if prob is not specified.  When given either
    a formula or an object of class ergm, theta0 are the parameters
    from which the sample is drawn.}
 \item{directed}{Whether the Bernoulli network should be directed or undirected.}
 \item{numedges}{If present, sample the network(s) conditional on this 
    number of edges (rather than independently with the 
    specified probability).}
 \item{n}{Size of the sample of networks to be randomly drawn
    from the given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.}
 \item{burnin}{The number of proposed proposals before any MCMC sampling
    is done.  Currently, there is no support for any check of the Markov
    chain mixing, so \code{burnin} should be set to a fairly large number.}
 \item{interval}{The number of proposals between sampled statistics.
    The program prints a warning if too few proposals are being
    accepted (if the number of proposals between sampled observations
    ever equals an integral multiple of 100(1+the number of proposals
    accepted)).}
 \item{randseed}{Random number integer seed.
   The default is \code{sample(10000000, size=1)}.}
 \item{sequential}{Should the returned draws use the prior draw 
   as the starting network or always use the initially passed network?
   For random draws the results should be similar (stochastically), but the
   \code{sequential=TRUE} option is useful for dynamic draws.
   }
 \item{summarizestats}{Print out a summary of the sufficient statistics of
  the generated network. This is useful as a diagnostic.}
 \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program, including (currently) some goodness of fit 
  statistics.}
  \item{\dots}{further arguments passed to or used by methods.}
}

\details{
  A sample of networks is randomly drawn from the specified model.  The
  model is either specified by the first argument of the function.  If
  the first argument is a formula then this defines the model.
  If this is not given as the
  first argument then a Bernoulli network is generated with the probability
  of ties defined by \code{prob} or \code{theta0}.

  Note that the first network is sampled after \code{burnin} + \code{interval}
  steps, and any subsequent networks are sampled each
  \code{interval} steps after the first.

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}
\value{\code{rergm} returns an object of class \code{network.series} that is a list
  consisting of the following elements:
  \item{\$formula}{The formula used to generate the sample.}
  \item{\$networks}{A list of the generated networks.}
  \item{\$stats}{The \eqn{n\times p}  matrix of network change statistics, 
   where \eqn{n} is the sample size and \eqn{p} is the number
   of network change statistics specified in the model.}

}
\seealso{ergm, network, print.network}
\examples{
#
# Let's draw from a Bernoulli model with 16 nodes
# and tie probability 0.1
#
g.use <- rergm(16,prob=0.1,directed=FALSE)
#
data(sampson)
gest <- ergm(samplike ~ latent(k=2))
summary(gest)
#
# Draw from the fitted model
#
g.sim <- rergm(gest,n=100,burnin=1000,interval=1000)
g.sim
}
\keyword{graphs}
