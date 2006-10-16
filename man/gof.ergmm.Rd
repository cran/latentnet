\name{gof.ergmm}
\alias{gof.ergmm}
\alias{gof.default}
\alias{gof}
\title{Conduct Goodness-of-Fit Diagnostics on a Latent Position or Latent Cluster Position Model}
\description{
  \code{\link{gof}} computes and plots posterior predictive checks for various
  networks statistics (e.g., geodesic
  distance, degree, and reachability) to
  diagnose the goodness-of-fit of  Latent Position or Latent Cluster Position
  models.  See \code{\link{ergmm}} for more information on these models.
}
\usage{
\method{gof}{default}(object,\dots)
\method{gof}{ergmm}(object, \dots,
         nsim=100,
         GOF=~degree+espartners+distance,
         verbose=FALSE)
}
\arguments{
 \item{object}{an ergmm object. 
   See documentation for \code{\link{ergmm}}.}
  \item{nsim}{The number of posterior draws to use for the posterior predictive
   checks.
   This is the size of the sample of networks to be randomly drawn
   from the posterior distribution specified by the 
   object.
   }
  \item{GOF}{formula; an \R formula object, of the form
   \code{~ <model terms>} specifying the
   statistics to use to diagnosis the goodness-of-fit of the model. 
   They do not need to be in the model formula specified in
   \code{formula}, and typically are not.
   Examples are the degree distribution ("degree"), 
   minimum geodesic distances ("distance"), and shared partner distributions
   ("espartners" and "dspartners").
   For the details on the possible
   \code{<model terms>}, see \code{\link[statnet]{ergm-terms}}.}
  \item{verbose}{Provide verbose information on the progress of the
    simulation.}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\details{
  A sample of graphs is randomly drawn from the posterior
  of the specified model.
  The first argument is 
  the output of a call to \code{\link{ergmm}} and the model
  used for that call is the one evalauted.

  A plot of the summary measures is plotted.
}
\value{\code{\link{gof}}, \code{\link{gof.ergmm}}, and \code{\link[statnet]{gof.formula}}
         return an object of class \code{gofobject}.
       This is a list of the tables of statistics and posterior
       ``\eqn{p}-values".
       This is typically plotted using \code{\link{plot.gofobject}}.
      }

}
\seealso{ergmm, network, rergm, summary.ergmm, plot.gofobject}
\examples{
\dontrun{
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
# Fit the two-dimensional latent social space model 
#
# This may take a few minutes ...
#
samp.fit <- ergmm(samplike ~ latent(k=2), burnin=10000,
                  MCMCsamplesize=2000, interval=30)
#
# Posterior Predictive Checks
gofsamplike <- gof.ergmm(samp.fit, GOF=~idegree + distance)
gofsamplike
#
# Place both on the same page
# with nice margins
#
par(mfrow=c(1,2))
par(oma=c(0.5,2,1,0.5))
#
plot(gofsamplike)
#
# And now the odds 
#
plot(gofsamplike, plotodds=TRUE)
#
# Using Sampson's Monk data, lets 
# fit the two-dimensional clustered latent social space model 
#
# The ngroups parameter fits 3 groups
#
# This may take a few minutes ...
#
samp.fit <- ergmm(samplike ~ latentcluster(k=2, ngroups=3), burnin=10000,
                  MCMCsamplesize=2000, interval=30)
#
# Posterior Predictive Checks
gofsamplike <- gof.ergmm(samp.fit, GOF=~idegree + distance)
gofsamplike
#
plot(gofsamplike)
#
# And now the odds 
#
plot(gofsamplike, plotodds=TRUE)
}
}
\keyword{models}
