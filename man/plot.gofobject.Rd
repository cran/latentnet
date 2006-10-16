\name{plot.gofobject}
\alias{plot.gofobject}
\title{Plot Goodness-of-Fit Diagnostics on a Latent Space Random Graph Model}
\description{
  \code{\link{plot.gofobject}} plots diagnostics such as the 
  degree distribution, geodesic distances, shared partner distributions,
  and reachability for the posterior predictive goodness-of-fit of Latent Space random graph
  models.  See \code{\link{ergmm}} for more information on these models.
}
\usage{
\method{plot}{gofobject}(x, \dots,
         cex.axis=0.7, plotodds=FALSE,
         main = "Goodness-of-fit diagnostics",
         normalize.reachability=FALSE,
         verbose=FALSE)
}
\arguments{
  \item{x}{an object of class \code{gofobject}, 
       typically produced by the
       \code{\link{gof.ergmm}} or \code{\link[statnet]{gof.formula}} functions.
       See the documentation for these.}
  \item{cex.axis}{Character expansion of the axis labels relative to that for
   the plot.}
  \item{plotodds}{Plot the odds of a dyad having given
  characteristics (e.g., reachability, minimum geodesic
  distance, shared partners). This is an alternative to the
  probability of a dyad having the same property.}
  \item{main}{Title for the goodness-of-fit plots.}
  \item{normalize.reachability}{Should the reachability proportion be 
   normalized to make it more comparable with the other geodesic distance 
   proportions.}
  \item{verbose}{Provide verbose information on the progress of the
    plotting.}
  \item{\dots}{Additional arguments, to be passed to the plot function.}
}

\details{
  \code{\link{gof.ergmm}} produces a sample of networks randomly drawn from the specified model.
  This function produces a plot of the summary measures.
}
\value{none}
\seealso{ergm, network, rergm.ergmm, summary.ergmm}
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
samp.fit <- ergmm(samplike ~ latent(k=2), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# Posterior Predictive Checks
gofsamplike <- gof.ergmm(samp.fit)
gofsamplike
#
# Place all three on the same page
# with nice margins
#
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
#
plot(gofsamplike)
#
# And now the odds 
#
plot(gofsamplike, plotodds=TRUE)
#
# Using Sampson's Monk data, lets fit a latent clustering model
#
samp.fit <- ergmm(samplike ~ latentcluster(k=2, ngroups=3), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# Posterior Predictive Checks
gofsamplike <- gof.ergmm(samp.fit)
gofsamplike
#
# Place all three on the same page
# with nice margins
#
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
#
plot(gofsamplike)
#
# And now the odds 
#
plot(gofsamplike, plotodds=TRUE)
}
}
\keyword{models}
