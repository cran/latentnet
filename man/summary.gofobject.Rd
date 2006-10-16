\name{summary.gofobject}
\alias{summary.gofobject}
\alias{print.gofobject}
\title{Summaries the Goodness-of-Fit Diagnostics on a Latent Space Graph Model}
\description{
  \code{\link{summary.gofobject}} summaries the diagnostics such as the 
  degree distribution, geodesic distances, shared partner distributions,
  and reachability for the goodness-of-fit of Latent Space random graph
  models.  See \code{\link{ergmm}} for more information on these models.
}
\usage{
\method{summary}{gofobject}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{gofobject}, 
       typically produced by the
       \code{\link{gof.ergmm}} or \code{\link[statnet]{gof.formula}} functions.
       See the documentation for these.}
  \item{\dots}{Additional arguments, to be passed to the plot function.}
}

\details{
  \code{\link{gof.ergmm}} produces a sample of networks randomly drawn from the specified model.
  This function produces a print out the summary measures.
}
\value{none}
}
\seealso{gof.ergmm, ergmm, network, rergm}
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
summary(gofsamplike)
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
summary(gofsamplike)
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
\keyword{graphs}
