\name{mcmc.diagnostics.ergm}
\alias{mcmc.diagnostics.ergm}
\alias{mcmc.diagnostics}
\alias{ergm.raftery.diag}
\title{Conduct MCMC diagnostics on an ergm fit}
\description{
  This function creates simple diagnostic plots for the MCMC sampled
  statistics produced from a fit. It also prints the Raftery-Lewis 
  diagnostics, indicates if they are sufficient,
  and suggests the run length required.
}
\usage{
\method{mcmc.diagnostics}{ergm}(object, sample = "sample",
                      smooth=TRUE, r = 0.0125, digits = 6,
                      maxplot = 1000, verbose = TRUE, 
                      mcmc.title = "Summary of MCMC samples", \dots)
ergm.raftery.diag(data, q = 0.025, rmargin = 0.005, s = 0.95, converge.eps = 0.001)
}
\arguments{
  \item{object}{An object.  See documentation for \code{\link{ergm}}.}
  \item{sample}{The component of \code{object} on which the
    diagnosis is based. The two usual ones are \code{thetasample} from the 
    auxiliary sample of the natural parameter and \code{sample} the (default)
    sample of the sufficient statistics from the model.}
  \item{smooth}{Draw a smooth line through trace plots}
  \item{r}{Percentile of the distribution to estimate}
  \item{digits}{Number of digits to print}
  \item{maxplot}{Maximum number of statistics to plot}
  \item{data}{an 'mcmc' object, typically the component of
        \code{object} on which the diagnosis is based.}
  \item{q}{the quantile to be estimated.}
  \item{rmargin}{the desired margin of error of the estimate.}
  \item{s}{the probability of obtaining an estimate in the interval (q-r,q+r).}
  \item{converge.eps}{Precision required for estimate of time to convergence.}
  \item{verbose}{If this is \code{TRUE}, print out more information about the MCMC runs including lag correlations.}
  \item{mcmc.title}{Figure title for the diagnostic plots.}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\details{
The plots produced are a trace of the sampled output
and a density estimate for each variable in the chain.

The Raftery-Lewis diagnostic is a run length control diagnostic based on a
criterion of accuracy of estimation of the quantile q. It
is intended for use on a short pilot run of a Markov chain.
The number of iterations required to estimate the quantile
q to within an accuracy of +/- r with probability p is
calculated. Separate calculations are performed for each
variable within each chain.

In fact, an \code{object} contains the matrix of
statistics from the MCMC run as component \code{$sample}.
This matrix is actually an object of class \code{mcmc} and
can be used directly in the \code{CODA} package to assess MCMC
convergence. \emph{Hence all MCMC diagnostic methods available
in \code{coda} are available directly.} See the examples and
%\url{http://www.mrc-bsu.cam.ac.uk/bugs/classic/coda04/readme.shtml}.
the \code{coda} package.

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}
\value{\code{mcmc.diagnostics.ergm} returns a table of Raftery-Lewis diagnostics.
      }

}
\seealso{ergm, network, coda}
\examples{
#
data(sampson)
#
# test the mcmc.diagnostics function
#
gest <- ergm(samplike ~ latent(k=2))
summary(gest)
#
# Plot the probabilities first
#
mcmc.diagnostics(gest)
#
# Use coda directly
#
library(coda)
#
plot(gest$sample, ask=FALSE)
ergm.raftery.diag(gest$sample, r=0.1)
#
# A full range of diagnostics are available 
# using codamenu()
#
}
\keyword{graphs}
