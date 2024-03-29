% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.ergmm.R
\name{summary.ergmm}
\alias{summary.ergmm}
\alias{print.summary.ergmm}
\alias{summary.ergmm.object}
\alias{bic.ergmm}
\title{ERGMM Fit Summaries}
\usage{
\method{summary}{ergmm}(
  object,
  point.est = c(if (!is.null(object[["mle"]])) "mle", if (!is.null(object[["sample"]]))
    c("pmean", "mkl")),
  quantiles = c(0.025, 0.975),
  se = "mle" \%in\% point.est,
  bic.eff.obs = c("ties", "dyads", "actors"),
  ...
)
}
\arguments{
\item{object}{An \code{\link[=ergmm.object]{ergmm}} object to be summarized.}

\item{point.est}{Point estimates to compute: a character vector with some
subset of \code{"mle"}, \code{"pmean"}, \code{"mkl"}, and \code{"pmode"}.
Defaults to a concatenation of \code{"mle"} (if fit), \code{"pmean"}, and
\code{"mkl"} (if MCMC was run).}

\item{quantiles}{Posterior quantiles (credible intervals) to compute.}

\item{se}{Whether to compute standard errors. Defaults to \code{TRUE} if MLE
was fit.}

\item{\dots}{Additional arguments.}

\item{eff.obs, bic.eff.obs}{What effective sample size to use for BIC
calculation?
\describe{
\item{\code{"ties"}}{the number of non-missing ties in the network. This is the approach recommended by Handcock et al. (2007) and the default. Not well-defined for valued networks.}
\item{\code{"dyads"}}{the number of non-missing dyads (potential ties) in the network.}
\item{\code{"actors"}}{the number of actors in the network. The default prior to 2.7.0.}
\item{a number}{to specify a specific sample size.}
\item{\code{NULL}}{Don't compute the BIC at all. Mostly for internal use.}
}}
}
\value{
For \code{summary}, an object of class
\code{\link[=summary.ergmm.object]{summary.ergmm}}. A print method is
available.

The BICs are available as the element "bic" of the object returned.

\code{bic.ergmm} returns the BIC for the model directly.
}
\description{
\code{summary.ergmm} prodcues a summary of an
\code{\link[=ergmm.object]{ergmm}} object, including point estimates,
standard errors, and BIC calculation.
}
\details{
Note that BIC computed for the random effects models uses the same
formualtion as Handcock et al., so it is likely correct, but has not been
peer-reviewed.

This BIC can be (reasonably) safely used to select the number of clusters or
which fixed effects to include in the model. It is not clear whether it is
appropriate to use this BIC to select the dimension of latent space and
whether or not to include random actor effects. These considerations are
independent of the bug described below.

Prior to version 2.7.0, there was a bug in BIC calculation that used \eqn{p
+ n(d+r+s)} as the number of parameters in the likelihood (where \eqn{p} is
the number of fixed effects, \eqn{n} the number of actors, \eqn{d}, the
latent space dimension, and \eqn{r} and \eqn{s} indicators of presence of
sender and receiver (or sociality) effects). This value should have been
just \eqn{p}.

The following applications could have produced different results:
\itemize{
\item{Using the BIC to select latent space dimension.}
\item{Using the BIC to decide whether or not to include random effects.}
}
The following applications could not (i.e., would be off by a constant):
\itemize{
\item{Using the BIC to select the number of clusters.}
\item{Using the BIC to select the fixed effects to be used.}
}
}
\examples{

\donttest{
data(sampson)
# Fit the model for cluster sizes 1 through 4:
fits<-list(
           ergmm(samplike~euclidean(d=2,G=1)),
           ergmm(samplike~euclidean(d=2,G=2)),
           ergmm(samplike~euclidean(d=2,G=3)),
           ergmm(samplike~euclidean(d=2,G=4))
           )

\dontrun{
# Optionally, plot all fits.
lapply(fits,plot)
}

# Compute the BICs for the fits and plot them:
(bics<-reshape(
    as.data.frame(t(sapply(fits,
                           function(x)c(G=x$model$G,unlist(bic.ergmm(x))[c("Y","Z","overall")])))),
    list(c("Y","Z","overall")),idvar="G",v.names="BIC",timevar="Component",
    times=c("likelihood","clustering","overall"),direction="long"
    ))

with(bics,interaction.plot(G,Component,BIC,type="b",xlab="Clusters", ylab="BIC"))

# Summarize and plot whichever fit has the lowest overall BIC:
bestG<-with(bics[bics$Component=="overall",],G[which.min(BIC)])
summary(fits[[bestG]])
plot(fits[[bestG]])
}

}
\references{
Chris Fraley and Adrian E. Raftery (2002). \emph{Model-based
clustering, discriminant analysis, and density estimation}. Journal of the
American Statistical Association, 97(458), 611-631.

Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
\emph{Model-Based Clustering for Social Networks}.  Journal of the Royal
Statistical Society: Series A (Statistics in Society), 170(2), 301-354.
}
\seealso{
\code{\link{ergmm.object}}, \code{\link{ergmm}}
}
\keyword{graphs}
\keyword{models}
\keyword{print}
