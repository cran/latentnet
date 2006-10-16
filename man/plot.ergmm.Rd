\name{plot.ergmm}
\alias{plot.ergmm}
\title{Plotting Method for class ERGMM}
\usage{
\method{plot}{ergmm}(x, \dots, mle=FALSE, comp.mat = NULL,
          label = NULL, label.col = "black",
          xlab, ylab, main, label.cex = 0.8, edge.lwd = 1,
          edge.col=1, al = 0.1,
          contours=0, density=FALSE, only.subdens = FALSE, 
          drawarrows=FALSE,
          contour.color=1, plotgraph=FALSE, pie = FALSE, piesize=0.07,
          vertex.col=1, vertex.pch=19, vertex.cex=2,
          mycol=c("black","red","green","blue","cyan",
                  "magenta","orange","yellow","purple"),
          mypch=15:19, mycex=2:10)

}
\description{
  \code{\link{plot.ergmm}} is the plotting method for \code{\link{ergmm}} objects.  For
  latent models, this plots the minimum Kullback-Leibler positions by 
  default.  The maximum likelihood positions
  can be used instead, or pie charts of the posterior
  probabilities of cluster membership can be shown.  It also plots the
  posterior density of the nodes or a series of contour plots, one for
  each node. See \code{\link{ergmm}} for more information on how to fit these
  models. Some plotting options depend on the \code{KernSmooth} package.
}
\arguments{
  \item{x}{an R object of class \code{\link{ergmm}}.
    See documentation for \code{\link{ergmm}}.}

  \item{mle}{Plots the network using the MLE of the positions for latent models.}

  \item{pie}{For latent clustering models, each node is drawn as a pie
    chart representing the probabilities of cluster membership.}
  \item{piesize}{The size of the pie charts.}
  
  \item{contours}{Plots a contours by contours array of the network with
    one contour per network corresponding to the posterior distribution of
    each of the nodes.}
 \item{contour.color}{Color of the contour lines.}

  \item{density}{If density=TRUE, plots the density of the posterior
    position of the nodes.  If density=c(nr,nc), plots a nr by nc array
    of density estimates for each cluster.}
  \item{only.subdens}{If density=c(nr,nc), only plots the densities of
    the clusters, not the overall density.}
 \item{drawarrows}{If density=TRUE, draws the ties on the density plot.}
 \item{plotgraph}{If density=c(nr,nc), a plot of the network is also shown.}
  
 \item{comp.mat}{For latent models, the positions are Procrustes
   transformed to look like comp.mat.}
 \item{label}{A vector of the same length as the number of nodes
   containing the labels of the nodes.}
 \item{label.col}{The color to be used for plotting the labels.}
 \item{label.cex}{The size of the node labels.}

 \item{xlab}{Title for the x axis.}
 \item{ylab}{Title for the y axis.}
 \item{main}{The main title for the graph.}

 \item{edge.lwd}{The line width for the arrows between nodes.}
 \item{edge.col}{The color of the arrows between nodes.}
 \item{al}{The length of the arrow heads.}

 \item{vertex.col}{The color of the nodes as defined by \code{mycol}.  Can be
   specified as an attribute of the network used in the model.}
 \item{vertex.pch}{The plotting character of the nodes as defined by
   \code{mypch}.  Can be specified as an attribute of the network used in
   the model.}
 \item{vertex.cex}{The size of the nodes as defined by \code{mycex}.  Can be
   specified as an attribute of the network used in the model.}
 
 \item{mycol}{Vector of colors to be used.  Defaults to:
   c("black","red","green","blue","cyan",
   "magenta","orange","yellow","purple")}
 \item{mypch}{Vector of plotting characters to be used.  Defaults to: }
 \item{mycex}{Vector of character expansion values.}
 \item{\dots}{Other optional arguments to be used by the plot function.}
}

\details{
  Plots the results of an ergmm fit.
  
  More information can be found by looking at the documentation of
  \code{\link{ergmm}}.
}
\value{NULL

}
\seealso{ergmm, network, plot.network, plot, ergmm.add.contours}
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
}
\keyword{graphs}
