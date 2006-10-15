pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# Read in the set of ties from Sampson monastery data 
# These are the cumulative ties over three periods.
# A tie  from monk A to monk B exists if A nominated 
# B as one of his three best friends at any of the last
# three time points.
# 
data(sampson)
#
# Summarize the graph
#
summary(samplike)
#
# Look at the names of the monks (followed by there order in the original)
#
network.vertex.names(samplike)
#
# Look at the cohesive sub-groups designated by Sampson 1969.
#
samplike %v% "group"
#
# Fit the two-dimensional latent social space model 
# This may take a few minutes ...
#
pause()
#
samplike.fit <- ergmm(samplike ~ latent(k=2),
   burnin=50000, MCMCsamplesize=1000, interval=1000)
#
# Plot some diagnostics
#
pause()
mcmc.diagnostics(samplike.fit)
#
# Summarize the fit
#
summary(samplike.fit)
pause()
#
samplike.fit$Z.mkl	#to list the minimum Kullback-Leibler positions
samplike.fit$Z.mle	#to list the MLE positions
pause()
#
# Create the plot symbols for the groups
#
oneL <- samplike %v% "group"
oneL[oneL=="Turks"] <- "T"
oneL[oneL=="outcasts"] <- "O"
oneL[oneL=="loyal"] <- "L"
oneL[c(1,7,15)] <- "W"
oneL
#
# Create colors
oneLcolors <- c("red","blue","black","green")[match(oneL,c("T","O","L","W"))]
#
# Plot the MLE positions
#
plot(samplike.fit,label=oneL,vertex.col=oneLcolors,mle=T,main="MLE positions")
title(sub="Color represents the estimated groups; Labels the Sampson's groups")
pause()
#
# Plot the MKL positions
#
plot(samplike.fit,label=oneL,vertex.col=oneLcolors)
title(sub="Color represents the estimated groups; Labels the Sampson's groups")
pause()
#
par(mfrow=c(2,2))
#
# Plot 24 versions of the positions in latent space to
# see how uncertain they are
#
number.to.plot <- 24
for(i in 1:number.to.plot){
 isamp <- 1 + (i-1)*round(dim(samplike.fit$Z)[3]/number.to.plot)
 aaa <- plot(samplike, vertex.col=oneLcolors, label="",
  arrowhead.length = 0.1, vertex.cex=2,
  coord=samplike.fit$Z[,,isamp],
  main=paste("Sample number",isamp))
 if(4*trunc(i/4)==i){pause()}
}
