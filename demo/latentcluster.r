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
# Fit the two-dimensional clustered latent social space model 
#
# The ngroups parameter fits 3 groups
#
# This may take a few minutes ...
pause()
#
samplike.fit3 <- ergmm(samplike ~ latentcluster(k=2, ngroups=3,
   Sigprior = qchisq(0.05,2), muSigprior = 2, alphaprior=2, b.prior.sd=2),
   burnin=5000, MCMCsamplesize=1000, interval=30)
#
# Look at the goodness-of-fit of the 4 group model
#
samplike.fit3$BIC
#
# Plot some diagnostics
#
pause()
mcmc.diagnostics(samplike.fit3)
#
# Let's try the fit of the 4 group model
#
# This may take a few minutes ...
#
pause()
samplike.fit4 <- ergmm(samplike ~ latentcluster(k=2, ngroups=4,
   Sigprior = qchisq(0.05,2), muSigprior = 2, alphaprior=2, b.prior.sd=2),
   burnin=5000, MCMCsamplesize=1000, interval=30)
#
# Look at the goodness-of-fit of the 4 group model
# The BIC (Bayesian Information Criterion) is a form
# of (negative) deviance penalized for the complexity of the model
# (as measured by the number of dimensions and the number of groups).
#
samplike.fit4$BIC
#
# The better model has the higher BIC, so we will go with the 
# 3 group model
pause()
#
# Summarize the fit
#
summary(samplike.fit3)
pause()
#
# Print out the probabilities of group membership 
# for each monk
#
# First for the 4 group model
#
round(samplike.fit4$qig,2)
#
# Now for the 3 group model
#
round(samplike.fit3$qig,2)
#
samplike.fit3$Z.mkl	#to list the minimum Kullback-Leibler positions
samplike.fit3$class	#to list the vector of posterior modal classes/clusters
samplike.fit3$Z.mle	#to list the MLE positions
samplike.fit3$Ki.mle	#to list the vector of maximum likelihood classes/clusters
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
# Plot the MLE positions
#
plot(samplike.fit3,label=oneL,mle=T,main="MLE positions")
title(sub="Color represents the estimated groups; Labels the Sampson's groups")
pause()
#
# Plot the MKL positions
#
plot(samplike.fit3,label=oneL)
title(sub="Color represents the estimated groups; Labels the Sampson's groups")
pause()
#
# Plot the MKL positions as pie charts
#
plot(samplike.fit3,label=oneL,pie=T)
title(sub="Color represents the estimated groups; Labels the Sampson's groups")
pause()
# Plot some densities
plot(samplike.fit3,density=c(2,2))
pause()

par(mfrow=c(2,2))
#
# Plot 24 versions of the positions in latent space to
# see how uncertain they are
#
number.to.plot <- 24
for(i in 1:number.to.plot){
 isamp <- 1 + (i-1)*round(dim(samplike.fit3$Z)[3]/number.to.plot)
 aaa <- plot(samplike, vertex.col=samplike.fit3$Ki[isamp,], label="",
  arrowhead.length = 0.1, vertex.cex=2,
  coord=samplike.fit3$Z[,,isamp],
  main=paste("Sample number",isamp))
 if(4*trunc(i/4)==i){pause()}
}
