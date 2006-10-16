pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# This requires the package "statnet"
# See http://www.csde.washngton.edu/statnet
#
library(statnet)
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
pause()
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
pause()
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
