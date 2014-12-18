#  File R/zzz.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

# Create a place for package-wide global variables.
.latentnetEnv <- new.env()

.onLoad <- function(lib, pkg){
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  .latentnetEnv$path.to.me <- file_path_as_absolute(lib)
  .latentnetEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
}

.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("latentnet",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}
