ergmm.getmodel.latent <- function (trms, g, drop=TRUE, expanded=FALSE) 
{
    v <- attr(trms, "variables")
    degindexcount <- 1
    degindex <- 0
    oldbottom <- 0
    oldtop <- 0
    newtop <- 0
    if (length(v) < 3) 
        stop(paste("No model specified for network ", trms[[2]]))
    ######Changed next line!
    m <- structure(list(node.attrib = NULL, coef.names = NULL,
          options = NULL, networkstats.0 = NULL, degreeinfo = NULL,
          latent=FALSE, cluster=FALSE),
        class = "model.ergmm")
    for (i in 3:length(v)) {
        if (is.call(v[[i]])) {
            v[[i]][[1]] <- as.name(paste("InitErgmm.", v[[i]][[1]], 
                sep = ""))
            for (j in length(v[[i]]):1) {
                v[[i]][[j + 2]] <- v[[i]][[j]]
                names(v[[i]])[j + 2] <- names(v[[i]])[j]
            }
            degindexcount <- degindexcount + 1
            #Note that subsetting has been changed -- before, it was
            #attempted without checking v[[i]][[2]] for length; this
            #caused problems for edgecov/dyadcov.  -CTB
            if ((length(v[[i]][[2]])>1)&&(v[[i]][[2]][[1]] == ":")) {
                degindexcount <- degindexcount + v[[i]][[2]][[3]] - 
                  v[[i]][[2]][[2]]
            }
        } else {
            degindexcount <- degindexcount + 1
            v[[i]] <- call(paste("InitErgmm.", v[[i]], sep = ""))
        }
        v[[i]][[2]] <- g
        names(v[[i]])[2] <-  ""
        v[[i]][[3]] <- m
        v[[i]][[length(v[[i]])+1]] <- drop
        m <- eval(v[[i]], sys.parent())
    }
    m$degreeinfo <- list(degindex = degindex, oldbottom = oldbottom, 
        oldtop = oldtop, newtop = newtop)
    m
}
