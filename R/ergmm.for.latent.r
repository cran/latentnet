ergmm <- function(formula, theta0=NULL, 
        burnin=1000, MCMCsamplesize=1000, interval=10,
        latent.control=list(maxit=40,penalty.sigma=c(10,0.5),MLEonly=FALSE),
        returnMCMCstats=TRUE, randseed=NULL,
        verbose=FALSE, ...)
{
    current.warn <- options()$warn
    options(warn=0)
    statsmatrix <- 0
    Clist <- 0

    verb <- match(verbose,
      c("FALSE","TRUE", "very"), nomatch=1)-1

    trms <- ergmm.getterms.latent(formula)
    termnames <- ergmm.gettermnames.latent(trms)
    g <- try(as.network(eval(trms[[2]],sys.parent())))
    if(inherits(g,"try-error")){
     stop("Invalid network. Is the left-hand-side of the formula correct?")
    }

    m <- ergmm.getmodel.latent(trms, g)
    Clist <- ergmm.Cprepare.latent(g, m)

    if(is.null(randseed)){randseed <- sample(10000000, size=1)}
    set.seed(as.integer(randseed))

    notobserved <- get.network.attribute(g,"design")
    if(is.null(notobserved)){
     mClist <- list(heads=0, tails=0, nedges=0, dir=is.directed(g))
    }else{
     mClist <- ergmm.Cprepare.latent(notobserved, m)
     cat("Design matrix:\n")
     summary(notobserved)
    }
    
    v <- latent.wrapper(theta0, trms, g, m, Clist, mClist,
                        MCMCsamplesize, burnin, interval,
                        formula,
                        latent.control, verbose)
    theta.original <- theta0

    v$network <- g
    if(returnMCMCstats){
        endrun <- burnin+interval*(MCMCsamplesize-1)
        attr(v$sample, "mcpar") <- c(burnin+1, endrun, interval)
        attr(v$sample, "class") <- "mcmc"
     }else{
        v$sample <- NULL
     }
#
#   Final ordering of results
#
    v$interval <- interval
    v$theta.original <- theta.original
    options(warn=current.warn)
    v
}
