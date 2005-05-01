ergmm.gettermnames.latent <- function (trms) 
{
    v <- attr(trms, "variables")
    vnames <- character(length=length(v)-2)
    if (length(v) < 3) 
        stop(paste("No model specified for network ", trms[[2]]))

    for (i in 3:length(v)) {
      vnames[i-2] <- as.character(v[[i]])[1]
    }
    vnames
}
