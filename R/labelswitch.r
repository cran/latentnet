ergmm.labelswitch <- function(Ki,per.to)
{
  Ki.new <- Ki
  for(i in 1:length(per.to))
    Ki.new[Ki == per.to[i]] <- i
  return(Ki.new)
}
