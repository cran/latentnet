#  File R/ergmm.permutation.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
ergmm.permutation <- function(n)
{
  if(n==2)
    return(matrix(c(1,2,2,1),2,2))
  temp <- ergmm.permutation(n-1)
  for(i in 0:(n-1))
  {
    if(i==0)
      temp2 <- cbind(temp,n)
    else{
      if (i==(n-1))
        temp2 <- rbind(temp2,cbind(n,temp))
      else
        temp2 <- rbind(temp2,cbind(temp[,seq_len(i)],n,temp[,(i+1):(n-1)]))
    }
  }
  colnames(temp2)<-seq_len(n)
  return(temp2)
}

which.perm.nearest<-function(K.ref,K){
  perms<-ergmm.permutation(max(c(K.ref,K)))
  order(as.numeric(perms[which.min(apply(perms[,K],1,function(K.perm)sum(K.perm!=K.ref))),]))
}
