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
        temp2 <- rbind(temp2,cbind(temp[,1:i],n,temp[,(i+1):(n-1)]))
    }
  }
  return(temp2)
}
