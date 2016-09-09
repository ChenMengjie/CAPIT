Thresholding <- function(element, gamma, selection = "hard"){
  hard.e <- element*as.integer(abs(element) > gamma)
  if(is.vector(element)){
    soft.e <- rep(0, length(element))
    for(i in 1:length(element)){
      soft.e[i] <- sign(element[i])*max(abs(element[i])-gamma, 0)
    }
  }
  else{
    soft.e <- sign(element)*max(abs(element)-gamma, 0)
  }
  res <- switch(selection, hard = hard.e, soft = soft.e) 
  return(res)
}
