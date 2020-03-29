
Pmat_product<- function(x, P, type = NULL){
  if(is.null(type)){type = attr(P,"type")}
  ## check x is matrix
  if(ncol(x)==1){return(Pmat_product(x = as.numeric(x), P = as.matrix(P), type = type))}
  ### check P is vector
  if(is.matrix(P)){
    if(ncol(P)!=1){warning("Only 1 transformation is allowed when x is a matrix.")}
    else{
      P <- as.numeric(P)}}
  ### check length
  if(nrow(x)!=length(P)){warning("x and P should have matching row length.")}
  ## product
  switch(type,
         "permutation" ={x <- matrix(x[as.numeric(P),],ncol = ncol(x))},
         "signflip" = {x <- matrix(x*as.numeric(P),ncol = ncol(x))})
  
  return(x)
}