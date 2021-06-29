#' Maximum normalization
#' 
#' @description Performs the maximum normalization dividing every element of the 
#' vector with the maximum absolute value.
#' 
#' @param x Numeric vector.
#' 
#' @return Normalized numeric vector.
#' 
#' @examples{
#' max_normalization(c(1,2,3,4))
#' }
#'
#' @export
max_normalization <- function(x=numeric(0)){
  if(!is.numeric(x)){
    stop('"x" must be numeric.')
  }
  if(length(x)==0){
    stop('Empty vector.')
  }
  return(x/max(abs(x)))
}


#' Negative maximum normalization
#' 
#' @description Performs the negative maximum normalization on a non-negative 
#' numeric vector. The resulting vector is both reversed and normalized.
#' 
#' @param x Non-negative numeric vector.
#' 
#' @return Normalized numeric vector.
#' 
#' @details
#' For negative maximum normalization of a vector, difference between the 
#' maximum and minimal value of the vector is added tothe inverse value of every 
#' element and the result is divided with the original maximum value. 
#' Used formula: `(max(x)-min(x)-x)/max(x)`.
#' 
#' If the negative maximum normalization is applied twice to a 
#' non-negative vector, obtained result is the same as if the maximum 
#' normalization was applied once.
#' 
#' @examples{
#' neg_max_normalization(c(1,2,3,4))
#' neg_max_normalization(neg_max_normalization(c(1,2,3,4)))
#' }
#'
#' @export
neg_max_normalization <- function(x=numeric(0)){
  if(!is.numeric(x)){
    stop('"x" must be numeric.')
  }
  if(length(x)==0){
    stop('Empty vector.')
  }
  if(any(x<0)){
    stop('No negative values allowed')
  }
  return((max(x)-x+min(x))/max(x))
}
