#'@export
max_normalization  <- function(x=numeric(0)){
  if(!is.numeric(x)){
    stop('"x" must be numeric.')
  }
  if(length(x)==0){
    stop('Empty vector.')
  }
  return(x/max(abs(x)))
}


# this function squared gives "max_normalization"
#'@export
neg_max_normalization  <- function(x=numeric(0), add_part=0.1){
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