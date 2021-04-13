#' Calculates the mean and standard deviation for correlation by random resampling
#'
#'@description For the given data matrix, the mean and standard deviation of correlation is 
#'calculated for every pair of variables. The mean and standard deviation are calculated
#'based on correlation values obtained after random subsampling.
#'
#'@params data Data matrix or data frame with samples in rows and variables in columns.
#'@params size_of_subsample Size of subsample taken in every repetition.
#'@params replace Defines whether the resampling is done with replacement. Defaults to FALSE.
#'@params n_repetitions Number of repetitions of subsampling. For large 'data' matrices, large
#'number of repetitions can take a long time.
#'@params zero_precisionC The zero precision passed to the C++ function. All numbers with absolute
#'value smaller than 'zero_precisionC' will be set to zero. Defaults to 10^(-6).
#'
#'@returns A data frame with names of variables in first two columns. In third columns is mean
#'and in the fourth column is standard deviation for the given pair of variables. In the rows, 
#'all pairwise correlations are listed.
#'
#'@details For 'n_repetitions' times, a random subsample is taken based on parameters 
#''size_of_subsample' and 'replace'. Then for every pair of variables a correlation is calculated
#'for a random subsample. At the end of repetitions, the mean and standard deviation of correlation
#'are calculated for every pair of variables.
#'
#'@export
cor_resamplingStats <- function(data, size_of_subsample, replace=FALSE, n_repetitions=100, zero_precisionC=0.000001){
  if(!inherits(data,"matrix")&!inherits(data,"data.frame")){
    stop("'data' must be either matrix or data frame.")
  }
  if(is.null(colnames(data))){
    variable_names <- 1:ncol(data)
  }else{
    variable_names <- colnames(data)
  }
  data <- as.matrix(data)
  if(any(is.na(data))){
    stop("NA values in 'data' are not allowed.")
  }
  if(!all(is.numeric(data))){
    stop("Values in 'data' must be numeric.")
  }
  #
  cor_summary <- cor_meanSd(data, size_of_subsample, replace, n_repetitions, zero_precisionC)
  first_variables <- unlist(sapply(1:(length(variable_names)-1), function(x) variable_names[1:x]))
  second_variables <- rep(variable_names[2:length(variable_names)], 1:(length(variable_names)-1))
  output_df <- data.frame(Var1=first_variables, Var2=second_variables, MeanCor=cor_summary[,1], SdCor=cor_summary[,2])
  return(output_df) #order of pairs of variables is the same as in upper triangular correlation matrix: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4) (1,5)...
}