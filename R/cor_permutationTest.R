#' Correlation permutation test
#'
#' @description For the given data matrix, the permutation test of correlation is calculated
#' for every pair of variables.
#'
#' @param data Data matrix or data frame with samples in rows and variables in columns.
#' @param n_repetitions Number of repetitions of subsampling. For large `data` matrices, large
#' number of repetitions can take a long time.
#' @param alternative Alternative hypothesis for the permutation test. Takes one of the values: 
#' `"less"`, `"greater"`, `"two_sided"` or `"two_sided_signed"`. See details.
#' @param zero_precisionC The zero precision passed to the C++ function. All numbers with absolute
#' value smaller than `zero_precisionC` will be set to zero. Defaults to 10^(-6).
#'
#' @return A data frame with names of variables in first two columns and p-value of correlation
#' test between the two variables in the third column. In the rows, all pairwise correlations are 
#' listed.
#'
#' @details 
#' For every pair of variables, for `n_repetitions` times, a permutation of the samples is 
#' randomly chosen and the correlation between the original samples for variable i and permuted 
#' samples for variable j is calculated and compared to the true correlation coefficient. The p-value
#' is obtained as percentage of times when the "permuted" correlation was more significant that the
#' "true" correlation. Significance is determined based on the `alternative` parameter.
#' If the `alternative` is
#' \itemize{
#'    \item{`"less"`}{We count number of times "permuted" correlation is smaller 
#'    than the "true" correlation.}
#'    \item{`"greater"`}{We count number of times "permuted" correlation 
#'    is greater than the "true" correlation.}
#'    \item{`"two_sided"`}{We count number of times absolute value of "permuted" 
#'    correlation is greater than the absolute value of the "true" correlation.}
#'    \item{`"two_sided_signed"`}{We count when "permuted" correlation is greater 
#'    than the "true" correlation for positive "true" correlation or when "permuted" 
#'    correlation is smaller than the "true" correlation for negative "true" correlation.}
#' }
#'
#' @export
cor_permutationTest <- function(data, n_repetitions=100, alternative="two_sided", zero_precisionC=0.000001){
  if(!is.character(alternative)){
    stop(paste0("'alternative'",' must have one of the values "less", "greater", "two_sided" or "two_sided_signed".'))
  }
  if(!(alternative %in% c("less", "greater", "two_sided", "two_sided_signed"))){
    stop(paste0("'alternative'",' must have one of the values "less", "greater", "two_sided" or "two_sided_signed".'))
  }
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
  cor_pval <- cor_prTest(data, n_repetitions, alternative, zero_precisionC)
  first_variables <- unlist(sapply(1:(length(variable_names)-1), function(x) variable_names[1:x]))
  second_variables <- rep(variable_names[2:length(variable_names)], 1:(length(variable_names)-1))
  output_df <- data.frame(Var1=first_variables, Var2=second_variables, prTest_pval=cor_pval)
  return(output_df) #order of pairs of variables is the same as in upper triangular correlation matrix: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4) (1,5)...
}
