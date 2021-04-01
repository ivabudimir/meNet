#@pars data has variables in columns and samples in rows
#'@export
cor_permutationTest <- function(data, n_repetitions=100, alternative="two_sided", zero_precisionC=0.000001){
  if(!is.character(alternative)){
    stop('"alternative" must have one of the values "less", "greater", "two_sided" or "two_sided_signed".')
  }
  if(!(alternative %in% c("less", "greater", "two_sided", "two_sided_signed"))){
    stop('"alternative" must have one of the values "less", "greater", "two_sided" or "two_sided_signed".')
  }
  if(!inherits(data,"matrix")&!inherits(data,"data.frame")){
    stop('"data" must be either matrix or data frame.')
  }
  if(is.null(colnames(data))){
    variable_names <- 1:ncol(data)
  }else{
    variable_names <- colnames(data)
  }
  data <- as.matrix(data)
  if(any(is.na(data))){
    stop('NA values in "data" are not allowed.')
  }
  if(!all(is.numeric(data))){
    stop('Values in "data" must be numeric.')
  }
  #
  cor_pval <- cor_prTest(data, n_repetitions, alternative, zero_precisionC)
  first_variables <- unlist(sapply(1:(length(variable_names)-1), function(x) variable_names[1:x]))
  second_variables <- rep(variable_names[2:length(variable_names)], 1:(length(variable_names)-1))
  output_df <- data.frame(Var1=first_variables, Var2=second_variables, prTest_pval=cor_pval)
  return(output_df) #order of pairs of variables is the same as in upper triangular correlation matrix: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4) (1,5)...
}
