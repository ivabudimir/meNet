#' Builds a CpG network where edges are based on correlation
#' 
#' @description Different options are available to determine which correlations 
#' are significant, such as a constant threshold for correlation, standard 
#' deviation of bootstrapping or p-value of the permutation test. More than one 
#' method can be used simultaneously. Significant correlations determine the 
#' edged of the network.  For better explanation of the significant criteria, 
#' read Details.
#' 
#' @param cor_matrix Correlation matrix of CpG sites.
#' @param data Data frame with CpGs in columns. Variables in rows are used to
#' calculate "cor_matrix".
#' @param cg_ids List of CpGs for which we reconstruct the network. If names of
#' CpGs are given as "cor_matrix" or "data" column names, "cg_ids" defines a subset
#' of CpGs which should be used in network. If omitted, all CpGs are used.
#' If CpG names are not given, "cg_ids" is a mandatory parameter which gives 
#' the names of CpGs.
#' @param cor_threshold Correlation threshold. Defaults to 0.2.
#' @param neg_cor_threshold Negative correlation threshold. This parameter is 
#' ignored if "cor_threshold" is not given. Defaults to NULL.
#' @param cor_stDev Threshold for the standard deviation of correlation.
#' Default to "cor_threshold"/3.
#' @param cor_alpha Significance level applied to the correlation permutation
#' test. Defaults to NULL.
#' @param n_repetitions Number of repetitions for resampling and/or for the 
#' correlation permutation test. Defaults to 100.
#' @param alternative Alternative hypothesis for the correlation permutation
#' test. Default to "two_sided".
#' @param normalization_function Normalization function applied to the weights
#' of the edges. By default, no normalization is applied and weights correspond
#' to correlation.
#' 
#' 
#' @return Weighted network as igraph object with CpGs as nodes and edges
#' representing significantly correlated pairs of CpGs. The weights of edges
#' represent (possibly normalized) correlation.
#' Isolated nodes are kept in the network.
#' 
#' @details A method for a significance of correlations is used if its parameters
#' are correctly specified. At the same time, more than one method can be used
#' in which case only correlations which are significant according to all methods
#' are kept.
#' For a threshold method to be used, "cor_threshold" has to have a numeric value.
#' Additionally, "neg_cor_threshold" can be specified. If only "cor_threshold"
#' is given, correlations with absolute value larger than "cor_threshold" are 
#' considered significant. If also "neg_cor_threshold" is given, correlations 
#' smaller than "neg_cor_threshold" or larger than "cor_threshold" are 
#' considered significant. This allows different penalization of negative
#' correlation values.
#' For a method based on standard deviation of bootstrapping to be used,
#' "cor_stDev" has to have a numeric value, "data" has to be provided and 
#' "n_repetitions" has to be correctly specified. Function
#' "meNet::cor_resamplingStats" is called for "data" with subsample size equal
#' to the number of rows in "data" and with replacement. Correlations for which
#' the calculated standard deviation of resampling is smaller than "cor_stDev"
#' are considered significant. This method tests variability of correlation.
#' For a method based on the permutation test to be used, "cor_alpha" has to 
#' have a numeric value, "data" has to be provided and parameters "n_repetitions"
#' and "alternative" have to be correctly specified. Parameter "alternative" has 
#' to have one of the values "less", "greater", "two_sided" or "two_sided_signed".
#' Function "meNet::cor_permutationTest" is called for "data". Correlations for 
#' which the p-value of the permutation test are smaller than "cor_alpha"
#' are considered significant.
#' 
#' @import igraph
#' 
#' @export
meNet_cor <- function(cor_matrix=NULL, data=NULL, cg_ids=NULL, cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=cor_threshold/3, cor_alpha=NULL,
                      n_repetitions=100, alternative="two_sided", normalization_fun=NULL){
  #
  matrix_preprocessing <- .check_corM_dataDF(cor_matrix, data)
  cor_matrix <- matrix_preprocessing[[1]]
  data <- matrix_preprocessing[[2]]
  #
  if(!is.null(cg_ids)){
    if(all(rownames(cor_matrix)==colnames(cor_matrix))){
      if(length(intersect(cg_ids, rownames(cor_matrix)))==0){
        warning('"cg_ids" elements are not represented in the "cor_matrix"/"data".')
        return(make_empty_graph(n=0,directed=FALSE))
      }
      cor_matrix <- cor_matrix[cg_ids, cg_ids]
      data <- data[,cg_ids]
    }else{
      tryCatch(
        expr = {
          rownames(cor_matrix) <- cg_ids
          colnames(cor_matrix) <- cg_ids
        },
        error = function(e){
          stop('If variable names are not given in the matrix, the length of "cg_ids" must be equal to the number of variables in "cor_matrix"/"data".')
        }
      )
    }
  }else if(any(rownames(cor_matrix)!=colnames(cor_matrix))){
    rownames(cor_matrix) <- colnames(cor_matrix)
    warning('"cor_matrix"/"data" column names will be used as node labels.')
  }
  # remove correlations below threshold (if given)
  adjacency_matrix <- cor_matrix
  diag(adjacency_matrix) <- 0
  if(is.numeric(cor_threshold)){
    if(is.numeric(neg_cor_threshold)){
      adjacency_matrix <- (adjacency_matrix>cor_threshold)*adjacency_matrix + (adjacency_matrix<neg_cor_threshold)*adjacency_matrix
    }else{
      adjacency_matrix <- (abs(adjacency_matrix)>cor_threshold)*adjacency_matrix
    }
  }
  if(sum(abs(adjacency_matrix))>0 & is.numeric(cor_stDev) & !is.null(data)){
    if(!is.numeric(n_repetitions)){
      stop('"n_repetitions" must be a positive integer.')
    }
    significance_matrix <- matrix(0, nrow=nrow(cor_matrix), ncol=nrow(cor_matrix))
    significance_matrix[upper.tri(significance_matrix,diag=FALSE)] <- cor_resamplingStats(data=as.matrix(data), size_of_subsample=nrow(data), replace=TRUE, n_repetitions=n_repetitions)[,4]
    adjacency_matrix <- (significance_matrix<cor_stDev)*adjacency_matrix
  }else if(!is.null(cor_stDev)&(sum(abs(adjacency_matrix))==0)){
    warning("Standard deviation significance testing wasn't performed.")
  }else if(!is.null(cor_stDev)){
    warning("Standard deviation significance testing wasn't performed. Check the parameters.")
  }
  if(sum(abs(adjacency_matrix))>0 & is.numeric(cor_alpha) & !is.null(data) & (alternative %in% c("less", "greater", "two_sided", "two_sided_signed"))){
    if(!is.numeric(n_repetitions)){
      stop('"n_repetitions" must be a positive integer. Check the parameters.')
    }
    significance_matrix <- matrix(0, nrow=nrow(cor_matrix), ncol=nrow(cor_matrix))
    significance_matrix[upper.tri(significance_matrix,diag=FALSE)] <- cor_permutationTest(data=as.matrix(data), n_repetitions=n_repetitions, alternative=alternative)[,3]
    adjacency_matrix <- (significance_matrix<cor_alpha)*adjacency_matrix
  }else if(!is.null(cor_alpha)&(sum(abs(adjacency_matrix))==0)){
    warning("Permutation test for the significance wasn't performed.")
  }else if(!is.null(cor_alpha)){
    warning("Permutation test for the significance wasn't performed. Check the parameters.")
  }
  if(sum(abs(adjacency_matrix))==0){
    g <- make_empty_graph(n=nrow(adjacency_matrix), directed=FALSE)
    g <- set_vertex_attr(g, "name", value=rownames(adjacency_matrix))
    return(g)
  }
  cg_rows <- unlist(sapply(1:(nrow(adjacency_matrix)-1), function(x) rownames(adjacency_matrix)[1:x]))
  cg_columns <- rep(rownames(adjacency_matrix)[2:nrow(adjacency_matrix)], 1:(nrow(adjacency_matrix)-1))
  edge_weights <- adjacency_matrix[upper.tri(adjacency_matrix, diag=FALSE)]
  edge_df <- data.frame(Node1=cg_rows[edge_weights!=0], Node2=cg_columns[edge_weights!=0], Cor=edge_weights[edge_weights!=0])
  node_df <- data.frame(IlmnID=rownames(adjacency_matrix))
  # normalization
  graph <- graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE)
  if(!is.null(normalization_fun)){
    tryCatch({
      graph <- set_edge_attr(graph, 'Cor', value=normalization_fun(edge_attr(graph, 'Cor')))
    }, error=function(e){
      warning('"normalization_fun" incorrectly specified. Normalization step is skipped.')
    })
  }
  #
  return(graph)
}

