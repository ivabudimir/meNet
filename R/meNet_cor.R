#' Builds a CpG network where edges are based on correlation
#' 
#' @description Different options are available to determine which correlations 
#' are kept in the network, such as a constant threshold for correlation, 
#' permutation test p-value or standard deviation of bootstrapping.
#' 
#' @param cor_matrix Correlation matrix of CpG sites.
#' @param data Data frame with CpGs in rows. Variables in columns are used to
#' calculate "cor_matrix".
#' @param cg_ids List of CpGs for which we reconstruct the network. If names of
#' CpGs are given as "cor_matrix" or "data" row names, "cg_ids" defines a subset
#' of CpGs which should be used in network. If omitted, all CpGs are used.
#' If CpG names are not given, "cg_ids" is a mandatory parameter which gives 
#' the names of CpGs.
#' @param cor_threshold
#' @param neg_cor_threshold
#' @param cor_stDev
#' cor_threshold/3
#' @param cor_alpha
#' @param n_repetitions
#' @param alternative
#' @param normalization_function
#' 
#' 
#' @return 
#' 
#' @details 
#' 
#' @import igraph
#' 
#' @export





# @details We don't remove isolated nodes
# cor_threshold=NA is the same as cor_threshold=0
# data samples in rows and variables in columns
#
# thress different methods for significance of correlation c("threshold", "permutation_test", "st_dev")
# if both cor_matrix and data are given, cor_matrix is ignored
#
# neg_cor_threshold can be used if we want to penalize stronger/less strong negative correlations
# e.g. if we use meNet_cor as one layer of meMultiplex then we may want to keep all negative correlations since they will influence other layers
#'@export
meNet_cor <- function(cor_matrix=NULL, data=NULL, cg_ids=NULL, cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=cor_threshold/3, cor_alpha=NULL, normalization_fun=NULL,
                      n_repetitions=100, alternative="two_sided"){
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
    significance_matrix[upper.tri(significance_matrix,diag=FALSE)] <- cor_meanSd(as.matrix(data), nrow(data), TRUE, n_repetitions)[,2]
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
    significance_matrix[upper.tri(significance_matrix,diag=FALSE)] <- cor_prTest(as.matrix(data), n_repetitions, alternative)
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

