# Functions which check validity of input
################################################################################
#' Checks if correlation matrix and data are correctly specified
#'
#' @description Function which checks whether provided correlation matrix and
#' data frame with data are correctly specified. Only one of two parameters must
#' be provided.
#' If "data" is given, function makes sure it is non-empty data frame with no 
#' NA values. For given "cor_matrix", function makes sure it is non-empty 
#' non-negative symmetric matrix.
#' 
#' @param cor_matrix Correlation matrix.
#' @param data Data matrix.
#' 
#' @return List with correlation matrix and data. In case "data" is provided, 
#' "cor_matrix" is calculated from "data" and given "cor_matrix" is ignored.
#' 
#' @noRd
.check_corM_dataDF <- function(cor_matrix, data){
  if(is.null(cor_matrix)&is.null(data)){
    stop(paste0('You must specify either "',deparse(substitute(cor_matrix)),'" or "',deparse(substitute(data)),'".'))
  }
  if(!is.null(data)){
    tryCatch(
      expr = {
        data <- data.frame(data)
      },
      error = function(e){
        stop(paste0('"',deparse(substitute(data)),'" must be convertible to the data frame.'))
      }
    )
    if(nrow(data)==0){
      stop(paste0('"',deparse(substitute(data)),'" specified, but empty. Set "',deparse(substitute(data)),'" to NULL to ignore it.'))
    }
    if(any(is.na(data))){
      stop(paste0('No NAs allowed in "',deparse(substitute(data)),'".'))
    }
    if(!all(is.numeric(as.matrix(data)))){
      stop(paste0('"',deparse(substitute(data)),'" must be numeric.'))
    }
    cor_matrix <- cor(data) # row and column names are taken from columns of data
  }else{ #cor_matrix is not null
    tryCatch(
      expr = {
        cor_matrix <- data.frame(cor_matrix)
      },
      error = function(e){
        stop(paste0('"',deparse(substitute(cor_matrix)),'" must be convertible to the data frame.'))
      }
    )
    if(dim(cor_matrix)[1]!=dim(cor_matrix)[2]){
      stop(paste0('"',deparse(substitute(cor_matrix)),'" must be a symmetric matrix.'))
    }else if(nrow(cor_matrix)==0){
      stop(paste0('"',deparse(substitute(cor_matrix)),'" is empty.'))
    }
    if(!isSymmetric(as.matrix(cor_matrix))){
      stop(paste0('"',deparse(substitute(cor_matrix)),'" must be a symmetric matrix.'))
    }
    if(any(is.na(cor_matrix))){
      stop(paste0('No NAs allowed in "',deparse(substitute(cor_matrix)),'".'))
    }
    if(!is.numeric(as.matrix(cor_matrix))){
      stop(paste0('"',deparse(substitute(cor_matrix)),'" must be numeric.'))
    }
  }
  #
  return(list(cor_matrix,data))
}



#' Checks if cg_meta is correctly specified
#' 
#' @description Function which checks whether provided data frame "cg_meta" is
#' correctly specified.
#' Function makes sure that "cg_meta" is data frame with all "cg_meta_cols" 
#' columns and that all "necessary_cols" are names of elements in "cg_meta_cols".
#' If "index_column" is given, function checks that column doesn't contain
#' duplicates in "cg_meta".
#'
#' 
#' @param cg_meta Data frame.
#' @param cg_meta_cols Named list where elements are columns of "cg_meta".
#' @param neccessary_cols Vector where elements are names of elements in
#' "cg_meta_cols".
#' @param index_column Name of index column in "cg_meta"
#' 
#' @return Returns "cg_meta" data frame. If "index_column" was provided, it is
#' used as row names in the resulting data frame.
#' 
#' @noRd
.check_cgMeta_wCols <- function(cg_meta, cg_meta_cols, necessary_cols, index_column){
  if(length(intersect(names(cg_meta_cols), necessary_cols))<length(necessary_cols)){
    message <- .words_listing(necessary_cols)
    if(grepl(" ", message)){
      stop(paste('Named list "',deparse(substitute(cg_meta_cols)),'" must contain elements with names',message,"."))
    }else{
      stop(paste('Named list "',deparse(substitute(cg_meta_cols)),'" must contain an element with name',message,"."))
    }
  }
  if(!inherits(cg_meta, "data.frame")){
    stop(paste0('"',deparse(substitute(cg_meta)),'" must be data frame.'))
  }
  if(!all(unlist(cg_meta_cols, use.names=FALSE) %in% colnames(cg_meta))){
    stop(paste0('"',deparse(substitute(cg_meta_cols)),'" are not columns of "',deparse(substitute(cg_meta)),'".'))
  }
  if(!is.null(index_column)){
    tryCatch(
      expr = {
        rownames(cg_meta) <- cg_meta[,cg_meta_cols[[index_column]]] # later calculations are easier with rownames
      },
      error = function(e){
        stop(paste0('Duplicate values in "',index_column,'" column of "',deparse(substitute(cg_meta)),'" data frame are not allowed.'))
      },
      warning = function(w){
        stop(paste0('Duplicate values in "',index_column,'" column of "',deparse(substitute(cg_meta)),'" data frame are not allowed.'))
      }
    )
  }
  return(cg_meta) #added rownames if index_column is given
}
################################################################################



################################################################################
# Functions which take gene region and check whether it is "Promoter", "Body" or "3'UTR" region
################################################################################
#' Checks the gene region
#' 
#' @description Function which takes a name of gene region and checks it the
#' region is "Promoter", "Body" or "3'UTR". Known gene regions are:
#' "TSS200", "TSS1500", "5'UTR", "1stExon", "Promoter".
#' 
#' @param gene_region Name of gene region.
#' 
#' @return Type of gene region: "Promoter", "Body", "3'UTR" or "" for unknown
#' region.
#' 
#' @noRd
.transform_one_gene_region <- function(gene_region){
  if(gene_region %in% c("TSS200", "TSS1500", "5'UTR", "1stExon", "Promoter")){
    return("Promoter")
  }else if(gene_region=="Body"){
    return(gene_region)
  }else if(gene_region=="3'UTR"){
    return(gene_region)
  }else{
    return("")
  }
}



#' Checks the vector of gene regions
#' 
#' @description Function which takes a vector of gene region names and return
#' their type without the repetition. Types of gene regions are "Promoter", 
#' "Body" or "3'UTR". Known gene regions are: "TSS200", "TSS1500", "5'UTR", 
#' "1stExon", "Promoter".
#' 
#' @param gene_region vector of gene region names.
#' 
#' @return Vector of gene region types. Only unique types are returned.
#' 
#' @noRd
.transform_gene_region <- function(gene_region){
  if(length(gene_region)>1){
    gene_region <- sapply(gene_region, .transform_one_gene_region)
    gene_region <- unique(gene_region)
    if(length(gene_region)==1){
      gene_region <- gene_region[1]
    }else{
      return(paste(gene_region,collapse=","))
    }
  }
  return(.transform_one_gene_region(gene_region))
}
################################################################################



################################################################################
# Additional functions
################################################################################
#' Intersects edges from two different "igraph" objects
#' 
#' @description Function which uses object from the "igraph" package. For two
#' graphs "g1" and "g2" and the corresponding subset of edges, "edgeL1" and
#' "edgeL2", function returns the intersect of edges as subset of "g2" edges.
#' 
#' @param g1 First graph.
#' @param g2 Second graph.
#' @param edgeL1 List of "g1" edges to use. Default to all edges.
#' @param edgeL2 List of "g2" edges to use. Default to all edges.
#' 
#' @return Subset of "g2" edges which are found both in "edgeL1" and "edgeL2"
#' 
#' @import igraph
#' 
#' @noRd
.intersect_edges <- function(g1, g2, edgeL1=E(g1), edgeL2=E(g2)){
  if(length(edgeL1)==0){
    return(E(g2)[0])
  }
  edge_list1 <- lapply(edgeL1, function(x) c(head_of(g1,x)$name, tail_of(g1,x)$name))
  edge_list2 <- sapply(edge_list1, function(x) ifelse((x[1]%in%V(g2)$name)&((x[2]%in%V(g2)$name)),get.edge.ids(g2, c(x[1], x[2])),0))
  edge_list2 <- intersect(E(g2)[edge_list2[edge_list2!=0]], edgeL2)
  return(E(g2)[edge_list2])
}



#' Creates data frame with all pairs of elements
#'
#' @description For a list of CpG sites, function returns data frame with all 
#' pairwise combinations of them.
#' 
#' @param cg_list Vector of CpGs.
#' @param col1_name Name of the first column in the resulting data frame.
#' Defaults to "Node1".
#' @param col2_name Name of the second column in the resulting data frame.
#' Defaults to "Node2".
#' 
#' @return Data frame with two columns. In every row, first column contains one
#' and the second column contains the other CpG in a pair. 
#' 
#' @noRd
.get_all_pairwiseCG <- function(cg_list, col1_name="Node1", col2_name="Node2"){
  cg_list <- unique(cg_list)
  df_temp <- data.frame(matrix(nrow=length(cg_list)*(length(cg_list)-1)/2,ncol=2))
  colnames(df_temp) <- c(col1_name, col2_name)
  if(length(cg_list)<=1){
    return(df_temp)
  }
  df_temp[,1] <- unlist(sapply(1:(length(cg_list)-1), function(x) cg_list[1:x]))
  df_temp[,2] <- unlist(rep(cg_list[2:length(cg_list)], 1:(length(cg_list)-1)))
  return(df_temp)
}



#' Compares chromosome names for different version of a name: e.g. 2, chr2, CHR2, Chr2
#'
#' @description Checks if two chromosome names are the same taking into account 
#' different naming conventions, e.g. "2", "chr2", "CHR2" and "Chr2".
#' 
#' @param chr_v1 First name.
#' @param chr_v2 Second name.
#' 
#' @return TRUE or FALSE depending whether two version of name describe the same
#' chromosome.
#' 
#' @details Function is slow when applied on large vector.
#' 
#' @noRd
.matchChr <- function(chr_v1, chr_v2){
  if(chr_v1==chr_v2){
    return(TRUE)
  }
  chr_v1 <- tolower(chr_v1)
  chr_v2 <- tolower(chr_v2)
  if(nchar(chr_v1)<4){
    chr_v1 <- paste0("chr", chr_v1)
  }
  if(nchar(chr_v2)<4){
    chr_v2 <- paste0("chr", chr_v2)
  }
  return(chr_v1==chr_v2)
}



#' Merges a list of words
#'
#' @description For a list of words, returns a single string of type 
#' "word_1, word_2 and word_3". The conjuction "and" can be changed.
#' 
#' @param list_of_words Vector of words.
#' @param conjuction Word used between the second last and the last word.
#' Defaults to "and".
#' 
#' @return String of merged words.
#'
#' @details Used by "".check_cgMeta_wCols" function.
#'
#' @noRd
.words_listing <- function(list_of_words, conjuction="and"){
  n_words <- length(list_of_words)
  if(n_words==0){
    return("")
  }else if(n_words==1){
    return(list_of_words[1])
  }else{
    return(paste(paste(list_of_words[1:(n_words-1)],collapse=", "),list_of_words[n_words],sep=paste0(" ",conjuction," ")))
  }
}
