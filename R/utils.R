# Validity of input
################################################################################
#' Checks if correlation matrix and data are correctly specified
#'
#' @description Function which checks whether provided correlation matrix and
#' data frame with data are correctly specified. Only one of two parameters must
#' be provided.
#' If `data` is given, function makes sure it is non-empty data frame with no 
#' NA values. For given `cor_matrix`, function makes sure it is non-empty 
#' non-negative symmetric matrix.
#' 
#' @param cor_matrix Correlation matrix.
#' @param data Data matrix. Correlation matrix is calculated with `cor(data)`.
#' 
#' @return List with correlation matrix and data. In case `data` is provided, 
#' `cor_matrix` is calculated from `data` and given `cor_matrix` is ignored.
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



#' Validity of CpG description file
#' 
#' @description Function which checks whether provided data frame `cg_meta` is
#' correctly specified.
#' Function makes sure that `cg_meta` is data frame with all `cg_meta_cols` 
#' columns and that all `necessary_cols` are names of elements in `cg_meta_cols`.
#' If `index_column` is given, function checks that column doesn't contain
#' duplicates in `cg_meta`.
#' 
#' @param cg_meta Data frame.
#' @param cg_meta_cols Named list where elements are columns of `cg_meta`.
#' @param neccessary_cols Vector where elements are names of elements in
#' `cg_meta_cols`.
#' @param index_column Name of index column in `cg_meta`.
#' 
#' @return Returns `cg_meta` data frame. If `index_column` was provided, it is
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
#' Type of gene region
#' 
#' @description Function which takes a name of gene region and checks it the
#' region is `"Promoter"`, `"Body"` or `"3'UTR"`. Known gene regions are:
#' `"TSS200"`, `"TSS1500"`, `"5'UTR"`, `"1stExon"`, `"Promoter"`.
#' 
#' @param gene_region Name of gene region.
#' 
#' @return Type of gene region: `"Promoter"`, `"Body"`, `"3'UTR"` or `""` for
#' an unknown region.
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



#' Types of multiple gene regions
#' 
#' @description Function which takes a vector of gene region names and return
#' their type without the repetition. Types of gene regions are `"Promoter"`, 
#' `"Body"` or `"3'UTR"`. Known gene regions are: `"TSS200"`, `"TSS1500"`, 
#' `"5'UTR"`, `"1stExon"`, `"Promoter"`.
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
# Help functions for plotting
################################################################################

#' Check validity of colors
#'
#' @description For a list of colors, function checks if given colors are well 
#' defined in R environment.
#' 
#' @param x Vector of colors.
#' 
#' @return Boolean vector.
#'
#' @details Uses `col2rgb` function with `alpha=TRUE`.
#'
#' @noRd
.is_color <- function(x){
  # we don't allow a data frame of colors
  if(is.data.frame(x)){
    return(FALSE)
  }
  
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X, alpha=TRUE)), 
             error = function(e) FALSE)}, USE.NAMES=FALSE)
  
}


#' Expand labels within a group
#' 
#' @description For given positions of labels new positions are returned so that
#' the distance between them is set to `dist`.
#' 
#' @param old_pos Given position of labels.
#' @param dist Distance between two consecutive labels.
#' 
#' @return Vector of new positions.
#' 
#' @details Starting from the middle labels are expanded so than the distance
#' between two consecutive labels is exactly equal to `dist`.
#'
#' @noRd
.new_position_within_group <- function(old_pos, dist)
{
  n <- length(old_pos)
  if(n==1)
  {
    return(old_pos)
  }
  new_pos <- rep(0,length(old_pos))
  
  # if there is an even number of labels
  if(n%%2==0){
    
    p <- old_pos[n/2+1]-old_pos[n/2] # distance between middle labels
    move <- (dist-p)/2
    new_pos[n/2+1] <- old_pos[n/2+1] + move
    new_pos[n/2] <- old_pos[n/2] - move
    if(n==2){
      return(new_pos)
    }
    for(i in 2:(n/2))
    {
      new_pos[n/2+i] <- new_pos[n/2+i-1] + dist
      new_pos[n/2-i+1] <- new_pos[n/2-i+2] - dist
    }
    
    # if there is an odd number of labels
  }else{''
    middle <- ceiling(n/2)
    new_pos[middle] <- old_pos[middle]
    for(i in 1:((n-1)/2))
    {
      new_pos[middle+i] <- new_pos[middle+i-1] + dist
      new_pos[middle-i] <- new_pos[middle-i+1] - dist
    }
  }
  return(new_pos)
  
}


#' Expand labels
#' 
#' @description  For given positions of labels new expanded positions are 
#' returned.
#' 
#' @param old_pos Given position of labels.
#' @param expand_factor Factor by which the labels are expanded. It represents
#' the percentage of the total axis length. Value of `x` represents `x%`.
#' 
#' @return Vector of new positions.
#' 
#' @details Firstly, all labels are separated into groups of cluttered labels.
#' Then for each group function `.new_positions_within_group` is called to
#' expand the labels. `expand_factor` represents the percentage of the total 
#' axis length (the maximum minus the minimum label position) by which two
#' consecutive labels should be separated.
#' 
#'
#' @noRd
.new_position_for_label <- function(old_pos, expand_factor)
{
  old_pos <- sort(old_pos)
  n <- length(old_pos)
  if(n==1){
    return(old_pos)
  }
  
  dist <- (old_pos[n]-old_pos[1])*expand_factor/100
  dd <- data.frame(old=old_pos,new=old_pos,group=0)
  dd[1,"group"] <- 1
  
  # separate labels into groups of "cluttered" labels
  count <- 1
  for(i in 2:n)
  {
    if(old_pos[i]-old_pos[i-1]>dist){
      count <- count+1
    }
    dd[i,"group"] <- count
  }
  
  # expand labels within each group
  n_groups <- length(unique(dd$group))
  first <- rep(0,n_groups)
  last <- rep(0,n_groups)
  for(i in 1:n_groups){
    new_values <- .new_position_within_group(dd[dd$group==i,"old"],dist)
    dd[dd$group==i,"new"] <- new_values
    first[i] <- min(new_values)
    last[i] <- max(new_values)
  }
  if(n_groups==1)
  {
    return(dd$new)
  }
  
  # if after expanding some groups clash, we expand them further
  if(n_groups%%2==0)
  {
    p <- first[n_groups/2+1] - last[n_groups/2]
    move <- (dist-p)/2 # if move<0 then groups are distant enough
    dd[dd$group==n_groups/2+1,"new"] <- dd[dd$group==n_groups/2+1,"new"] + max(0,move)
    dd[dd$group==n_groups/2,"new"] <- dd[dd$group==n_groups/2,"new"] - max(0,move)
    first[n_groups/2+1] <- first[n_groups/2+1] + max(0,move)
    last[n_groups/2] <- last[n_groups/2] - max(0,move)
    if(n_groups/2>=2){
      for(i in 2:(n_groups/2))
      {
        p <- first[n_groups/2+i]-last[n_groups/2+i-1]
        if(p<dist)
        {
          dd[dd$group==n_groups/2+i,"new"] <- dd[dd$group==n_groups/2+i,"new"] + dist - p
          first[n_groups/2+i] <- first[n_groups/2+i] + dist - p
          last[n_groups/2+i] <- last[n_groups/2+i] + dist - p
        }
        p <- first[n_groups/2-i+2]-last[n_groups/2-i+1]
        if(p<dist)
        {
          dd[dd$group==n_groups/2-i+1,"new"] <- dd[dd$group==n_groups/2-i+1,"new"] - dist + p
          first[n_groups/2-i+1] <- first[n_groups/2-i+1] - dist + p
          last[n_groups/2-i+1] <- last[n_groups/2-i+1] - dist + p
        }
      }
    }
  }else
  {
    middle <- ceiling(n_groups/2)
    for(i in 1:((n_groups-1)/2))
    {
      p <- first[middle+i]-last[middle+i-1]
      if(p<dist)
      {
        dd[dd$group==middle+i,"new"] <- dd[dd$group==middle+i,"new"] + dist - p
        first[middle+i] <- first[middle+i] + dist - p
        last[middle+i] <- last[middle+i] + dist - p
      }
      p <- first[middle-i+1]-last[middle-i]
      if(p<dist)
      {
        dd[dd$group==middle-i,"new"] <- dd[dd$group==middle-i,"new"] - dist + p
        first[middle-i] <- first[middle-i] - dist + p
        last[middle-i] <- last[middle-i] - dist + p
      }
    }
  }
  
  return(dd$new)
}

################################################################################


################################################################################
# Additional functions
################################################################################
#' Edge intersection for two `igraph` objects
#' 
#' @description Function which uses object from the `igraph` package. For two
#' graphs `g1` and `g2` and the corresponding subset of edges, `edgeL1` and
#' `edgeL2`, function returns the intersect of edges as subset of `g2` edges.
#' 
#' @param g1 First graph.
#' @param g2 Second graph.
#' @param edgeL1 List of `g1` edges to use. Default to all edges.
#' @param edgeL2 List of `g2` edges to use. Default to all edges.
#' 
#' @return Subset of `g2` edges which are found both in `edgeL1` and `edgeL2`
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



#' All pairwise combinations of elements
#'
#' @description For a list of CpG sites, function returns data frame with all 
#' pairwise combinations of them.
#' 
#' @param cg_list Vector of CpGs.
#' @param col1_name Name of the first column in the resulting data frame.
#' Defaults to `"Node1"`.
#' @param col2_name Name of the second column in the resulting data frame.
#' Defaults to `"Node2"`.
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



#' Equality of chromosomes
#' 
#' @description Checks if two chromosome names are the same taking into account 
#' different naming conventions, e.g. `2`, `"chr2"`, `"CHR2"`, `"Chr2"`.
#' 
#' @param chr_v1 First name.
#' @param chr_v2 Second name.
#' 
#' @return `TRUE` or `FALSE` depending whether two version of name describe the same
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



#' Creating sentence from words
#'
#' @description For a list of words, returns a single string of type 
#' `"word_1, word_2 and word_3"`. The conjuction `"and"` can be changed.
#' 
#' @param list_of_words Vector of words.
#' @param conjuction Word used between the second last and the last word.
#' Defaults to `"and"`.
#' 
#' @return String of merged words.
#'
#' @details Used by `.check_cgMeta_wCols` function.
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