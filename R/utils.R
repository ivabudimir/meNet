#' Match chromosome names
#
.matchChr <- function(chr_v1, chr_v2){
  # recognize different version of name, e.g. 2, chr2, CHR2, Chr2
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



# must be called with the variables
# checks if cor_matrix and data are correctly specified. In case of data!=null, cor_matrix is ignored
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


.words_listing <- function(list_of_words, quote='"'){
  n_words <- length(list_of_words)
  if(n_words==0){
    return("")
  }else if(n_words==1){
    return(paste0(quote,list_of_words[1],quote))
  }else{
    list_of_words <- paste0(quote,list_of_words,quote)
    return(paste(paste(list_of_words[1:(n_words-1)],collapse=", "),list_of_words[n_words],sep=" and "))
  }
}

# check if cg_meta is correctly specified
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

.clrs_for_region <- function(gene_region){
  if(gene_region=="Promoter"){
    return("cadetblue3")
  }else if(gene_region=="Body"){
    return("gold")
  }else if(gene_region=="3'UTR"){
    return("yellowgreen")
  }else{
    return("gray50")
  }
}


.max_normalization  <- function(x=numeric(0)){
  if(!is.numeric(x)){
    stop('"x" must be numeric.')
  }
  if(length(x)==0){
    stop('Empty vector.')
  }
  return(x/max(abs(x)))
}

# composition of this function is qual to .max_normalization
.neg_max_normalization  <- function(x=numeric(0), add_part=0.1){
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

#' Intersects edge lists of two graphs
.intersect_edges <- function(g1, g2, edgeL1=E(g1), edgeL2=E(g2)){
  if(length(edgeL1)==0){
    return(E(g2)[0])
  }
  edge_list1 <- lapply(edgeL1, function(x) c(head_of(g1,x)$name, tail_of(g1,x)$name))
  edge_list2 <- sapply(edge_list1, function(x) ifelse((x[1]%in%V(g2)$name)&((x[2]%in%V(g2)$name)),get.edge.ids(g2, c(x[1], x[2])),0))
  edge_list2 <- intersect(E(g2)[edge_list2[edge_list2!=0]], edgeL2)
  return(E(g2)[edge_list2])
}


.coord_sort <- function(cg_list, cg_meta=cg_anno450k, cg_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", cg_chr="CHR")){
  cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "cg_coord", "cg_chr"), index_column="cg_id")
  if(!all(cg_list %in% cg_meta[,cg_meta_cols$cg_id])){
    stop('Elements of "cg_list" not found in "cg_meta$cg_id".')
  }
  cg_meta <- cg_meta[cg_meta[,cg_meta_cols$cg_id]%in%cg_list,]
  cg_meta <- cg_meta[order(cg_meta[,cg_meta_cols$cg_chr], cg_meta[,cg_meta_cols$cg_coord]),]
  permutation <- sapply(cg_list, function(x) which(cg_meta[,cg_meta_cols$cg_id]==x), USE.NAMES=FALSE)
  return(permutation) #used by igraph::permute
}


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

