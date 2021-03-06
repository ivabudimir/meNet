#' CpG network: edges are based on gene affiliation
#' 
#' @description For every present gene, function "meNet_singleGene" is called 
#' and the resulting network is a union of the single-gene networks. 
#' (Optional) weights are chromosomal distances.
#' 
#' @param cg_list
#' 
#' @return 
#' 
#' @details 
#' 
#' 
#' if both cg_list and gene_list are given, we merge both lists;
# we leave also isolated CpGs
# if expand_cg_list=TRUE, we add all CpGs from genes
# gene list can be a list genes or transcripts (gene_id)
# we use the data frame constructed and saved (load with data in package)
#    * we, could also use Promoter=TSS200+TSS1500+5'UTR+1stExon, as in
#     "Validation of a DNA methylation microarray for 450,000 CpG sites 
#'     in the human genome" by J. Sandoval et al.
#' 
#' 
#' 
#' @import igraph
#' 
#' @export
meNet_gene <- function(cg_list=NULL, gene_list=NULL, weighted=TRUE, link_method="twoLyr_clust",
                       cor_matrix=NULL, data=NULL, cor_normalization_fun=max_normalization, dist_normalization_fun=neg_max_normalization,
                       cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=NULL, cor_alpha=NULL, n_repetitions=1000, alternative="two_sided",
                       infomap_call="infomap", folder="./meNet/", file_basename="meNet_CGI_infomap", relaxation_rate=0.15,
                       cgGene_meta=data("CpG_genes", package="meNet"), cgGene_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", gene_id="UCSC_RefGene_Name", gene_region="UCSC_RefGene_Group"),
                       gene_regions=c("Promoter", "Body", "3'UTR"), expand_cg_list=FALSE, normalization_fun=NULL, save_all_files=FALSE, delete_files=FALSE){
  #
  if(is.null(cg_list)&is.null(gene_list)){
    stop('You must specify either "cg_list" or "gene_list".')
  }
  #
  if(link_method%in%c("twoLyr_clust", "twoLyr_clust_wRegion")){
    error_message <- 'Infomap command not working. Change the "link_method" or check installation of "Infomap" and parameter "infomap_call". For more details, see https://www.mapequation.org/infomap/.'
    tryCatch({system(paste(infomap_call,"--help"))}, error=function(e){stop(error_message)}, warning=function(w){stop(error_message)})
  }
  #
  if(!(link_method%in%c("full", "clust", "twoLyr_clust", "twoLyr_clust_wRegion"))){
    stop('"link_method" must have one of the values "full", "clust" or "twoLyr_clust".')
  }
  #
  cgGene_meta <- .check_cgMeta_wCols(cgGene_meta, cgGene_meta_cols, necessary_cols=c("cg_id", "gene_id", "gene_region"), index_column=NULL)
  if((weighted|(link_method%in%c("twoLyr_clust","twoLyr_clust_wRegion")))&!("cg_coord"%in%names(cgGene_meta_cols))){
    stop('If distances among CpGs are used in network construction, element "cg_coord" must be provided in the named list "cgGene_meta_cols".')
  }
  #
  if(link_method!="full"){
    matrix_preprocessing <- .check_corM_dataDF(cor_matrix, data)
    cor_matrix <- matrix_preprocessing[[1]]
    data <- matrix_preprocessing[[2]]
    if(!all(rownames(cor_matrix)==colnames(cor_matrix))){
      stop(paste('"cor_matrix" row names',"don't match the column names."))
    }
  }else{
    cor_matrix <- NULL
    data <- NULL
  }
  # gene regions
  gene_regions <- intersect(gene_regions, c("Promoter", "TSS200", "TSS1500", "5'UTR", "1stExon", "Body", "3'UTR"))
  if(length(gene_regions)==0){
    stop('"gene_regions" incorrectly specified.')
  }
  if("Promoter" %in% gene_regions){
    gene_regions <- union(gene_regions, c("TSS200", "TSS1500", "5'UTR", "1stExon"))
    gene_regions <- gene_regions[gene_regions!="Promoter"]
  }
  #
  if(!is.null(cg_list)){
    cg_list <- intersect(cg_list, cgGene_meta[,cgGene_meta_cols$cg_id])
    if(!is.null(cor_matrix)){
      cg_list <- intersect(cg_list, colnames(cor_matrix))
    }
    gene_to_add <- unique(cgGene_meta[(cgGene_meta[,cgGene_meta_cols$cg_id]%in%cg_list)&(cgGene_meta[,cgGene_meta_cols$gene_region]%in%gene_regions),cgGene_meta_cols$gene_id])
    if(!is.null(gene_list)){
      gene_list <- intersect(gene_list, cgGene_meta[,cgGene_meta_cols$gene_id])
    }else{
      gene_list <- c()
    }
  }else{
    cg_list <- c()
    gene_list <- unique(intersect(gene_list, cgGene_meta[,cgGene_meta_cols$gene_id]))
    gene_to_add <- c()
  }
  #
  if((length(cg_list)==0)&(length(gene_list)==0)){
    return(make_empty_graph(n=0,directed=FALSE))
  }
  # create nodes and edge data frame
  if(weighted){
    edge_df <- data.frame(matrix(nrow=0, ncol=3))
    colnames(edge_df) <- c("Node1", "Node2", "Dist")
  }else{
    edge_df <- data.frame(matrix(nrow=0, ncol=2))
    colnames(edge_df) <- c("Node1", "Node2")
  }
  node_df <- data.frame(matrix(nrow=0,ncol=2))
  colnames(node_df) <- c("IlmnID", "Region")
  #
  # loop over all genes: add CGs
  for(i in 1:length(gene_list)){
    cg_gene_i <- unique(cgGene_meta[(cgGene_meta[,cgGene_meta_cols$gene_id]==gene_list[i])&(cgGene_meta[,cgGene_meta_cols$gene_region]%in%gene_regions),cgGene_meta_cols$cg_id])
    if(!is.null(cor_matrix)){
      cg_gene_i <- intersect(cg_gene_i, colnames(cor_matrix))
      corM_i <- cor_matrix[cg_gene_i,cg_gene_i]
      if(!is.null(data)){
        dataM_i <- data[,cg_gene_i]
      }else{
        dataM_i <- NULL
      }
    }else{
      corM_i <- NULL
    }
    if(length(cg_gene_i)==0){
      next
    }else if(length(cg_gene_i)==1){
      cg_list <- unique(c(cg_list, cg_gene_i))
      next
    }
    if(save_all_files){
      file_basename <- cgi_list[i]
      delete_files <- FALSE
    }else{
      file_basename <- file_basename
    }
    graph_gene <- meNet_singleGene(gene_list[i], cor_matrix=corM_i, data=dataM_i, link_method=link_method, weighted=weighted,
                                   cor_normalization_fun=cor_normalization_fun, dist_normalization_fun=dist_normalization_fun,
                                   cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold, cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative,
                                   infomap_call=infomap_call, folder=folder, file_basename=file_basename, relaxation_rate=relaxation_rate,
                                   cgGene_meta=cgGene_meta, cgGene_meta_cols=cgGene_meta_cols, gene_regions=gene_regions, check_matrices=FALSE, delete_files=delete_files)
    edge_df_i <- igraph::as_data_frame(graph_gene)
    colnames(edge_df_i)[1:2] <- c("Node1", "Node2")
    edge_df <- rbind(edge_df, edge_df_i)
    cg_list <- unique(c(cg_list, cg_gene_i))
    node_df_i <- igraph::as_data_frame(graph_gene, what="vertices")
    colnames(node_df_i) <- c("IlmnID", "Region")
    node_df <- rbind(node_df, node_df_i)
  }
  #
  # we separate them so that we can control the addition of CpGs to the list of nodes (when gene is given, all CpGs are added)
  gene_to_add <- setdiff(gene_to_add, gene_list)
  for(i in 1:length(gene_to_add)){
    cg_gene_i <- unique(cgGene_meta[(cgGene_meta[,cgGene_meta_cols$gene_id]==gene_to_add[i])&(cgGene_meta[,cgGene_meta_cols$gene_region]%in%gene_regions),cgGene_meta_cols$cg_id])
    if(!is.null(cor_matrix)){
      cg_gene_i <- intersect(cg_gene_i, colnames(cor_matrix))
      corM_i <- cor_matrix[cg_gene_i,cg_gene_i]
      if(!is.null(data)){
        dataM_i <- data[,cg_gene_i]
      }else{
        dataM_i <- NULL
      }
    }else{
      corM_i <- NULL
    }
    if(length(cg_gene_i)<=1){
      next
    }
    if(save_all_files){
      file_basename <- cgi_list[i]
      delete_files <- FALSE
    }else{
      file_basename <- file_basename
    }
    graph_gene <- meNet_singleGene(gene_to_add[i], cor_matrix=corM_i, data=dataM_i, link_method=link_method, weighted=weighted,
                                   cor_normalization_fun=cor_normalization_fun, dist_normalization_fun=dist_normalization_fun,
                                   cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold, cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative,
                                   infomap_call=infomap_call, folder=folder, file_basename=file_basename, relaxation_rate=relaxation_rate,
                                   cgGene_meta=cgGene_meta, cgGene_meta_cols=cgGene_meta_cols, gene_regions=gene_regions, check_matrices=FALSE, delete_files=delete_files)
    if(expand_cg_list){
      cg_list <- unique(c(cg_list, cg_gene_i))
    }else{
      nodes_to_delete <- setdiff(V(graph_gene)$name,cg_list)
      graph_gene <- delete_vertices(graph_gene, V(graph_gene)[nodes_to_delete])
    }
    edge_df_i <- igraph::as_data_frame(graph_gene)
    colnames(edge_df_i)[1:2] <- c("Node1", "Node2")
    edge_df <- rbind(edge_df, edge_df_i)
    node_df_i <- igraph::as_data_frame(graph_gene, what="vertices")
    colnames(node_df_i) <- c("IlmnID", "Region")
    node_df <- rbind(node_df, node_df_i)
  }
  # final construction
  if(nrow(edge_df)==0){
    graph <- make_empty_graph(n=0, directed=FALSE)
    graph <- add_vertices(graph, length(cg_list), name=cg_list)
    return(graph)
  }
  edge_df[,1:2] <- t(apply(edge_df[,1:2], 1, sort)) #in case nodes are in different order
  edge_df <- edge_df[!duplicated(edge_df),] # just in case some CpGs are shared among multiple genes
  #
  node_df <- node_df[!duplicated(node_df),] # if all are duplicated
  node_df[duplicated(node_df$IlmnID)|duplicated(node_df$IlmnID, fromLast=TRUE),"Region"] <- "mix"
  node_df <- node_df[!duplicated(node_df),]
  cg_to_add <- setdiff(cg_list, node_df$IlmnID)
  if(length(cg_to_add)>0){
    node_df_temp <- data.frame(matrix("",nrow=length(cg_to_add),ncol=2))
    colnames(node_df_temp) <- c("IlmnID", "Region")
    node_df_temp$IlmnID <- cg_to_add
    node_df <- rbind(node_df, node_df_temp)
  }
  # normalization
  graph <- graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE)
  if(!is.null(normalization_fun)){
    tryCatch({
      graph <- set_edge_attr(graph, 'Dist', value=normalization_fun(edge_attr(graph, 'Dist')))
    }, error=function(e){
      warning('"normalization_fun" incorrectly specified. Normalization step is skipped.')
    })
  }
  #
  return(graph)
}

