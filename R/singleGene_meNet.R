#' 
#' @description Builds a CpG network for a single gene with (optional) weights 
#' representing chromosomal distances between CpGs. 
#' Edge significance is determined similarly as in singleCGI_meNet with an 
#' additional clustering method “twoLyr_clust_wRegion” which additionally 
#' incorporates the gene-region-information (e.g. “Promoter”, “Body”) 
#' in the distance layer prior to the Infomap multiplex clustering. 
#' 
#' @param cor_layer
#' 
#' @return 
#' 
#' @details 
#' 
#' @import igraph
#' 
#' @export

# return also attribute "Region"
#'@export
singleGene_meNet <- function(gene_name, cor_matrix=NULL, data=NULL, link_method="twoLyr_clust", weighted=TRUE,
                             cor_normalization_fun=max_normalization, dist_normalization_fun=neg_max_normalization,
                             cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=NULL, cor_alpha=NULL, n_repetitions=1000, alternative="two_sided",
                             infomap_call="infomap", folder="./meNet/", file_basename="meNet_gene_infomap", relaxation_rate=0.15,
                             cgGene_meta=data("CpG_genes", package="meNet"), cgGene_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", gene_id="UCSC_RefGene_Name", gene_region="UCSC_RefGene_Group"),
                             gene_regions=c("Promoter", "Body", "3'UTR"), check_matrices=TRUE, delete_files=FALSE){
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
  if(check_matrices){
    cgGene_meta <- .check_cgMeta_wCols(cgGene_meta, cgGene_meta_cols, necessary_cols=c("cg_id", "gene_id", "gene_region"), index_column=NULL)
    if((weighted|(link_method%in%c("twoLyr_clust","twoLyr_clust_wRegion")))&!("cg_coord"%in%names(cgGene_meta_cols))){
      stop('If distances among CpGs are used in network construction, element "cg_coord" must be provided in the named list "cgGene_meta_cols".')
    }
  }
  if(link_method!="full"&check_matrices){
    matrix_preprocessing <- .check_corM_dataDF(cor_matrix, data)
    cor_matrix <- matrix_preprocessing[[1]]
    data <- matrix_preprocessing[[2]]
    if(!all(rownames(cor_matrix)==colnames(cor_matrix))){
      stop(paste('"cor_matrix" row names',"don't match the column names."))
    }
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
  # cg_list
  cgGene_meta <- cgGene_meta[cgGene_meta[,cgGene_meta_cols$gene_id]==gene_name,]
  if(nrow(cgGene_meta)==0){
    stop(paste0('Gene "',gene_name,'" is not found in "cgGene_meta".'))
  }
  if(weighted|(link_method%in%c("twoLyr_clust","twoLyr_clust_wRegion"))){
    dist_df <- cgGene_meta[,c(cgGene_meta_cols$cg_id,cgGene_meta_cols$cg_coord)]
    dist_df <- dist_df[!duplicated(dist_df),]
    rownames(dist_df) <- dist_df[,cgGene_meta_cols$cg_id]
  }
  row_indices <- sapply(cgGene_meta[,cgGene_meta_cols$gene_region], function(x) x %in% gene_regions, USE.NAMES=FALSE)
  cg_list <- unique(cgGene_meta[row_indices,cgGene_meta_cols$cg_id]) #there may be more regions for the same gene
  if(length(cg_list)==0){
    stop(paste0('No CpGs found in gene "', gene_name,'".'))
  }
  if(link_method!="full"){
    cg_list <- intersect(cg_list, colnames(cor_matrix))
    if(length(cg_list)==0){
      warning(paste0('Column names of "cor_matrix"/"data" don',"'t contain the CpGs found in gene ",'"',gene_name ,'".'))
      return(make_empty_graph(n=0, directed=FALSE))
    }
  }
  if(length(cg_list)==1){
    output_graph <- make_empty_graph(n=0, directed=FALSE)
    output_graph <- add_vertices(output_graph, length(cg_list), name=cg_list)
    graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
    return(output_graph)
  }
  # construction of the output graph
  node_df <- data.frame(IlmnID=cg_list)
  #
  # the link method:
  if(link_method=="full"){ #fully connected graph
    edge_df <- data.frame(Node1=unlist(sapply(1:(length(cg_list)-1), function(x) cg_list[1:x])),
                          Node2=rep(cg_list[2:length(cg_list)], 1:(length(cg_list)-1)))
    if(weighted){
      distM <- as.matrix(dist(matrix(dist_df[cg_list,cgGene_meta_cols$cg_coord], nrow=length(cg_list))))
      edge_df$Dist <- distM[upper.tri(distM, diag=FALSE)]
    }
    #
  }else if(link_method=="clust"){ #igraph infomap clustering
    graphGene_cor <- meNet_cor(cor_matrix=cor_matrix, data=data, cg_ids=cg_list, cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold,
                               cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative)
    # remove negative edges
    graphGene_cor <- delete_edges(graphGene_cor, E(graphGene_cor)[edge_attr(graphGene_cor,"Cor")<0])
    if(length(E(graphGene_cor))==0){
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
      return(graph)
    }
    # clustering
    communities_list <- igraph::membership(igraph::cluster_infomap(graphGene_cor, e.weights=E(graphGene_cor)$Cor))
    cg_list_temp <- names(communities_list)
    cg_membership <- as.numeric(communities_list)
    communities <- as.numeric(names(table(cg_membership))[as.numeric(table(cg_membership))>1])
    if(length(communities)>0){
      edge_df <- dplyr::bind_rows(lapply(communities, function(x) .get_all_pairwiseCG(cg_list_temp[cg_membership==x],col1_name="Node1",col2_name="Node2")))
      if(weighted){
        distM <- as.matrix(dist(matrix(dist_df[cg_list,cgGene_meta_cols$cg_coord], nrow=length(cg_list), dimnames=list(cg_list))))
        edge_df$Dist <- sapply(1:nrow(edge_df), function(x) distM[edge_df[x,1],edge_df[x,2]])
      }
    }else{
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
      return(graph)
    }
  }else{ #2-layer infomap clustering (calls infomap on the command line)
    # cor layer
    graphGene_cor <- meNet_cor(cor_matrix=cor_matrix, data=data, cg_ids=cg_list, cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold,
                               cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative)
    if(length(E(graphGene_cor)[edge_attr(graphGene_cor, "Cor")>0])==0){ #if no edges or only negative edges are left
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
      return(graph)
    }
    if(link_method=="twoLyr_clust"){ #only difference between "twoLyr_clust" and "twoLyr_clust_wRegion" is in the distance layer
      # dist layer
      edgeDist_df <- data.frame(Node1=unlist(sapply(1:(length(cg_list)-1), function(x) cg_list[1:x])),
                                Node2=rep(cg_list[2:length(cg_list)], 1:(length(cg_list)-1)))
      distM <- as.matrix(dist(matrix(dist_df[cg_list,cgGene_meta_cols$cg_coord], nrow=length(cg_list))))
      edgeDist_df$Dist <- distM[upper.tri(distM, diag=FALSE)]
      #
    }else{
      edgeDist_df <- data.frame(matrix(nrow=0,ncol=3))
      colnames(edgeDist_df) <- c("Node1", "Node2", "Dist")
      cgGene_meta$region_temp <- sapply(cgGene_meta[,cgGene_meta_cols$gene_region], function(x) .transform_gene_region(x))
      for(i in 1:3){
        region <- c("Promoter", "Body", "3'UTR")[i] # there may be duplicates
        cg_list_temp <- intersect(cgGene_meta[cgGene_meta$region_temp==region,cgGene_meta_cols$cg_id], cg_list)
        if(length(cg_list_temp)<=1){
          next
        }
        edgeDist_df_temp <- data.frame(Node1=unlist(sapply(1:(length(cg_list_temp)-1), function(x) cg_list_temp[1:x])),
                                       Node2=rep(cg_list_temp[2:length(cg_list_temp)], 1:(length(cg_list_temp)-1)))
        distM <- as.matrix(dist(matrix(dist_df[cg_list_temp,cgGene_meta_cols$cg_coord], nrow=length(cg_list_temp))))
        edgeDist_df_temp$Dist <- distM[upper.tri(distM, diag=FALSE)]
        edgeDist_df <- rbind(edgeDist_df, edgeDist_df_temp)
      }
      edgeDist_df[,1:2] <- t(apply(edgeDist_df[,1:2], 1, sort)) #in case nodes are in different order
      edgeDist_df <- edgeDist_df[!duplicated(edgeDist_df),] # just in case some CpGs are shared among multiple regions
      # not sure what happens it there are no edges (calls infomap)
    }
    graphGene_dist <- graph_from_data_frame(d=edgeDist_df, vertices=node_df, directed=FALSE)
    communities_df <- meMultiplex_communities(graphGene_cor, graphGene_dist, physical_nodes=FALSE, layer=1, folder=folder, file_basename=file_basename,
                                              relaxation_rate=relaxation_rate, delete_files=delete_files, infomap_call=infomap_call,
                                              cor_weighted=TRUE, supp_weighted=TRUE, cor_normalization_fun=.max_normalization, supp_normalization_fun=.neg_max_normalization,
                                              inter_cor_supp=NULL, inter_supp_cor=NULL, infomap_seed=NULL)
    #
    cg_list_temp <- communities_df[,1]
    cg_membership <- communities_df[,2]
    communities <- as.numeric(names(table(cg_membership))[as.numeric(table(cg_membership))>1])
    if(length(communities)>0){
      edge_df <- dplyr::bind_rows(lapply(communities, function(x) .get_all_pairwiseCG(cg_list_temp[cg_membership==x],col1_name="Node1",col2_name="Node2")))
      if(weighted){
        distM <- as.matrix(dist(matrix(dist_df[cg_list,cgGene_meta_cols$cg_coord], nrow=length(cg_list), dimnames=list(cg_list))))
        edge_df$Dist <- sapply(1:nrow(edge_df), function(x) distM[edge_df[x,1],edge_df[x,2]])
      }
    }else{
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
      return(graph)
    }
  }
  #
  graph <- graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE)
  graph <- set_vertex_attr(graph, 'Region', index=V(graph), sapply(V(graph)$name, function(x) .transform_gene_region(cgGene_meta[cgGene_meta[,cgGene_meta_cols$cg_id]==x,cgGene_meta_cols$gene_region]), USE.NAMES=FALSE))
  return(graph)
}
