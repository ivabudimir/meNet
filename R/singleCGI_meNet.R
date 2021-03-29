# check_matrices=FALSE when it is called inside other function since we already checked them
# clustering base on correlation
# return weighted (attribute Dist)
# much faster if instead of the whole corM you just give the part corresponding to the island (or just smaller)
#'@export
singleCGI_meNet <- function(cg_island, cor_matrix=NULL, data=NULL, link_method="twoLyr_clust", weighted=TRUE,
                            cor_normalization_fun=.max_normalization, dist_normalization_fun=.neg_max_normalization,
                            cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=NULL, cor_alpha=NULL, n_repetitions=1000, alternative="two_sided",
                            infomap_call="infomap", folder="./meNet/", file_basename="meNet_CGI_infomap", relaxation_rate=0.15,
                            cg_meta=cg_anno450k, cg_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", island_name="UCSC_CpG_Islands_Name", island_region="Relation_to_UCSC_CpG_Island"),
                            include_regions=c(), check_matrices=TRUE, delete_files=FALSE){
  #
  if(link_method=="twoLyr_clust"){
    error_message <- 'Infomap command not working. Change the "link_method" or check installation of "Infomap" and parameter "infomap_call". For more details, see https://www.mapequation.org/infomap/.'
    tryCatch({system(paste(infomap_call,"--help"))}, error=function(e){stop(error_message)}, warning=function(w){stop(error_message)})
  }
  #
  if(!(link_method%in%c("full", "clust", "twoLyr_clust"))){
    stop('"link_method" must have one of the values "full", "clust" or "twoLyr_clust".')
  }
  #
  if(check_matrices){
    cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "island_name", "island_region"), index_column="cg_id")
    if((weighted|(link_method=="twoLyr_clust"))&!("cg_coord"%in%names(cg_meta_cols))){
      stop('If distances among CpGs are used in network construction, element "cg_coord" must be provided in the named list "cg_meta_cols".')
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
  # cg_list
  cg_list <- CGinIsland(cg_island, cg_meta, cg_meta_cols, include_regions)
  if(length(cg_list)==0){
    stop(paste0('No CpGs found in CpG island "',cg_island ,'".'))
  }
  if(link_method!="full"){
    cg_list <- intersect(cg_list, colnames(cor_matrix))
    if(length(cg_list)==0){
      warning(paste0('Column names of "cor_matrix"/"data" don',"'t contain the CpGs found in CpG island ",'"',cg_island ,'".'))
      return(make_empty_graph(n=0, directed=FALSE))
    }
  }
  if(length(cg_list)==1){
    output_graph <- make_empty_graph(n=0, directed=FALSE)
    output_graph <- add_vertices(output_graph, length(cg_list), name=cg_list)
    return(output_graph)
  }
  # construction of the output graph
  node_df <- data.frame(IlmnID=cg_list)
  # the link method:
  if(link_method=="full"){ #fully connected graph
    edge_df <- data.frame(Node1=unlist(sapply(1:(length(cg_list)-1), function(x) cg_list[1:x])),
                          Node2=rep(cg_list[2:length(cg_list)], 1:(length(cg_list)-1)))
    if(weighted){
      distM <- as.matrix(dist(matrix(cg_meta[cg_list,cg_meta_cols$cg_coord], nrow=length(cg_list))))
      edge_df$Dist <- distM[upper.tri(distM, diag=FALSE)]
    }
    #
  }else if(link_method=="clust"){ #igraph infomap clustering
    graphCGI_cor <- meNet_cor(cor_matrix=cor_matrix, data=data, cg_ids=cg_list, cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold,
                              cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative)
    # remove negative edges
    graphCGI_cor <- delete_edges(graphCGI_cor, E(graphCGI_cor)[edge_attr(graphCGI_cor,"Cor")<0])
    if(length(E(graphCGI_cor))==0){
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      return(graph)
    }
    # clustering
    communities_list <- igraph::membership(igraph::cluster_infomap(graphCGI_cor, e.weights=E(graphCGI_cor)$Cor))
    cg_list_temp <- names(communities_list)
    cg_membership <- as.numeric(communities_list)
    communities <- as.numeric(names(table(cg_membership))[as.numeric(table(cg_membership))>1])
    if(length(communities)>0){
      edge_df <- dplyr::bind_rows(lapply(communities, function(x) .get_all_pairwiseCG(cg_list_temp[cg_membership==x],col1_name="Node1",col2_name="Node2")))
      if(weighted){
        distM <- as.matrix(dist(matrix(cg_meta[cg_list,cg_meta_cols$cg_coord], nrow=length(cg_list), dimnames=list(cg_list))))
        edge_df$Dist <- sapply(1:nrow(edge_df), function(x) distM[edge_df[x,1],edge_df[x,2]])
      }
    }else{
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      return(graph)
    }
  }else{ #2-layer infomap clustering (calls infomap on the command line)
    # cor layer
    graphCGI_cor <- meNet_cor(cor_matrix=cor_matrix, data=data, cg_ids=cg_list, cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold,
                              cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative)
    if(length(E(graphCGI_cor)[edge_attr(graphCGI_cor, "Cor")>0])==0){ #if no edges or only negative edges are left
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      return(graph)
    }
    # dist layer
    edgeDist_df <- data.frame(Node1=unlist(sapply(1:(length(cg_list)-1), function(x) cg_list[1:x])),
                              Node2=rep(cg_list[2:length(cg_list)], 1:(length(cg_list)-1)))
    distM <- as.matrix(dist(matrix(cg_meta[cg_list,cg_meta_cols$cg_coord], nrow=length(cg_list))))
    edgeDist_df$Dist <- distM[upper.tri(distM, diag=FALSE)]
    graphCGI_dist <- graph_from_data_frame(d=edgeDist_df, vertices=node_df, directed=FALSE)
    #
    communities_df <- meMultiplex_infomapCommunities(graphCGI_cor, graphCGI_dist, physical_nodes=FALSE, layer=1, folder=folder, file_basename=file_basename,
                                                     relaxation_rate=relaxation_rate, delete_files=delete_files, infomap_call=infomap_call,
                                                     cor_weighted=TRUE, supp_weighted=TRUE, cor_normalization_fun=.max_normalization, supp_normalization_fun=.neg_max_normalization,
                                                     inter_cor_supp=NULL, inter_supp_cor=NULL, infomap_seed=NULL)
    cg_list_temp <- communities_df[,1]
    cg_membership <- communities_df[,2]
    communities <- as.numeric(names(table(cg_membership))[as.numeric(table(cg_membership))>1])
    if(length(communities)>0){
      edge_df <- dplyr::bind_rows(lapply(communities, function(x) .get_all_pairwiseCG(cg_list_temp[cg_membership==x],col1_name="Node1",col2_name="Node2")))
      if(weighted){
        distM <- as.matrix(dist(matrix(cg_meta[cg_list,cg_meta_cols$cg_coord], nrow=length(cg_list), dimnames=list(cg_list))))
        edge_df$Dist <- sapply(1:nrow(edge_df), function(x) distM[edge_df[x,1],edge_df[x,2]])
      }
    }else{
      graph <- make_empty_graph(n=0, directed=FALSE)
      graph <- add_vertices(graph, length(cg_list), name=cg_list)
      return(graph)
    }
  }
  #
  return(graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE))
}
