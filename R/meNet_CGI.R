#' Builds a CpG network where edges are based on CpG island affiliation
#'
# @details if both cg_list and cgi_list are given, we merge both lists;
# we leave also isolated CpGs
# if expand_cg_list=TRUE, we add all CpGs from islands
# weights are distances
# edge_method = "full", "clust", "twoLyr_clust" (if clustering with distance check infomap installation)
# weights = "dist" or unweighted
#'@export
meNet_CGI <- function(cg_list=NULL, cgi_list=NULL, weighted=TRUE, link_method="twoLyr_clust",
                      cor_matrix=NULL, data=NULL, cor_normalization_fun=.max_normalization, dist_normalization_fun=.neg_max_normalization,
                      cor_threshold=0.2, neg_cor_threshold=NULL, cor_stDev=NULL, cor_alpha=NULL, n_repetitions=1000, alternative="two_sided",
                      infomap_call="infomap", folder="./meNet/", file_basename="meNet_CGI_infomap", relaxation_rate=0.15,
                      cg_meta=cg_anno450k, cg_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", island_name="UCSC_CpG_Islands_Name", island_region="Relation_to_UCSC_CpG_Island"),
                      include_regions=character(0), expand_cg_list=FALSE, normalization_fun=NULL, save_all_files=FALSE, delete_files=FALSE){
  #
  if(is.null(cg_list)&is.null(cgi_list)){
    stop('You must specify either "cg_list" or "cgi_list".')
  }
  #
  if(!(link_method%in%c("full", "clust", "twoLyr_clust"))){
    stop('"link_method" must have one of the values "full", "clust" or "twoLyr_clust".')
  }
  if(link_method=="twoLyr_clust"){
    error_message <- 'Infomap command not working. Change the "link_method" or check installation of "Infomap" and parameter "infomap_call". For more details, see https://www.mapequation.org/infomap/.'
    tryCatch({system(paste(infomap_call,"--help"))}, error=function(e){stop(error_message)}, warning=function(w){stop(error_message)})
  }
  #
  cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "island_name", "island_region"), index_column="cg_id")
  if((weighted|(link_method=="twoLyr_clust"))&!("cg_coord"%in%names(cg_meta_cols))){
    stop('If distances among CpGs are used in network construction, element "cg_coord" must be provided in the named list "cg_meta_cols".')
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
  #
  if(!is.null(cg_list)){
    cg_list <- intersect(cg_list, cg_meta[,cg_meta_cols$cg_id])
    if(!is.null(cor_matrix)){
      cg_list <- intersect(cg_list, colnames(cor_matrix))
    }
    island_regions <- c("Island")
    island_regions <- unique(c(island_regions, intersect(include_regions, c('N_Shore','S_Shelf','S_Shore','N_Shelf'))))
    cgi_to_add <- unique(cg_meta[(cg_meta[,cg_meta_cols$cg_id]%in%cg_list)&(cg_meta[,cg_meta_cols$island_region]%in%island_regions),cg_meta_cols$island_name])
    # there won't be "" (open sea CpGs) because we specified island_region
    if(!is.null(cgi_list)){
      cgi_list <- unique(intersect(cgi_list, cg_meta[,cg_meta_cols$island_name]))
    }else{
      cgi_list <- c()
    }
  }else{
    cg_list <- c()
    cgi_list <- unique(intersect(cgi_list, cg_meta[,cg_meta_cols$island_name]))
    cgi_to_add <- c()
  }
  if((length(cg_list)==0)&(length(cgi_list)==0)){
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
  # loop over all islands: add CGs
  for(i in 1:length(cgi_list)){
    cg_island_i <- CGinIsland(cgi_list[i], cg_meta, cg_meta_cols, include_regions)
    if(!is.null(cor_matrix)){
      cg_island_i <- intersect(cg_island_i, colnames(cor_matrix))
      corM_i <- cor_matrix[cg_island_i,cg_island_i]
      if(!is.null(data)){
        dataM_i <- data[,cg_island_i]
      }else{
        dataM_i <- NULL
      }
    }else{
      corM_i <- NULL
    }
    if(length(cg_island_i)==0){
      next
    }else if(length(cg_island_i)==1){
      cg_list <- unique(c(cg_list, cg_island_i))
      next
    }
    if(save_all_files){
      file_basename <- cgi_list[i]
      delete_files <- FALSE
    }else{
      file_basename <- file_basename
    }
    graphCGI <- singleCGI_meNet(cgi_list[i], cor_matrix=corM_i, data=dataM_i, link_method=link_method, weighted=weighted,
                                cor_normalization_fun=cor_normalization_fun, dist_normalization_fun=dist_normalization_fun,
                                cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold, cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative,
                                infomap_call=infomap_call, folder=folder, file_basename=file_basename, relaxation_rate=relaxation_rate,
                                cg_meta=cg_meta, cg_meta_cols=cg_meta_cols, include_regions=include_regions, check_matrices=FALSE, delete_files=delete_files)
    edge_df_i <- igraph::as_data_frame(graphCGI)
    colnames(edge_df_i)[1:2] <- c("Node1", "Node2")
    edge_df <- rbind(edge_df, edge_df_i)
    cg_list <- unique(c(cg_list, cg_island_i))
  }
  # we separate them so that we can control the addition of CpGs to the list of nodes
  cgi_to_add <- setdiff(cgi_to_add, cgi_list)
  for(i in 1:length(cgi_to_add)){
    cg_island_i <- CGinIsland(cgi_to_add[i], cg_meta, cg_meta_cols, include_regions)
    if(!is.null(cor_matrix)){
      cg_island_i <- intersect(cg_island_i, colnames(cor_matrix))
      corM_i <- cor_matrix[cg_island_i,cg_island_i]
      if(!is.null(data)){
        dataM_i <- data[,cg_island_i]
      }else{
        dataM_i <- NULL
      }
    }else{
      corM_i <- NULL
    }
    if(length(cg_island_i)<=1){
      next
    }
    if(save_all_files){
      file_basename <- cgi_list[i]
      delete_files <- FALSE
    }else{
      file_basename <- file_basename
    }
    graphCGI <- singleCGI_meNet(cgi_to_add[i], cor_matrix=corM_i, data=dataM_i, link_method=link_method, weighted=weighted,
                                cor_normalization_fun=cor_normalization_fun, dist_normalization_fun=dist_normalization_fun,
                                cor_threshold=cor_threshold, neg_cor_threshold=neg_cor_threshold, cor_stDev=cor_stDev, cor_alpha=cor_alpha, n_repetitions=n_repetitions, alternative=alternative,
                                infomap_call=infomap_call, folder=folder, file_basename=file_basename, relaxation_rate=relaxation_rate,
                                cg_meta=cg_meta, cg_meta_cols=cg_meta_cols, include_regions=include_regions, check_matrices=FALSE, delete_files=delete_files)
    #
    if(expand_cg_list){
      cg_list <- unique(c(cg_list, cg_island_i))
    }else{
      nodes_to_delete <- setdiff(V(graphCGI)$name,cg_list)
      graphCGI <- delete_vertices(graphCGI, V(graphCGI)[nodes_to_delete])
    }
    edge_df_i <- igraph::as_data_frame(graphCGI)
    colnames(edge_df_i)[1:2] <- c("Node1", "Node2")
    edge_df <- rbind(edge_df, edge_df_i)
  }
  # final construction
  node_df <- data.frame(IlmnID=cg_list)
  if(nrow(edge_df)==0){
    graph <- make_empty_graph(n=0, directed=FALSE)
    graph <- add_vertices(graph, length(cg_list), name=cg_list)
    return(graph)
  }
  edge_df[,1:2] <- t(apply(edge_df[,1:2], 1, sort)) #in case nodes are in different order
  edge_df <- edge_df[!duplicated(edge_df),] # just in case some CpGs are shared among multiple islands (if we include shores and shelves)
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
