#' Multiplex: `meNet` file from two `igraph` layers
#'
#' @description Constructs a multiplex network from two `igraph` layers and 
#' writes the multiplex structure to file. Layers should share some, but
#' not necessary all, nodes. Function interprets the first layer as correlation 
#' layer. The second layer, if weighted, should have only non.negative weights and 
#' can be interpreted as distance layer. Multiplex structure can be written as
#' `Infomap` of `multinet` object.
#' 
#' @param cor_layer Correlation layer of the multiplex as `igraph` object.
#' @param supplementary_layer Supplementary layer of the multiples as `igraph`
#' object.
#' @param output_file Name of the output file to which multiplex structure
#' is written.
#' @param cor_weighted Whether the correlation layer is weighted. If `FALSE`,
#' all weights are set to `1`. Defaults to `TRUE`.
#' @param supp_weighted Whether the supplementary layer is weighted. If `FALSE`,
#' all weights are set to `1`. Defaults to `TRUE`.
#' @param cor_normalization_fun Normalization function applied on the weights
#' of the correlation layer, if the layer is weighted. If `NULL`, no normalization
#' is applied. Defaults to `max_normalization`.
#' @param supp_normalization_fun Normalization function applied on the weights
#' of the supplementary layer, if the layer is weighted. If `NULL`, no normalization
#' is applied. Defaults to `neg_max_normalization`.
#' @param output_type Structure of the output file. Has to be one of the values
#' "infomap"`` or `"multinet"`. Default value is `"infomap"`.
#' @param inter_cor_supp Weight of the inter-layer edges from the correlation
#' layer to the supplementary layer. Defaults to `NULL`.
#' @param inter_supp_cor Weight of the inter-layer edges from the supplementary
#' layer to the correlation layer.
#' By default, the values is equal to `inter_cor_supp`.
#' 
#' @return A vector of removed nodes. Invisibly, structure of the created
#' multiplex is saved to the `output_file`.
#'
#' @details 
#' Layers could previously be constructed using some of the functions 
#' `meNet_cor`, `meNet_CGI` or `meNet_gene`.
#' If layers are weighted, function searches for the edge attributes with
#' names`"Cor"` and `"Dist"`. If the attribute doesn't exist, it takes the first
#' edge attribute instead. The function handles weight of 
#' the first layer as correlations and weights of the second layer as distances. 
#'
#' Function checks if all edge weights of `supplementary_layer` are positive.
#' For negative weights of `correlation_layer`, corresponding edge is removed
#' from both layers. In this step more power is given to significantly negative 
#' correlations in community breakage.
#' After the removal of negative edges, nodes which are isolated in both layers 
#' are removed from the multiplex.
#' 
#' Additional parameters can be set, such as the weight of inter-layer edges.
#' `inter_cor_supp` defines the weight of all inter-layer edges from the 
#' correlation layer to the supplementary layer, while `inter_supp_cor` 
#' defines the weight of all inter-layer edges from the supplementary
#' layer to the correlation layer. By default, both values are `NULL`.
#' If `NULL`, for `Infomap` output inter-layer edges are not added to file while
#' for `multinet` output the type of network is defined as `multiplex` instead
#' of `multilayer`.
#' 
#' If a layer is weighted, normalization function can be used on the
#' edge weights. For correlation layer the default is maximum normalization using 
#' `max_normalization` function. For supplementary layer, the default is
#' negative maximum normalization with `neg_max_normalization`. The difference
#' in normalization function exists because weights of correlation layer can be 
#' interpreted as "strengths" while the weights of supplementary layer can be
#' interpreted as "cost". That is, large correlation has opposite impact
#' than the large distance for a pair of nodes.
#' 
#' Resulting multiplex structure is written to a file either in the `Infomap` 
#' \insertCite{infomap}{meNet} format or in the `multinet` 
#' \insertCite{multinet}{meNet} format, depending on the `output_type`.
#' `Infomap`-style file can be used by `Infomap` algorithm online or from the 
#' command line upon installation.
#' `multinet`-style file can be converted to `multinet` object using
#' `meMultiplex_to_multinet` function.
#' 
#' @references
#'       \insertAllCited{}
#' 
#' @import igraph
#' 
#' @export
meMultiplex <- function(cor_layer, supplementary_layer, output_file, cor_weighted=TRUE, supp_weighted=TRUE,
                        cor_normalization_fun=meNet::max_normalization, supp_normalization_fun=meNet::neg_max_normalization, 
                        output_type="infomap", inter_cor_supp=NULL, inter_supp_cor=inter_cor_supp){
  if(!(output_type) %in% c("infomap", "multinet")){
    stop('"output_type" must be either "infomap" or "multinet".')
  }
  #
  if((!is.null(inter_cor_supp))&(!is.numeric(inter_cor_supp))){
    stop('"inter_cor_supp" must be numeric.')
  }
  if((!is.null(inter_supp_cor))&(!is.numeric(inter_supp_cor))){
    stop('"inter_supp_cor" must be numeric.')
  }
  if(is.null(inter_supp_cor)&(!is.null(inter_cor_supp))){
    inter_supp_cor <- inter_cor_supp
  }
  if((!is.null(inter_supp_cor))&is.null(inter_cor_supp)){
    inter_cor_supp <- inter_supp_cor
  }
  #
  if(!is.character(output_file)){
    stop('"output_file" incorrectly specified.')
  }
  if(!dir.exists(dirname(output_file))){
    stop('"output_file incorrectly specified."')
  }
  #
  if(!inherits(cor_layer,"igraph")){
    stop('"cor_layer" must be an "igraph" object.')
  }
  if(!inherits(supplementary_layer,"igraph")){
    stop('"supplementary_layer" must be an "igraph" object.')
  }
  if(length(E(cor_layer))==0){
    stop('No edges in "cor_layer".')
  }
  if(length(E(supplementary_layer))==0){
    stop('No edges in "supplementary_layer".')
  }
  # check for attributes in cor_layer
  cor_attr <- "Cor"
  if(!(cor_attr %in% edge_attr_names(cor_layer))){
    if(cor_weighted&(length(edge_attr_names(cor_layer))==0)){
      stop("No weight attribute given.")
    }else if((!cor_weighted)&(length(edge_attr_names(cor_layer))==0)){
      cor_layer <- set_edge_attr(cor_layer, cor_attr, value=1)
    }else{
      cor_attr <- edge_attr_names(cor_layer)[1]
    }
  }
  if(!is.numeric(edge_attr(cor_layer, cor_attr))){
    stop('Correlation attribute of "cor_layer" must be numeric.')
  }
  # check for attributes in supplementary_layer (we can't set unweighted to 1 because we still need to check negative edges)
  dist_attr <- "Dist"
  if(!(dist_attr %in% edge_attr_names(supplementary_layer))){
    if(supp_weighted&(length(edge_attr_names(supplementary_layer))==0)){
      stop("No weight attribute given.")
    }else if((!supp_weighted)&(length(edge_attr_names(supplementary_layer))==0)){
      supplementary_layer <- set_edge_attr(supplementary_layer, dist_attr, value=1)
    }else{
      dist_attr <- edge_attr_names(supplementary_layer)[1]
    }
  }
  if(!is.numeric(edge_attr(supplementary_layer, dist_attr))){
    stop('Distance attribute of "supplementary_layer" must be numeric.')
  }
  if(any(edge_attr(supplementary_layer, dist_attr)<0)){
    stop('Distance attribute of "supplementary_layer" must be non-negative.')
  }
  # remove negative correlations in both layers
  cor_edges_to_remove <- E(cor_layer)[edge_attr(cor_layer, cor_attr)<0]
  supp_edges_to_remove <- .intersect_edges(cor_layer, supplementary_layer, cor_edges_to_remove)
  cor_layer <- delete_edges(cor_layer, cor_edges_to_remove)
  supplementary_layer <- delete_edges(supplementary_layer, supp_edges_to_remove)
  # save and remove list of nodes isolated in both layers
  removed_nodes <- intersect(V(cor_layer)[which(degree(cor_layer)==0)]$name,V(supplementary_layer)[which(degree(supplementary_layer)==0)]$name)
  if(length(removed_nodes)!=0){
    cor_layer <- delete_vertices(cor_layer, V(cor_layer)[removed_nodes])
    supplementary_layer <- delete_vertices(supplementary_layer, V(supplementary_layer)[removed_nodes])
  }
  if(length(E(cor_layer))==0){
    stop('After the removal of negative edges, no edges left in "cor_layer".')
  }
  if(length(E(supplementary_layer))==0){
    stop('After the removal of negative edges, no edges left in "supplementary_layer".')
  }
  # change weights: weighted/unweighted + normalization
  if(!cor_weighted){
    cor_layer <- set_edge_attr(cor_layer, cor_attr, value=1)
  }else if(!is.null(cor_normalization_fun)){
    cor_layer <- set_edge_attr(cor_layer, cor_attr, value=cor_normalization_fun(edge_attr(cor_layer, cor_attr)))
  }
  if(!supp_weighted){
    supplementary_layer <- set_edge_attr(supplementary_layer, dist_attr, value=1)
  }else if(!is.null(supp_normalization_fun)){
    supplementary_layer <- set_edge_attr(supplementary_layer, dist_attr, value=supp_normalization_fun(edge_attr(supplementary_layer, dist_attr)))
  }
  # output file
  list_of_lines <- character()
  if(output_type=="infomap"){
    # output based on "infomap" https://www.mapequation.org/infomap/
    list_of_lines <- c("# A multilayer methylation network")
    # vertices
    physical_nodes <- union(V(cor_layer)$name,V(supplementary_layer)$name)
    list_of_lines <- c(list_of_lines, paste("*Vertices", length(physical_nodes)), "# node_id name", paste0(1:length(physical_nodes),' "',physical_nodes,'"'))
    # intra-layer edges
    list_of_lines <- c(list_of_lines, "*Intra", "# layer_id node_id node_id weight")
    cor_layer_node1_in_edge <- sapply(E(cor_layer), function(x) which(physical_nodes==head_of(cor_layer,x)$name))
    cor_layer_node2_in_edge <- sapply(E(cor_layer), function(x) which(physical_nodes==tail_of(cor_layer,x)$name))
    list_of_lines <- c(list_of_lines, paste("1", cor_layer_node1_in_edge, cor_layer_node2_in_edge, edge_attr(cor_layer, cor_attr)), paste("1", cor_layer_node2_in_edge, cor_layer_node1_in_edge, edge_attr(cor_layer, cor_attr)))
    supplementary_layer_node1_in_edge <- sapply(E(supplementary_layer), function(x) which(physical_nodes==head_of(supplementary_layer,x)$name))
    supplementary_layer_node2_in_edge <- sapply(E(supplementary_layer), function(x) which(physical_nodes==tail_of(supplementary_layer,x)$name))
    list_of_lines <- c(list_of_lines, paste("2", supplementary_layer_node1_in_edge, supplementary_layer_node2_in_edge, edge_attr(supplementary_layer, dist_attr)), paste("2", supplementary_layer_node2_in_edge, supplementary_layer_node1_in_edge, edge_attr(supplementary_layer, dist_attr)))
    # inter-layer edges
    # we add inter-layer edges only if inter_cor_supp and inter_supp_cor are specified
    if(!is.null(inter_cor_supp)){
      list_of_lines <- c(list_of_lines, "*Inter", "# layer_id node_id layer_id weight")
      intersection_nodes <- which(physical_nodes %in% intersect(V(cor_layer)$name,V(supplementary_layer)$name))
      list_of_lines <- c(list_of_lines, paste("1", intersection_nodes, "2", inter_cor_supp), paste("2", intersection_nodes, "1", inter_supp_cor))
    }
    #######################################################################
  }else{
    # output based on "multinet" R package
    if(is.null(inter_cor_supp)){
      layers_multinet <- c("Correlation", "Supplementary")
      list_of_lines <- c("#TYPE", "multiplex")
      list_of_lines <- c(list_of_lines, "#LAYERS", paste(layers_multinet,"UNDIRECTED",sep=","))
      if(cor_weighted | supp_weighted){
        list_of_lines <- c(list_of_lines, "#EDGE ATTRIBUTES", paste(layers_multinet[c(cor_weighted,supp_weighted)],"weight","NUMERIC",sep=","))
      }
      # edges
      list_of_lines <- c(list_of_lines, "#EDGES")
      cor_layer_node1_in_edge <- sapply(E(cor_layer), function(x) head_of(cor_layer,x)$name)
      cor_layer_node2_in_edge <- sapply(E(cor_layer), function(x) tail_of(cor_layer,x)$name)
      if(cor_weighted){
        list_of_lines <- c(list_of_lines, paste(cor_layer_node1_in_edge, cor_layer_node2_in_edge, layers_multinet[1], edge_attr(cor_layer, cor_attr), sep=","))
      }else{
        list_of_lines <- c(list_of_lines, paste(cor_layer_node1_in_edge, cor_layer_node2_in_edge, layers_multinet[1]))
      }
      supplementary_layer_node1_in_edge <- sapply(E(supplementary_layer), function(x) head_of(supplementary_layer,x)$name)
      supplementary_layer_node2_in_edge <- sapply(E(supplementary_layer), function(x) tail_of(supplementary_layer,x)$name)
      if(supp_weighted){
        list_of_lines <- c(list_of_lines, paste(supplementary_layer_node1_in_edge, supplementary_layer_node2_in_edge, layers_multinet[2], edge_attr(supplementary_layer, dist_attr), sep=","))
      }else{
        list_of_lines <- c(list_of_lines, paste(supplementary_layer_node1_in_edge, supplementary_layer_node2_in_edge, layers_multinet[2]))
      }
    }else{
      list_of_lines <- c("#TYPE", "multilayer")
      list_of_lines <- c(list_of_lines, "#LAYERS", paste("Correlation","Correlation","UNDIRECTED",sep=","), paste("Supplementary","Supplementary","UNDIRECTED",sep=","), paste("Correlation","Supplementary","DIRECTED",sep=","), paste("Supplementary","Correlation","DIRECTED",sep=","))
      list_of_lines <- c(list_of_lines, "#EDGE ATTRIBUTES", "weight,NUMERIC") #global attribute
      # edges
      list_of_lines <- c(list_of_lines, "#EDGES")
      cor_layer_node1_in_edge <- sapply(E(cor_layer), function(x) head_of(cor_layer,x)$name)
      cor_layer_node2_in_edge <- sapply(E(cor_layer), function(x) tail_of(cor_layer,x)$name)
      list_of_lines <- c(list_of_lines, paste(cor_layer_node1_in_edge, layers_multinet[1], cor_layer_node2_in_edge, layers_multinet[1], edge_attr(cor_layer, cor_attr), sep=","))
      supplementary_layer_node1_in_edge <- sapply(E(supplementary_layer), function(x) head_of(supplementary_layer,x)$name)
      supplementary_layer_node2_in_edge <- sapply(E(supplementary_layer), function(x) tail_of(supplementary_layer,x)$name)
      list_of_lines <- c(list_of_lines, paste(supplementary_layer_node1_in_edge, layers_multinet[2], supplementary_layer_node2_in_edge, layers_multinet[2], edge_attr(supplementary_layer, dist_attr), sep=","))
      # inter-layer edges
      intersection_nodes <- intersect(V(cor_layer)$name,V(supplementary_layer)$name)
      list_of_lines <- c(list_of_lines, paste(intersection_nodes, layers_multinet[1], intersection_nodes, layers_multinet[2], inter_cor_other, sep=","))
      list_of_lines <- c(list_of_lines, paste(intersection_nodes, layers_multinet[2], intersection_nodes, layers_multinet[1], inter_other_cor, sep=","))
    }
  }
  #
  # save to file
  fileConn <- file(output_file)
  writeLines(list_of_lines, fileConn)
  close(fileConn)
  #
  return(removed_nodes)
}

