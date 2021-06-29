#' Multiplex plot
#'
#' @description Plots a multiplex network from two `igraph` layers. It is a 
#' corresponding plotting function for `meMultiplex` function
#' which plots two layers next to each other preserving the location of 
#' shared nodes.
#' As in `meMultiplex`, the first layer can contain negative edge weights while 
#' the second layer must have only non-negative weights.
#' 
#' @param first_layer The first layer of the multiplex as `igraph` object.
#' @param second_layer The second layer of the multiplex as `igraph` object.
#' @param layout_layer Defines the plotting method used to determines the layout
#' of nodes. Can take values `"first"`, `"second"`, `"module"`, `"aggregated"` or 
#' `"coord"`. See 'Details'.
#' @param module_df Data frame with community structure. The first column 
#' contains names of nodes and the second column contains corresponding
#' community. It is mandatory parameter if `layout_layer="module"` or
#' `plot_module=TRUE`.
#' @param first_weighted Whether the first layer is weighted. If `FALSE`,
#' all weights are set to `1`. Default to `TRUE`.
#' @param second_weighted Whether the second layer is weighted. If `FALSE`,
#' all weights are set to `1`. Default to `TRUE`.
#' @param first_normalization_fun Normalization function applied on the weights
#' of the first layer, if the layer is weighted. If `NULL`, no normalization
#' is applied. Defaults to `max_normalization`.
#' @param second_normalization_fun Normalization function applied on the weights
#' of the second layer, if the layer is weighted. If `NULL`, no normalization
#' is applied. Defaults to `neg_max_normalization`.
#' @param node_sort_fun Sort function used on the nodes of the network.
#' Default to `.coord_sort`, which sorts CpGs based on their genomic location.
#' @param layout_fun Algorithm for the layout of the nodes used on 
#' the layout graph.
#' In choosing of the function, attention should be paid whether weights of 
#' the layout graph represent "strength" or "cost".
#' Default to `igraph::layout_with_fr`.
#' @param layout_weighted Whether the edge weights of the layout graph are
#' used by `layout_fun` to determine the layout of nodes.
#' Default to `FALSE`.
#' @param layout_weight_transformation_fun Normalization function applied on the 
#' weights of the layout graph, if the graph is weighted. If `NULL`, no normalization
#' is applied. Defaults to `max_normalization`.
#' @param main_title Main title of the plot. By default, there is no main title.
#' @param first_title Title of the first subplot with the first layer of 
#' multiplex. Default value is `"Correlation"`.
#' @param second_title Title of the second subplot with the second layer of 
#' multiplex. Default value is `"CpG island"`.
#' @param node_col Color of the nodes. Color name if all nodes have the same
#' color or data frame with node names in the first column and corresponding 
#' color in the second column. If some nodes are missing from the given data
#' frame, their color is set to `"black"`.
#' Defaults value of `node_col` is `"black"`.
#' @param node_size Size of the nodes. Defaults to `2`.
#' @param pos_edge_col Color of positive edges. Default value is `"grey70"`.
#' @param neg_edge_col Color of negative edges. Default value is `"blue"`.
#' @param plot_module Whether to plot the border around the communities defined
#' in `module_df`.
#' Defaults to `FALSE`.
#' @param module_border_col Color of the border around the communities. Passed
#' to `plot.igraph` as `mark.border` parameter.
#' Defaults to `"black"`.
#' @param module_border_expand Size of the border around the communities. Passed
#' to `plot.igraph` as `mark.expand` parameter.
#' Defaults to `10`.
#' @param plot_layout_graph Whether to add a plot of layout graph. If `TRUE`, third
#' plot is added to the right of the multiplex with the same layout of nodes and
#' with edges from the layout graph. Default to `FALSE`.
#' 
#' @return A vector of removed nodes.
#' 
#' @details 
#' Plots two-layer multiplex from two `igraph` layers. Layers are plotted
#' next to each other with the same layout of nodes. Plotted edges represent
#' edges found in a certain layer of multiplex. Parameters `main_title`,
#' `first_title`, `second_title`, `node_col`, `node_size`, `pos_edge_col` and
#' `neg_edge_col` help in customization of the plot.
#'
#' Nodes which are isolated in both layers are removed from the multiplex. Nodes
#' which are found only in one layer are added to the other layer as isolated 
#' nodes. This equalizes the set of nodes between two layers.
#' 
#' If `module_df` is given, borders can be plotted around the communities of
#' nodes setting `plot_module` to `TRUE`.
#' If the same node belong to multiple communities, only the first 
#' membership will be used in plotting. For customization, parameters
#' `module_border_col` and `module_border_expand` can be used.
#' 
#' To determine the layout of nodes, new layout graph is constructed and used
#' in the `layout_fun`. The nodes of the layout graph are the same as in
#' the first and the second layer. 
#' If `layout_weighted`, the weight of edges are used by
#' `layout_fun`. There are several option for the construction of the layout
#' graph, depending on the value of `layout_layer`. If `layout_layer` has value
#' \itemize{
#'    \item{`"first"`}{layout graph is the first layer with removed negative
#'    edges;}
#'    \item{`"second"`}{Layout graph is the second layer.}
#'    \item{`"module"`}{Layout graph is constructed from `module_df`. In layout
#'    graph edges exist between nodes which are found in the same community.
#'    Weight of an edge is the number of communities in which a selected pair
#'    of nodes is found. Weight can be larger than `1` only if the same node
#'    can be a member of multiple communities in `module_df`.}
#'    \item{`"aggregated"`}{Layout graph is obtained aggregating the first and
#'    the second layer. Weight of an edge in layout graph is the sum of weights
#'    of the same edge in the first and the second layer.}
#'    \item{`"coord"`}{Layout graph is disconnected graph with no edges. This
#'    only makes sense if `layout_fun` uses circular layout since then
#'    nodes will be ordered by `node_sort_fun`.}
#' }
#' Structure of the layout graph can be plotted as third layer setting
#' `plot_layout_graph` to `TRUE`. It if not advisable to use this in the final
#' plot, but this feature can be useful to find the appropriate `layout_fun`.
#' 
#' If layers are weighted, function searches for the edge attributes with
#' names`"Cor"` and `"Dist"`. If the attribute doesn't exist, it takes the first
#' edge attribute instead. The function handles weight of 
#' the first layer as correlations and weights of the second layer as distances. 
#' Function checks if all edge weights of `second_layer` are positive.
#' 
#' 
#' @import igraph
#' 
#' @export
plot_multiplex <- function(first_layer, second_layer, layout_layer="first", module_df=NULL, first_weighted=TRUE, second_weighted=TRUE,
                           first_normalization_fun=.max_normalization, second_normalization_fun=.neg_max_normalization, include_negative=TRUE,
                           node_sort_fun=.coord_sort, layout_fun=igraph::layout_with_fr, layout_weighted=FALSE,
                           layout_weight_transformation_fun=.max_normalization,
                           main_title=NULL, first_title="Correlation", second_title="CpG island",
                           node_col="black", node_size=2, pos_edge_col="grey70", neg_edge_col="blue",
                           plot_module=FALSE, module_border_col="black", module_border_expand=10, plot_layout_graph=FALSE){
  if(!(layout_layer) %in% c("first", "second", "module", "aggregated", "coord")){
    stop('"layout_layer" must have one of the values: "first", "second", "module", "aggregated" or "coord".')
  }
  #
  if(!inherits(first_layer,"igraph")){
    stop('"first_layer" must be an "igraph" object.')
  }
  if(!inherits(second_layer,"igraph")){
    stop('"second_layer" must be an "igraph" object.')
  }
  if(length(E(first_layer))==0){
    stop('No edges in "first_layer".')
  }
  if(length(E(second_layer))==0){
    stop('No edges in "first_layer".')
  }
  # check for module
  if(layout_layer=="module" | plot_module){
    if(is.null(module_df)){
      if(layout_layer=="module"){
        stop('For the "module" layout of nodes, "module_df" must be provided.')
      }else{
        stop('To plot module borders, "module_df" must be provided.')
      }
    }
    if(!inherits(module_df, "data.frame")){
      stop('"module_df" must be data frame.')
    }
    if(ncol(module_df)!=2){
      stop('"module_df" must have two columns, the first with node IDs and the second with corresponding community membership.')
    }
    if(length(union(intersect(V(first_layer)$name, module_df[,1]),intersect(V(second_layer)$name, module_df[,1])))==0){
      stop('The first column of "module_df" has no shared elements with the nodes. Provide different data frame or adjust other function parameters so that "module_df" is not required.')
    }
  }
  # check node_col
  if(!is.character(node_col)){
    if(is.null(node_col)){
      warning('"node_col" misspecified, it is set to "black".')
      node_col <- "black"
    }else if(!inherits(node_col, "data.frame")){
      warning('If node color is not the same for all nodes, "node_col" must be a data frame with node IDs in the first column and corresponding color in the second column. Node color is set to "black".')
      node_col <- "black"
    }else if(ncol(node_col)!=2){
      warning('If node color is not the same for all nodes, "node_col" must have two columns: the first with node IDs and the second with corresponding color. Node color is set to "black".')
      node_col <- "black"
    }
  }
  ## check for attributes in cor_layer
  first_attr <- "Cor"
  if(!(first_attr %in% edge_attr_names(first_layer))){
    if(first_weighted&(length(edge_attr_names(first_layer))==0)){
      stop("No weight attribute given.")
    }else if((!first_weighted)&(length(edge_attr_names(first_layer))==0)){
      cor_layer <- set_edge_attr(first_layer, first_attr, value=1)
    }else{
      first_attr <- edge_attr_names(first_layer)[1]
    }
  }
  if(!is.numeric(edge_attr(first_layer, first_attr))){
    stop('The first attribute of "first_layer" must be numeric.')
  }
  # check for attributes in supplementary_layer (we can't set unweighted to 1 because we still need to check negative edges)
  second_attr <- "Dist"
  if(!(second_attr %in% edge_attr_names(second_layer))){
    if(second_weighted&(length(edge_attr_names(second_layer))==0)){
      stop("No weight attribute given.")
    }else if((!second_weighted)&(length(edge_attr_names(second_layer))==0)){
      second_layer <- set_edge_attr(second_layer, second_attr, value=1)
    }else{
      second_attr <- edge_attr_names(second_layer)[1]
    }
  }
  if(!is.numeric(edge_attr(second_layer, second_attr))){
    stop('The first attribute of "second_layer" must be numeric.')
  }
  if(any(edge_attr(second_layer, second_attr)<0)){
    stop('The first attribute of "second_layer" must be non-negative.')
  }
  # remove nodes which are isolated in both layers
  removed_nodes <- intersect(V(first_layer)[which(degree(first_layer)==0)]$name,V(second_layer)[which(degree(second_layer)==0)]$name)
  if(length(removed_nodes)!=0){
    first_layer <- delete_vertices(first_layer, V(first_layer)[removed_nodes])
    second_layer <- delete_vertices(second_layer, V(second_layer)[removed_nodes])
  }
  # make set of nodes equal for both layers (and in the same order)
  second_min_first <- setdiff(V(second_layer)$name, V(first_layer)$name)
  first_min_second <- setdiff(V(first_layer)$name, V(second_layer)$name)
  first_layer <- add_vertices(first_layer, length(second_min_first), name=second_min_first)
  second_layer <- add_vertices(second_layer, length(first_min_second), name=first_min_second)
  first_layer <- igraph::permute(first_layer, node_sort_fun(V(first_layer)$name))
  second_layer <- igraph::permute(second_layer, node_sort_fun(V(second_layer)$name))
  # set color attribute
  if(!is.character(node_col)){
    if(length(setdiff(V(first_layer)$name, node_col[,1]))>0){
      warning(paste0('The first column of "node_col" ',"doesn't", ' contain all nodes. Node color for missing nodes is set to "black".'))
      missing_nodes <- setdiff(V(first_layer)$name, node_col[,1])
      tryCatch({
        node_col_add <- data.frame(matrix(nrow=length(missing_nodes),ncol=2))
        colnames(node_col_add) <- colnames(node_col)
        node_col_add[,1] <- missing_nodes
        node_col_add[,2] <- "black"
        node_col_df <- rbind(node_col, node_col_add)
      }, error=function(e){
        stop('"node_col" data frame is misspecified. The first column should contain node IDs and the second column corresponding color.')
      })
    }else{
      node_col_df <- node_col
    }
    if(sum(duplicated(node_col_df[,1]))>0){
      warning('"node_col" has duplicated node IDs. Only the first entry is used to set color.')
      node_col_df <- node_col_df[!duplicated(node_col_df[,1]),]
    }
    rownames(node_col_df) <- node_col_df[,1]
  }else{
    node_col_df <- data.frame(matrix(nrow=length(V(first_layer)), ncol=2))
    node_col_df[,1] <- V(first_layer)$name
    node_col_df[,2] <- node_col
    rownames(node_col_df) <- node_col_df[,1]
  }
  # change weights: weighted/unweighted + normalization
  if(!first_weighted){
    first_layer <- set_edge_attr(first_layer, first_attr, value=1)
  }else if(!is.null(first_normalization_fun)){
    first_layer <- set_edge_attr(first_layer, first_attr, value=first_normalization_fun(edge_attr(first_layer, first_attr)))
  }
  if(!second_weighted){
    second_layer <- set_edge_attr(second_layer, second_attr, value=1)
  }else if(!is.null(second_normalization_fun)){
    second_layer <- set_edge_attr(second_layer, second_attr, value=second_normalization_fun(edge_attr(second_layer, second_attr)))
  }
  # LAYOUT
  if(layout_layer=="first"){
    layout_graph <- delete.edges(first_layer, which(edge_attr(first_layer, first_attr)<0))
    layout_attr <- first_attr
  }else if(layout_layer=="second"){
    layout_graph <- second_layer
    layout_attr <- second_attr
  }else if(layout_layer=="module"){
    module_df <- module_df[module_df[,1]%in%V(first_layer)$name,]
    module_df <- module_df[!duplicated(module_df),]
    modules <- unique(module_df[,2])
    edge_df <- data.frame(matrix(nrow=0,ncol=2))
    colnames(edge_df) <- c("Node1", "Node2")
    for(i in 1:length(modules)){
      cg_list_i <- unique(module_df[module_df[,2]==modules[i],1])
      edge_df <- rbind(edge_df, .get_all_pairwiseCG(cg_list_i, "Node1", "Node2"))
    }
    edge_df[,1:2] <- t(apply(edge_df[,1:2], 1, sort)) #in case nodes are in different order
    edge_df <- plyr::ddply(edge_df,.(Node1,Node2),nrow) # we count in how many modules is the pair present
    colnames(edge_df) <- c("Node1", "Node2", "weight")
    node_df <- data.frame(IlmnID=V(first_layer)$name)
    layout_graph <- graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE)
    layout_attr <- "weight"
  }else if(layout_layer=="aggregated"){
    # new graph which aggregates two layers
    edge_df1 <- as.data.frame(t(sapply(E(first_layer), function(x) c(head_of(first_layer,x)$name, tail_of(first_layer,x)$name))))
    edge_df1 <- cbind(edge_df1,edge_attr(first_layer,first_attr))
    colnames(edge_df1) <- c("Node1", "Node2", "weight")
    edge_df2 <- as.data.frame(t(sapply(E(second_layer), function(x) c(head_of(second_layer,x)$name, tail_of(second_layer,x)$name))))
    edge_df2 <- cbind(edge_df2,edge_attr(second_layer,second_attr))
    colnames(edge_df2) <- c("Node1", "Node2", "weight")
    edge_df <- rbind(edge_df1, edge_df2)
    edge_df[,1:2] <- t(apply(edge_df[,1:2], 1, sort)) # to put nodes in the same order
    edge_df_intersect <- merge(edge_df[duplicated(edge_df[,1:2]),],edge_df[duplicated(edge_df[,1:2], fromLast=TRUE),], by=c("Node1", "Node2"))
    edge_df_intersect$weight <- edge_df_intersect[,3] + edge_df_intersect[,4]
    edge_df_intersect <- edge_df_intersect[,-c(3,4)]
    edge_df <- rbind(edge_df[!duplicated(edge_df[,1:2]),],edge_df_intersect)
    node_df <- data.frame(IlmnID=V(first_layer)$name)
    layout_graph <- graph_from_data_frame(d=edge_df, vertices=node_df, directed=FALSE)
    layout_attr <- "weight"
  }else{
    # disconnected graph, only makes sense if circular layout is used since vertices are ordered by coordinates
    layout_graph <- make_empty_graph(n=0, directed=FALSE)
    layout_graph <- add_vertices(layout_graph, lenth(V(first_layer)), name=V(first_layer)$name)
    layout_attr <- ""
  }
  if(!is.null(layout_weight_transformation_fun)&(layout_attr!="")&layout_weighted){
    layout_graph <- set_edge_attr(layout_graph, layout_attr, value=layout_weight_transformation_fun(edge_attr(layout_graph, layout_attr)))
  }
  # calculate layout
  if(layout_weighted&(layout_layer!="coord")){
    tryCatch({
      l <- layout_fun(layout_graph, weights=edge_attr(layout_graph, layout_attr))
    }, error=function(e){
      warning(paste0('"layout_fun" ', "doesn't accept weights."))
      l <- layout_fun(layout_graph)
    })
  }else{
    tryCatch({
      l <- layout_fun(layout_graph, weights=NA) #NA skips weight even if the graph has "weight" argument
    }, error=function(e){
      l <- layout_fun(layout_graph)
    })
  }
  # PLOTTING
  if(plot_layout_graph){
    par(mfrow=c(1,3))
  }else{
    par(mfrow=c(1,2))
  }
  if(!is.null(main_title)){
    par(oma=c(0,0,2,0))
    mtext(main_title, line=0, side=3, outer=TRUE, cex=2)
  }
  if(plot_module){
    # clean module_df
    module_df <- module_df[module_df[,1]%in%V(first_layer)$name,]
    if(sum(duplicated(module_df[,1]))>0){
      warning('Some nodes belong to multiple modules. For plotting of modules, only the first membership is used.')
      module_df <- module_df[!duplicated(module_df[,1]),]
    }
    module_df[,2] <- as.numeric(factor(module_df[,2])) # so that module names are integers starting from 1
    no_module_nodes <- setdiff(V(first_layer)$name, module_df[,1])
    if(length(no_module_nodes)>0){
      module_df_add <- data.frame(matrix(nrow=length(no_module_nodes), ncol=2))
      colnames(module_df_add) <- colnames(module_df)
      module_df_add[,1] <- no_module_nodes
      module_df_add[,2] <- NA
      module_df <- rbind(module_df, module_df_add)
    }
    rownames(module_df) <- module_df[,1]
    # create igraph communities objects
    first_communities <- make_clusters(first_layer, membership=module_df[V(first_layer)$name,2], modularity=FALSE)
    second_communities <- make_clusters(second_layer, membership=module_df[V(second_layer)$name,2], modularity=FALSE)
    # plot
    plot(first_communities, first_layer, vertex.label=NA, layout=l, vertex.size=node_size, mark.expand=module_border_expand, mark.border=module_border_col, mark.col=NA, col=node_col_df[V(first_layer)$name,2], edge.color=c(pos_edge_color, neg_edge_color)[(edge_attr(first_layer, first_attr)<0)+1], main=first_title)
    plot(second_communities, second_layer, vertex.label=NA, layout=l, vertex.size=node_size, mark.expand=module_border_expand, mark.border=module_border_col, mark.col=NA, col=node_col_df[V(second_layer)$name,2], edge.color=pos_edge_color, main=second_title)
    if(plot_layout_graph){
      layout_communities <- make_clusters(layout_layer, membership=module_df[V(layout_layer)$name,2], modularity=FALSE)
      plot(layout_communities, layout_graph, vertex.label=NA, layout=l, vertex.size=node_size, mark.expand=module_border_expand, mark.border=module_border_col, mark.col=NA, col=node_col_df[V(layout_graph)$name,2], edge.color=pos_edge_color, main="Layout")
    }
  }else{
    plot(first_layer, vertex.label=NA, layout=l, vertex.size=node_size, vertex.color=node_col_df[V(first_layer)$name,2], edge.color=c(pos_edge_color, neg_edge_color)[(edge_attr(first_layer, first_attr)<0)+1], main=first_title)
    plot(second_layer, vertex.label=NA, layout=l, vertex.size=node_size, vertex.color=node_col_df[V(second_layer)$name,2], edge.color=pos_edge_color, main=second_title)
    if(plot_layout_graph){
      plot(layout_graph, vertex.label=NA, layout=l, vertex.size=node_size, vertex.color=node_col_df[V(layout_graph)$name,2], edge.color=pos_edge_color, main="Layout")
    }
  }
  return(removed_nodes)
}



