# only first layer can have negative values
# aggregated layout sums the weights of shared links (keep in mind if changing normalization functions)
# module layout counts in how many modules two nodes are linked
# default parameters set for cor&dist weighted combination
#'@export
plot_multiplex <- function(first_layer, second_layer, module_df=NULL, cg_list=NULL, first_weighted=TRUE, second_weighted=TRUE,
                           first_normalization_fun=.max_normalization, second_normalization_fun=.neg_max_normalization, include_negative=TRUE,
                           include_module=FALSE, node_sort_fun=.coord_sort, layout_layer="first", layout_fun=igraph::layout_with_fr, layout_weighted=FALSE, layout_weight_transformation_fun=.max_normalization,
                           main_title=NULL, first_title="Correlation", second_title="CpG island", pos_link_col="grey70", neg_link_col="red", vertex_size=2, plot_layout_graph=FALSE){
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
  if(layout_layer=="module" | include_module){
    if(is.null(module_df)){
      stop('Based the current setting of the function parameters "module_df" must be provided.')
    }
    if(!inherits(module_df, "data.frame")){
      stop('"module_df" must be data frame.')
    }
    if(length(union(intersect(V(first_layer)$name, module_df[,1]),intersect(V(second_layer)$name, module_df[,1])))==0){
      warning('The first column of "module_df" had no shared elements with the nodes.')
      module_df <- NULL
      include_module <- FALSE
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
  # change weights: weighted/unweighted + normalization
  if(!first_weighted){
    first_layer <- set_edge_attr(first_layer, first_attr, value=1)
  }else{
    first_layer <- set_edge_attr(first_layer, first_attr, value=first_normalization_fun(edge_attr(first_layer, first_attr)))
  }
  if(!second_weighted){
    second_layer <- set_edge_attr(second_layer, second_attr, value=1)
  }else{
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
    if(is.null(module_df)){
      stop('"module_df" not provided. Change layout setting.')
    }
    module_df <- module_df[module_df[,1]%in%V(first_layer)$name,]
    module_df <- module_df[!duplicated(module_df),]
    modules <- unique(module_df[,2])
    edge_df <- data.frame(matrix(nrow=0,ncol=2))
    colnames(edge_df) <- c("Node1", "Node2")
    for(i in 1:length(modules)){
      cg_list_i <- unique(module_df[module_df[,2]==modules[i],1])
      edge_df <- rbind(edge_df, .get_all_pairwiseCG(cg_list_i, "Node1", "Node2"))
    }
    edge_df <- edge_df[!duplicated(edge_df),]
    edge_df <- plyr::ddply(edge_df,.(Node1,Node2),nrow)
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
  if(layout_attr!=""){
    layout_graph <- set_edge_attr(layout_graph, layout_attr, value=layout_weight_transformation_fun(edge_attr(layout_graph, layout_attr)))
  }
  # calculate layout
  if(layout_weighted&(layout_layer!="coord")){
    tryCatch({
      l <- layout_fun(layout_graph, weights=edge_attr(layout_graph, layout_attr))
    }, error=function(x){
      warning(paste0('"layout_fun" ', "doesn't accept weights."))
      l <- layout_fun(layout_graph)
    })
  }else{
    tryCatch({
      l <- layout_fun(layout_graph, weights=NA) #NA skips weight even if the graph has "weight" argument
    }, error=function(x){
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
  plot(first_layer, vertex.label=NA, layout=l, vertex.size=vertex_size, vertex.color="black", edge.color=c(pos_link_col, neg_link_col)[(edge_attr(first_layer, first_attr)<0)+1], main=first_title)
  plot(second_layer, vertex.label=NA, layout=l, vertex.size=vertex_size, vertex.color="black", edge.color="grey70", main=second_title)
  if(plot_layout_graph){
    plot(layout_graph, vertex.label=NA, layout=l, vertex.size=vertex_size, vertex.color="black", edge.color="grey70", main="Layout")
  }
}

