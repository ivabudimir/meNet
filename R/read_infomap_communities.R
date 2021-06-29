#' Multiplex community structure from file
#' 
#' @description 
#' Reads `Infomap` community structure for multiplex from two files.
#' The first file contains `Infomap`-style multiplex structure and the second
#' file is an output from the `Infomap` algorithm containing the communities.
#' Function his returns data frame with nodes and their corresponding community.
#' 
#' @param multiplex_file File with the `Infomap` structure of the multiplex.
#' @param clu_file File with the communities obtained as output of the 
#' `Infomap` algorithm.
#' @param physical_nodes Whether the `clu_file` contains physical of state
#' nodes. Default to `FALSE`.
#' @param layer The layer from which the community structure is obtained.
#' Only used if `clu_file` contains state nodes.
#' 
#' @return Data frame with two columns. The first column contains names of nodes
#' and the second column contains corresponding `Infomap` community.
#' 
#' @details 
#' The first file can be obtained with `meMultiplex` function with the
#' output type set to `"Infomap"`. The second file is a result of `Infomap` 
#' algorithm for community detection.
#' `Infomap` algorithm finds community structure of the multiplex using
#' `multiplex_file` and outputs it in `clu_file`. 
#' This function unites the two files in one tidy data frame with multiplex 
#' nodes and their corresponding community.
#' 
#' `Infomap` returns two different types of communities, communities based on 
#' physical nodes (normal nodes) and communities based on state nodes. State 
#' nodes are abstract nodes created by `Infomap` to label physical nodes which 
#' exist in multiple layers. For example, for a physical node which is shared 
#' between two layers, two state nodes are created. Thus community structure of 
#' the state nodes contains two community structures, one existing in each layer 
#' of multiplex and both od them influenced by the whole multiplex structure.
#' If `physical_nodes` is `FALSE`, `clu_file` must contain state node output.
#' In this case `layer` must be specified and the community structure is taken
#' from that layer.
#' If `physcal_nodes` is `TRUE`, `clu_file` must contain physical node output.
#' This output overlaps two state node community structures existing in different
#' layers and thus the same physical node can exist in multiple communities.
#' For more information, visit `Infomap` website \insertCite{infomap}{meNet}.
#' 
#' @references
#'       \insertAllCited{}
#'
# '@export
read_infomap_communities <- function(multiplex_file, clu_file, physical_nodes=FALSE, layer=NULL){
  if(!is.character(multiplex_file)){
    stop('"multiplex_file" incorrectly specified.')
  }
  if(!file.exists(multiplex_file)){
    stop(paste0('"multiplex_file"', " doesn't exist."))
  }
  if(!is.character(clu_file)){
    stop('"clu_file" incorrectly specified.')
  }
  if(!file.exists(clu_file)){
    stop(paste0('"clu_file"', " doesn't exist."))
  }
  if(!physical_nodes){
    if(!is.numeric(layer)){
      stop('In case of state nodes, layer must be specified.')
    }
  }
  #
  fileConn <- file(multiplex_file, "r")
  multiplex_df <- data.frame(matrix(nrow=0, ncol=2))
  colnames(multiplex_df) <- c("node_id", "IlmnID")
  read_more = TRUE
  while(read_more){
    line = readLines(fileConn, n = 1)
    if(length(line)==0){
      next
    }
    if(grepl("*Vertices", line))
    {
      n_nodes <- as.numeric(strsplit(line, " ")[[1]][2])
      multiplex_df <- data.frame(node_id=1:n_nodes, IlmnID="")
      line_after <- readLines(fileConn, n=1)
      while(grepl("#", line_after)){ #check if there are comments
        line_after <- readLines(fileConn, n=1)
      }
      multiplex_df[1,] <- c(as.numeric(strsplit(line_after, " ")[[1]][1]), strsplit(strsplit(line_after, " ")[[1]][2], '"')[[1]][2])
      node_lines <- readLines(fileConn, n=n_nodes-1)
      multiplex_df[2:n_nodes,1] <- sapply(node_lines, function(x) as.numeric(strsplit(x, " ")[[1]][1]))
      multiplex_df[2:n_nodes,2] <- sapply(node_lines, function(x) strsplit(strsplit(x, " ")[[1]][2], '"')[[1]][2])
      read_more <- FALSE
    }
  }
  close(fileConn)
  #
  fileConn <- file(clu_file, "r")
  all_lines <- readLines(fileConn)
  all_lines <- all_lines[!grepl("#", all_lines)]
  if(physical_nodes){
    # node_id module flow
    clu_df <- data.frame(node_id=sapply(all_lines, function(x) as.numeric(strsplit(x, " ")[[1]][1]), USE.NAMES=FALSE),
                         Module=sapply(all_lines, function(x) as.numeric(strsplit(x, " ")[[1]][2]), USE.NAMES=FALSE))
  }else{
    # state_id module flow node_id layer_id
    clu_df <- data.frame(node_id=sapply(all_lines, function(x) as.numeric(strsplit(x, " ")[[1]][4]), USE.NAMES=FALSE),
                         Module=sapply(all_lines, function(x) as.numeric(strsplit(x, " ")[[1]][2]), USE.NAMES=FALSE),
                         layer_id=sapply(all_lines, function(x) as.numeric(strsplit(x, " ")[[1]][5]), USE.NAMES=FALSE))
    clu_df <- clu_df[clu_df$layer_id==layer,]
    clu_df <- clu_df[,-c(3)]
  }
  close(fileConn)
  #
  result_df <- merge(multiplex_df, clu_df, by="node_id", how="inner")
  result_df <- result_df[,c("IlmnID", "Module")]
  return(result_df)
}
