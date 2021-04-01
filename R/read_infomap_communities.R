# physical or state nodes
#'@export
read_infomap_communities <- function(multiplex_file, clu_file, physical_nodes=TRUE, layer=NULL){
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
