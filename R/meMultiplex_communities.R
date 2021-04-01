# unite creation of multiplex for infomap and and reading communities from file  (need installed infomap)
#'@export
meMultiplex_communities <- function(cor_layer, supplementary_layer, physical_nodes=FALSE, layer=1, folder="./meNet/", file_basename="meNet_infomap", infomap_call="infomap",
                                    cor_weighted=TRUE, supp_weighted=TRUE, cor_normalization_fun=max_normalization, supp_normalization_fun=neg_max_normalization,
                                    inter_cor_supp=NULL, inter_supp_cor=inter_cor_supp, infomap_seed=NULL, relaxation_rate=0.15, delete_files=FALSE){
  if(!is.character(folder)){
    stop('"folder" incorrectly specified.')
  }
  dir.create(folder, showWarnings = FALSE)
  if(!is.character(file_basename)){
    stop('"file_basename" incorrectly specified.')
  }
  # check infomap command
  error_message <- 'Infomap command not working. Check installation of "Infomap" and parameter "infomap_call". For more details, see https://www.mapequation.org/infomap/.'
  tryCatch({system(paste(infomap_call,"--help"))}, error=function(e){stop(error_message)}, warning=function(w){stop(error_message)})
  # meMultiplex to file
  meNet_file <- paste0(folder,file_basename,".txt")
  meMultiplex(cor_layer, supplementary_layer, meNet_file, cor_weighted, supp_weighted, cor_normalization_fun, supp_normalization_fun, output_type="infomap", inter_cor_supp, inter_supp_cor)
  infomap_command <- paste(infomap_call, "-i multilayer -o clu")
  if(is.numeric(infomap_seed)){
    infomap_command <- paste(infomap_command, "--seed", infomap_seed)
  }
  if(is.numeric(relaxation_rate)){
    infomap_command <- paste(infomap_command, "--multilayer-relax-rate", relaxation_rate)
  }
  infomap_command <- paste(infomap_command, meNet_file, folder)
  error_message <- "Infomap call not executed correctly."
  tryCatch({system(infomap_command)}, error=function(e){stop(error_message)}, warning=function(w){stop(error_message)})
  #
  
  CLU_file <- paste0(folder, file_basename, ".clu")
  statesCLU_file <- paste0(folder, file_basename, "_states.clu")
  if(physical_nodes){
    communities <- read_infomap_communities(meNet_file, CLU_file, physical_nodes)
  }else{
    communities <- read_infomap_communities(meNet_file, statesCLU_file, physical_nodes, layer)
  }
  if(delete_files){
    is_deleted <- FALSE
    #
    tryCatch({ #linux
      system(paste("rm", meNet_file, CLU_file, statesCLU_file))
      is_deleted <- TRUE}, error=function(e){}, warning=function(w){})
    tryCatch({ #windows
      system(paste("del",meNet_file, CLU_file, statesCLU_file))
      is_deleted <- TRUE}, error=function(e){}, warning=function(w){})
    #
    if(length(setdiff(dir(folder, all.files=TRUE), c(".", "..")))==0){ # only if folder is empty, delete it
      tryCatch({ #linux
        system(paste0("rmdir ",'"',folder,'"'))}, error=function(e){}, warning=function(w){})
      tryCatch({ #windows
        system(paste0("Rmdir ",'"',folder,'"'))}, error=function(e){}, warning=function(w){})
    }
    if(!is_deleted){
      warning("Files couldn't be deleted.")
    }
    #
  }
  return(communities)
}

