#' Multiplex community structure
#'
#' @description 
#' For two `igraph` layers with shared nodes, the function finds `Infomap`
#' communities of the constructed multiplex. Function interprets the first layer 
#' as correlation layer. The second layer, if weighted, should have only 
#' positive weights and can be interpreted as distance layer.
#' 
#' `Infomap` has to be locally installed for this function to work. 
#' Help with installation and alternative ways to run the `Infomap` algorithm 
#' are explained in 'Details'.
#' 
#' @param cor_layer Correlation layer of the multiplex as igraph object.
#' @param supplementary_layer Supplementary layer of the multiples as igraph
#' object.
#' @param physical_nodes Should physical or state nodes be used for community
#' structure nodes. Default to `FALSE`, i.e. to state nodes.
#' @param layer The layer from which the community structure is obtained.
#' Only used for state nodes.
#' @param folder Folder in which the files will be saved. It is automatically
#' created if not already present. Default value is `"./meNet/"`.
#' @param file_basename Base name of the creates files. For different files,
#' different suffix is added to the base name. 
#' Default value is `"meNet_infomap"`.
#' @param infomap_call Path to the `Infomap` on user's system.
#' Defaults to `"infomap"`.
#' @param cor_weighted Whether the correlation layer is weighted. Passed to
#' `meMultiplex`. Defaults to `TRUE`.
#' @param supp_weighted Whether the supplementary layer is weighted. Passed to
#' `meMultiplex`. Defaults to `TRUE`.
#' @param cor_normalization_fun Normalization function for the correlation layer.
#' Passed to `meMultiplex`. Defaults to `max_normalization`.
#' @param supp_normalization_fun Normalization function for the supplementary
#' layer. Passed to `meMultiplex`. Defaults to `neg_max_normalization`.
#' @param inter_cor_supp Weight of the inter-layer edges. Passed to
#' `meMultiplex`. Defaults to `NULL`.
#' @param inter_supp_cor Weight of the inter-layer edges. Passed to
#' `meMultiplex`. By default, the values is equal to `inter_cor_supp`.
#' @param infomap_seed `seed` parameter in `Infomap` call.
#' By default, no seed is set.
#' @param relaxation_rate `multilayer-relax-rate` parameter in `Infomap` call.
#' Default to `0.15`.
#' @param delete_files Should created files be automatically deleted
#' from the user's system.
#' Default to `FALSE`. Changing the parameter to `TRUE` should be done with 
#' caution since it will allow the function to delete files from user's system.
#' 
#' @return Data frame with two columns. The first column contains names of nodes
#' and the second column contains corresponding `Infomap` community.
#' 
#' @details 
#' Function finds communities of nodes (CpGs) performing a sequence of three  
#' action. Firstly, `meMultiplex` function is called and the structure of
#' multiplex if written to a file. The file is used in a call of to
#' `Infomap` which performs the clustering of nodes. The clusters are saved to
#' a file. In the end, function `read_infomap_communities` reads the community
#' structure of the multiplex from the two files and returns a tidy data frame
#' with nodes and their corresponding communities.
#' 
#' For multiplex reconstruction, parameters `cor_layer`, `supplementary_layer`,
#' `cor_weighted`, `supp_weighted`, `cor_normalization_fun`, `supp_normalization_fun`,
#' `inter_cor_supp` and `inter_supp_cor` are passed to the `meMultiplex` with
#' `output_type` set to `"infomap"`.
#' 
#' `Infomap` is called from the command line using the path given in
#' `infomap_call` with flags `-i multilayer -o clu`.
#' If specified, `infomap_seed` is added as `--seed` flag and `relaxation-rate`
#' is added as `--multilayer-relax-rate` flag.
#' For help with the installation of the `Infomap` algorithm, visit
#' `Infomap` website: https://www.mapequation.org/infomap/.
#' If the installation fails, the output of the `meMultiplex` can be manually
#' uploaded online on the `Infomap` website and the resulting '.clu` files can 
#' be used as input to `read_infomap_community` function. For futher details,
#' refer to `Infomap` paper \insertCite{infomap}{meNet}.
#' 
#' To read multiplex communities, function `read_infomap_communities` is called
#' with the location of two files from the previous steps. Parameters 
#' `physical_nodes` and `layer` are passed to the function.
#' 
#' In the first two steps of the function algorithm, files are created and
#' saved in the `folder` using `file_basename`. These files are by-product of
#' the function and can be deleted if the user doesn't want to examine their
#' details. If `delete_files` is `TRUE`, these files will be automatically 
#' deleted from the user's system. If the folder contains only these files,
#' it will also be deleted.
#' 
#' @references
#'       \insertAllCited{}
#'       
#' @import igraph
#' 
#' @export
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

