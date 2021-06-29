#' Sorting CpGs by coordinate
#' 
#' @description Sorts vector of CpG sites based on their location. 
#' All CpGs are firstly sorted by the chromosome and then
#' by the position on the chromosome. 
#' If working with Illumina Infinium HM450K, `cg_meta` and `cg_meta_cols` 
#' parameters should keep their default values.
#' 
#' @param cg_list Vector of CpG sites.
#' @param cg_meta Data frame with CpGs in the rows and their location in the 
#' columns. The columns of `cg_meta` should include CpG IDs, corresponding 
#' chromosome and coordinate.
#' By default, function uses `CpG_anno450K` data frame which contains Illumina 
#' Infinium HumanMethylation450 manifest file.
#' @param cg_meta_cols Named list with `cg_meta` column names. The list must include:
#' `cg_id` naming the column with unique CpG names, `cg_chr` naming the column 
#' with chromosomes and `cg_coord` naming the column with coordinates.
#' Default value shouldn't be changed if `cg_meta` keeps it default value.
#' 
#' @return Sorted vector of CpGs.
#' 
#' @export
coord_sort <- function(cg_list, cg_meta=data("CpG_anno450K", package="meNet"), cg_meta_cols=list(cg_id="IlmnID", cg_coord="MAPINFO", cg_chr="CHR")){
  cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "cg_coord", "cg_chr"), index_column="cg_id")
  if(!all(cg_list %in% cg_meta[,cg_meta_cols$cg_id])){
    stop('Elements of "cg_list" not found in "cg_meta$cg_id".')
  }
  cg_meta <- cg_meta[cg_meta[,cg_meta_cols$cg_id]%in%cg_list,]
  cg_meta <- cg_meta[order(cg_meta[,cg_meta_cols$cg_chr], cg_meta[,cg_meta_cols$cg_coord]),]
  permutation <- sapply(cg_list, function(x) which(cg_meta[,cg_meta_cols$cg_id]==x), USE.NAMES=FALSE)
  return(permutation) #used by igraph::permute
}