# sorts list of CpGs based on their coordinate (CHR.MAPINFO)
#'@export
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