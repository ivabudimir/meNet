#' Find CpG sites which are associated with a given CpG island
#
#'@export
CpG_in_CGI <- function(cg_island, cg_meta=data("CpG_anno450K", package="meNet"),
                       cg_meta_cols=list(cg_id="IlmnID", island_name="UCSC_CpG_Islands_Name", island_region="Relation_to_UCSC_CpG_Island"),
                       include_regions=c()){
  if(!inherits(cg_meta, "data.frame")){
    stop('"cg_meta" must be data frame.')
  }
  if(!all(unlist(cg_meta_cols, use.names=FALSE) %in% colnames(cg_meta))){
    stop('"cg_meta_cols" are not columns of "cg_meta".')
  }
  # CpGs which are part of an island based on UCSC annotation
  # by default, we only report CpGs inside the island, but we may add additional regions such as shores and shelves
  island_regions <- c("Island")
  island_regions <- unique(c(island_regions, intersect(include_regions, c('N_Shore','S_Shelf','S_Shore','N_Shelf'))))
  return(cg_meta[(cg_meta[,cg_meta_cols$island_name]==cg_island)&(cg_meta[,cg_meta_cols$island_region]%in%island_regions),cg_meta_cols$cg_id])
}
