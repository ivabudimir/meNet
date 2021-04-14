#' Find CpG sites which are associated with a given CpG island
#
#'@export
CpG_in_CGI <- function(cg_island, cg_meta=cg_anno450k,
                       cg_meta_cols=list(cg_id="IlmnID", cg_chr="CHR", cg_coord="MAPINFO", island_name="UCSC_CpG_Islands_Name", island_region="Relation_to_UCSC_CpG_Island"),
                       include_regions=c(), sort_by_coord=FALSE){
  if(sort_by_coord){
    cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "cg_chr", "cg_coord", "island_name", "island_region"), index_column="cg_id")
    cg_meta <- cg_meta[order(cg_meta[,cg_meta_cols$cg_chr], cg_meta[,cg_meta_cols$cg_coord]),]
  }else{
    cg_meta <- .check_cgMeta_wCols(cg_meta, cg_meta_cols, necessary_cols=c("cg_id", "island_name", "island_region"), index_column="cg_id")
  }
  # CpGs which are part of an island based on UCSC annotation
  # by default, we only report CpGs inside the island, but we may add additional regions such as shores and shelves
  island_regions <- c("Island")
  island_regions <- unique(c(island_regions, intersect(include_regions, c('N_Shore','S_Shelf','S_Shore','N_Shelf'))))
  return(cg_meta[(cg_meta[,cg_meta_cols$island_name]==cg_island)&(cg_meta[,cg_meta_cols$island_region]%in%island_regions),cg_meta_cols$cg_id])
}