#' Finds CpG sites which are associated with a given CpG island
#'
#' @description For a given CpG island (CGI), function returns the list of CpGs
#' within the island. Island shores and shelves can be included in the search and
#' returned list can be sorted by the CpG coordinate. If working with Illumina Infinium
#' HM450K, 'cg_meta' and 'cg_meta_cols' parameters should keep their default values.
#' 
#' @param cg_island Name of a CpG island.
#' @param cg_meta Data frame which defines relationship between CpG sites and 
#' CpG islands. It has CpGs in the rows and their description in the columns.
#' If CpG is part of an CGI, "cg_meta" gives the name of the island and the region
#' of an island in which the CpG is found.
#' By default, function uses "CpG_anno450K" data frame which contains Illumina 
#' Infinium HumanMethylation450 manifest file.
#' @param cg_meta_cols Named list with 'cg_meta' column names. The list must include:
#' 'cg_id' naming the column with unique CpG names, 'island_name' naming the column 
#' with CGI names and 'island_region' naming the column with the CGI region in which
#' the CpG is located.
#' Default value shouldn't be changed if 'cg_meta' keeps it default value.
#' @param include_regions Regions of CGI which should be searched for CpGs besides the 
#' island itself. Can take any subset of values "S_Shelf", "S_Shore", 
#' "N_Shore", "N_Shelf". By default, only the CpGs inside the island are reported.
#' @param sort_by_coord If TRUE, returned list of CpGs will be ordered by their
#' coordinates. In this case, 'cg_meta' must contain two additional columns and
#' 'cg_meta_cols' must include two additional elements: 'cg_chr' naming column with
#' chromosome and 'cg_coord' naming column with coordinate of CpG. Defaults to FALSE.
#' 
#' @return A vector with CpG sites located within the given CpG island.
#' 
#' @export
CpG_in_CGI <- function(cg_island, cg_meta=meNet::CpG_anno450K,
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