#' Beta value error plot with CpG coordinates
#'
#' @details Plots error bars and confidence intervals for beta values of given 
#' CpGs. If list of CpGs is given, they are plotted equidistantly. 
#' If a range of coordinates and a chromosome are given, then all CpGs inside 
#' the range are plotted and positioned on x-axis based on their chromosomal 
#' coordinate. If the latter case, positions of unmeasured CpGs inside the 
#' region can be added to the plot. 
#'
#' @export
plotCI_meCpG <- function(beta_values, sample_groups,
                         cg_list=NULL, chr=NULL, first_coord=NULL,last_coord=NULL,
                         add_cg=FALSE, add_genes=FALSE, add_exons=FALSE){
  
  cg_meta <- data("CpG_anno450K", package="meNet")
  genes_meta <- data("UCSC_genes", package="meNet")
  
  return(0)
}