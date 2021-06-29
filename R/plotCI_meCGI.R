#' Beta value error plot for a CpG island
#' 
#' Description
#'
#' @details Plots error bars and confidence intervals for beta values of CpGs 
#' belonging to a given CGI. This function calls the plotCI_meCGs function 
#' where the range is given by the CGI start and end coordinates. 
#' Additionally, position of genes or exons can be plotted on the x-axis. 
#' 
#' @export
plotCI_meCGI <- function(beta_values, sample_groups,
                         cg_island=NULL, chr=NULL, include_regions=c(),
                         left_bp=0, right_bp=0,
                         add_cg=FALSE, add_genes=FALSE, add_exons=FALSE){
  
  
  cg_meta <- data("CpG_anno450K", package="meNet")
  genes_meta <- data("UCSC_genes", package="meNet")
  
  return(0)
}