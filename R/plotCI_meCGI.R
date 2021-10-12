#' Beta value error plot for a CpG island
#' 
#' @description Plots methylation beta values for the CpGs within the given CpG
#' island in the form of error bars. For every group of samples, the mean 
#' methylation value is plotted together with the confidence interval. 
#' X-axis position of CpG sites represents their true genomic coordinates. 
#' Additional information can be included in the plot such as position of genes
#' and exons.
#' 
#' @param beta_values Data frame with methylation beta values. CpG sites should
#' be listed in the rows with the row names being the Illumina CpG identifiers.
#' Samples should be listed in the columns.
#' @param sample_groups Data frame which defines the sample groups. It has to 
#' have two columns, the first column with sample IDs and the second column 
#' with the corresponding group.
#' @cg_island Name of a CpG island.
#' @param include_regions Regions of CGIs which should be included in the plot 
#' besides the island itself. Can take any subset of values `"S_Shelf"`, `"S_Shore"`, 
#' `"N_Shore"`, `"N_Shelf"`. By default, only the CpGs inside the island are reported.
#' @param upstream_bp Number of additional base pairs added to the plot to the
#' left. Default value is `0`.
#' @param downstream_bp Number of additional base pairs added to the plot to the
#' right. Default value is `0`.
#' @param ... Additional plotting parameters passed to the function 
#' `meNet::plotCI_meCpG`.
#' 
#' @return Vector with names of plotted CpG sites, in the same order as in the
#' plot.
#'
#' @details Plots error bars and confidence intervals for beta values of CpGs 
#' belonging to a given CGI. This function calls the `meNet::plotCI_meCpG` 
#' function where the plotted range is given as the start and the end coordinate
#' of the given CpG island. Additionally, CpG island shores and CpG island
#' shelves may be included in the plot. CpG island shores are genomic regions
#' spanning up to 2 kb upstream/downstream from the CpG island. Similarly,
#' CpG island shelves are genomic regions located from 2 to 4 kb 
#' upstream/downstream from the CpG island \insertCite{bibikova2011high}{meNet}.
#' If a shelf is included in the plot, also the shore which separates it from 
#' the CpG island is included.
#' Furthermore, the plot can be expanded upstream or downstream with parameters
#' `upstream_bp` and `downstream_bp`.
#' 
#' More information can be included on the x-axis, such as the exact position
#' of the CpG island, position of genes and exons or positions of unmeasured CpG
#' sites. For further details, see documentation of the function
#' `meNet::plotCI_meCpG`.
#' 
#' @references
#'       \insertAllCited{}
#' 
#' @export
plotCI_meCGI <- function(beta_values, sample_groups,
                         cg_island=NULL, include_regions=c(),
                         upstream_bp=0, downstream_bp=0, ...){
  
  
  cg_meta <- data("CpG_anno450K", package="meNet")
  genes_meta <- data("UCSC_genes", package="meNet")
  
  island_regions <- c("Island")
  island_regions <- unique(c(island_regions, intersect(include_regions, c('N_Shore','S_Shelf','S_Shore','N_Shelf'))))
  
  return(0)
}