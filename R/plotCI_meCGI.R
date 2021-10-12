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
#' @param cg_island Name of a CpG island. It should be given as a UCSC CpG island name,
#' stating chromosome, the first and the last coordinate (e.g. `"chrX:1-2"`).
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
#' of the given CpG island (from 5' end to 3' end). 
#' Additionally, CpG island shores and CpG island shelves may be included in the 
#' plot. CpG island shores are genomic regions spanning up to 2 kb upstream/downstream 
#' from the CpG island. Similarly, CpG island shelves are genomic regions located 
#' from 2 to 4 kb  upstream/downstream from the CpG island \insertCite{bibikova2011high}{meNet}.
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
                         cg_island=NULL, include_regions=NULL,
                         upstream_bp=0, downstream_bp=0, ...){
  
  data("UCSC_CGI", package="meNet")
  
  # check validity of given parameters
  ####################################
  # beta values
  if(!is.data.frame(beta_values)){
    stop('"beta_values" has to be a data frame.')
  }
  if(nrow(beta_values)==0|ncol(beta_values)==0){
    stop('"beta_values" missing data.')
  }
  if(!all(rownames(beta_values)%in%CpG_annoHM450K[,"IlmnID"])){
    stop('Row names of "beta_values" are not CpGs measured by Illumina HM450K microarray.')
  }
  # sample groups
  if(!is.data.frame(sample_groups)){
    stop('"sample_groups" has to be a data frame.')
  }else if(ncol(sample_groups)<2 ){
    stop('"sample_groups" has to have two columns, the first with sample IDs and the second with the corresponding group.')
  }else if(!all(colnames(beta_values)%in%sample_groups[,1])){
    stop('First column of "sample_groups" has to contain all column names of "beta_values"')
  }
  # included regions
  if(!is.null(include_regions)){
    if(!is.character(include_regions)){
      warning('If provided, "include_regions" has to be a vector of strings. Given value is discarded.')
      include_regions <- NULL
    }else{
      island_regions <- intersect(include_regions, c('N_Shore','S_Shelf','S_Shore','N_Shelf'))
      if(length(island_regions)<length(include_regions)){
        warning('Check the value of parameter "include_regions".')
      }
      if(length(island_regions)==0){
        island_regions <- NULL
      }
    }
  }
  # additional bps
  if(!is.numeric(upstream_bp)){
    stop('"upstream_bp" has to have a numeric value.')
  }else if(upstream_bp<0){
    stop('"upstream_bp" has to have a non-negative value.')
  }
  #
  if(!is.numeric(downstream_bp)){
    stop('"downstream_bp" has to have a numeric value.')
  }else if(downstream_bp<0){
    stop('"downstream_bp" has to have a non-negative value.')
  }
  
  # cg island
  ###########
  # remove "chr" if string starts with it
  cg_island <- sub("chr", "", cg_island)
  # check if an island is an UCSC island
  if(!(cg_island%in%UCSC_CGI$UCSC_CpG_Islands_Name)){
    warning("Given CpG island is not an UCSC CpG island.")
  }
  
  # calculate plotting values
  ###########################
  # find start and end of the island
  chr <- strsplit(cg_island, ":")[[1]][1]
  start_position <- strsplit(strsplit(cg_island, ":")[[1]][2], "-")[[1]][1]
  end_position <- strsplit(strsplit(cg_island, ":")[[1]][2], "-")[[1]][2]
  # check if the CpG island name is correctly written
  error_message <- "CpG island name is not in the correct format."
  tryCatch(
    {
      chr <- as.numeric(chr)
      start_position <- as.numeric(start_position)
      end_position <- as.numeric(end_position)
    }, 
    error=function(e){stop(error_message)}, 
    warning=function(w){stop(error_message)}
    )
  
  # adding base pairs to plot
  if("N_Shelf"%in%include_regions){
    start_position <- start_position - 4000
  }else if("N_Shore"%in%include_regions){
    start_position <- start_position - 2000
  }
  if("S_Shelf"%in%include_regions){
    end_position <- end_position + 4000
  }else if("S_Shore"%in%include_regions){
    end_position <- end_position + 2000
  }
  start_position <- start_position - upstream_bp
  end_position <- end_position + downstream_bp
  
  # PLOTTING
  plotCI_meCpG(beta_values=beta_values, sample_groups=sample_groups,
               chr=chr, first_coord=start_position+1, last_coord=end_position,...)

}
