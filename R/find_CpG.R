#' CpG coordinates on a chromosome
#'
#' @description For a given chromosome, returns positions of all "CG" 
#' dinucleotides (CpG sites) based on a human genome assembly hg19. If minimum
#' or maximum position is given, only CpGs in the given range are returned.
#' 
#' @param chromosome Name of a human chromosome in the format `"chr?"`, where `"?"`
#' must either be a number from `1` to `22` or `"X"` or `"Y"`.
#' @param min_coord Minimum returned coordinate.
#' @param max_coord Maximum returned coordinate.
#' 
#' @return Vector of CpG coordinates.
#' 
#' @details Returned position is a position of "C". The same convention is used
#' by Illumina Infinium HM450K manifest file.
#' Function is based on R package `"BSgenome.Hsapiens.UCSC.hg19"`
#' \insertCite{team2020bsgenome}{meNet}.
#' 
#' @references
#'       \insertAllCited{}
#' 
#' @import BSgenome.Hsapiens.UCSC.hg19
#'
#' @export
find_CpG <- function(chromosome, min_coord=NULL, max_coord=NULL){
  all_chromosomes <- c(paste0("chr",1:22), "chrX", "chrY")
  
  if(!is.character(chromosome)){
    stop('Chromosome must be given in the "chr?" format.')
  }else if(!(chromosome %in% all_chromosomes)){
    stop('Chromosome must be given in the "chr?" format, where "?" is either number from 1 to 22 or "X" or "Y".')
  }
  if(!is.null(min_coord) & !is.numeric(min_coord)){
    stop('"min_coord" must be a number.')
  }
  if(!is.null(max_coord) & !is.numeric(max_coord)){
    stop('"max_coord" must be a number.')
  }
  
  coordinates <- start(matchPattern('CG', Hsapiens[[chromosome]]))
  
  if(!is.null(min_coord)){
    coordinates <- coordinates[coordinates >= min_coord]
  }  
  if(!is.null(max_coord)){
    coordinates <- coordinates[coordinates <= max_coord]
  }  
  
  return(coordinates)
  
}
