#' Multiplex: `meNet` file to `multinet` object
#' 
#' @description Converts file with multiplex structure written in `multinet`-style
#' to `multinet` object. The file can be created with `meNet::meMultiplex`
#' function with the output type set to `"multinet"`. This allows the use
#' of different `multinet` function for multilayer networks 
#' \insertCite{multinet}{meNet}.
#' 
#' @param input_file File with the `multinet` structure of the multiplex.
#' 
#' @return Multiplex network as `multinet` object.
#' 
#' @references
#'       \insertAllCited{}
#'
#' @export
meMultiplex_to_multinet <- function(input_file){
  multinet_object <- multinet::read_ml(input_file)
  return(multinet_object)
}
