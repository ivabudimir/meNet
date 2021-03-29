#'@export
meMultiplex_to_multinet <- function(input_file){
  multinet_object <- multinet::read_ml(input_file)
  return(multinet_object)
}
