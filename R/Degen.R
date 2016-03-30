#' run the Degen script on a character vector
#'
#' @param infile
#' @param method
#' @export
#'
makedegenerate <- function(infile, method="S") {
  if (!method %in% c("S","Z","SZ")) stop("Degeneracy method is only one one of 'S','Z'. or 'SZ'")
  table        <- paste0("-t=","\'",method,"\'")
  degescript   <- system.file("Degen/Degen_v1_4.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl", degescript, infile, table)
  print(cli)
  system(cli)
}

paste0("'","S","'")
