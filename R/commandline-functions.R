#' run the DegePrime.pl script
#'
#' DEGEPRIME is a command-line perl program to find degenerate primers from a multiple sequence alignment. It attemps
#' to find the maximumcoverage ofthe MSA with the minimum degeneracy where the maximum degeneracy is bounded by a user-provided input/
#' This is a comman-line wrapper that shells ou tht ework to DEGEPRIME.
#'
#' @param alignmentfile Required. Input alignment path.
#' @param outfile Required. Output primer table path.
#' @param oligolength Required. An integer of how long the primers to search for should be.
#' @param maxdegeneracy Required. An integer specifying the maximum degeneracy.
#' @param minimumdepth Optional. Default \code{1}.
#' @param skiplength Optional. Default \code{20}.
#' @param number_iterations Optional. Default \code{100}.
#' @export
degePrime <- function(alignmentfile, outfile,  oligolength, maxdegeneracy, minimumdepth=1, skiplength=20, number_iterations=100) {
  if (!is.numeric(oligolength)) stop("Oligolength must be a number")
  if (!is.numeric(maxdegeneracy)) stop("maxdegenderacy must be a number")

  degescript   <- system.file("DEGEPRIME/DegePrime.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl",   degescript,
                        "-i",     alignmentfile,
                        "-o",     outfile,
                        "-l",     oligolength,
                        "-d",     maxdegeneracy,
                        "-depth", minimumdepth,
                        "-skip",  skiplength,
                        "-iter",  number_iterations)
  system(cli)
}
#' trim a multiplesequence alignmentfile of gaps
#'
#' DEGEPRIME recommends trimming alignments using thir trimming script
#' which will elimiante gaps from MSAs but note their absence with lowercase
#' letters. this is a simple wrapper function for the TrimAlignment.pl
#'
#' @param infile Required. Inputfile.
#' @param outfile Required. Outputfile.
#' @param minoccupancy Optional. Default \code{0}.
#' @param refsequence Optional. Default \code{NULL}.
#' @param trailgap Optional. Default \code{FALSE}.
#' @export
#'
trimAlignment <- function(infile, outfile, minoccupancy=0, refsequence = NULL, trailgap = FALSE) {
  minoccupancy <- ifelse(minoccupancy > 0, paste0("-min ", minoccupancy), " ")
  refsequence  <- ifelse(is.null(refsequence), "", paste0("-ref ", refsequence))
  trailgap     <- ifelse(trailgap==FALSE, "","-trailgap")
  degescript   <- system.file("DEGEPRIME/TrimAlignment.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl", degescript, "-i", infile, "-o", outfile, minoccupancy, refsequence, trailgap)
  system(cli)
}
