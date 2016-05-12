#' trim a multiplesequence alignmentfile of gaps
#'
#' @param infile
#' @param outfile
#' @param minoccupancy
#' @pram refsequence
#' @param trailgap
#' @export
#'
trimAlignment <- function(infile, outfile, minoccupancy=0, refsequence = NULL, trailgap = FALSE) {
  minoccupancy <- ifelse(minoccupancy > 0, paste0("-min ", minoccupancy), " ")
  refsequence  <- ifelse(is.null(refsequence), "", paste0("-ref ", refsequence))
  trailgap     <- ifelse(trailgap==FALSE, "","-trailgap")
  degescript   <- system.file("DEGEPRIME/TrimAlignment.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl", degescript, "-i", infile, "-o", outfile, minoccupancy, refsequence, trailgap)
  print(cli)
  system(cli)
}
#' Take An MSA Object and trim it
#'
#' use trimaligment.pl to trim an MSA object
#'
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @export
trimMSA <- function(msa) {
  temp1  <- tempfile()
  temp2 <- tempfile()

  dnaSS  <- as(msa, "DNAStringSet")
  writeXStringSet(dnaSS, file=temp1)
  trimAlignment(temp1,temp2)

  trimmed <- readDNAStringSet(temp2)
  return(trimmed)
}
#' run the DegePrime.pl script
#'
#' @param alignmentfile
#' @param outfile
#' @param oligolength
#' @pram maxdegeneracy
#' @param minimumdepth
#' @param skiplength
#' @param number_iterations
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
  print(cli)
  system(cli)
}
#' Run and load degeprime and return the dataframe of results
#'
#' @importFrom parallel mclapply
#' @export
rundegeprime <- function(alignmentfile, oligolength, degeneracyrange=c(1,4,100,400,1000),
                         minimumdepth=1, skiplength=20, number_iterations=100, ncpus=1) {

  # use degeneracy range to determine the nubmer of jobs
  drange <- seq(degeneracyrange)
  tempfiles <- lapply(drange, function(x) {tempfile()})

  degendata <- mclapply(drange, function(x) {
    #get per-run data
    outputfile    <- tempfiles[[x]]
    maxdegeneracy <-  degeneracyrange[[x]]
    #calculate degeneracies
    degePrime(alignmentfile=alignmentfile, outfile=outputfile,
              oligolength=oligolength, maxdegeneracy = maxdegeneracy,
              minimumdepth=minimumdepth, skiplength=skiplength, number_iterations=number_iterations)
    #load the file and return a dataframe
    df <- loadDegePrimeDF(outputfile)
    df$degeneracy <- maxdegeneracy
    return(df)
  }, mc.cores=ncpus)

  #aggregate the data
  aggdata <- Reduce(rbind, degendata)

  #add coverage information
  aggdata$coverage <- aggdata$PrimerMatching/aggdata$TotalSeq

  return(aggdata)
}
