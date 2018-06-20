87645879#' Design Primers
#'
#' A convenience function to find degenerate primers using the
#' seqeuences of the degenerateprimers object.  Return teh resutls to
#' the primerdata slot.
#'
#' @param dgprimer Required. A degeprimer object.
#' @param oligolength Optional. Default is \code{21}. The length of primers to design.
#' @param maxdegeneracies Optional. Default is \code{c(1,4,100,400,1000)}. 
#' @param minimumdepth Optional. Default is \code{1}. 
#' @param skiplength Optional. Default is \code{20}. 
#' @param number_iterations Optional. Default is \code{100}. The nubmer of
#'  iterations of the DEGEPRIMER algorithm.
#' @param ncpus Optional. efault is \code{1}. The number of cpus to use.
#' @param force Optional. Default is \code{FALSE}. Overwrite if exists?
#' 
#' @importFrom parallel mclapply
#' @export
design_primers <- function(dgprimer, oligolengths=21, maxdegeneracies=c(1,4,100,400,1000),
                           minimumdepth=1, skiplength=20, number_iterations=100, ncpus=1,
                           force=FALSE) {
  
  # check if degeprime data ia associated with the class. if so, require force=TRUE
  if (!is.null(dgprimer@primerdata) & force==FALSE) {
    stop("Degenerate primer data already exists. To overwrite please use force=TRUE")
  }
  
  # check if mutliple sequence alignment data ia associated with the class. if so, require force=TRUE
  if (is.null(dgprimer@msa)) {
    stop("Primer Calculation requires a multiple sequence alignment")
  }
  
  # use degeneracy range to determine the nubmer of jobs
  if (length(maxdegeneracies) == 1 & is.numeric(maxdegeneracies)) {
    maxdegeneracies <- c(maxdegeneracies)
  }
  
  drange    <- seq(maxdegeneracies)
  tempfiles  <- lapply(drange, function(x) {tempfile()})
  
  #create the trimmed file for DEGEPRIMER input
  #write alignments to disk and trim the alignment
  alignfile   <- tempfile()
  trimmedfile <- tempfile()
  writeXStringSet( as(dgprimer@msa,"DNAStringSet"), alignfile)
  trimAlignment(infile=alignfile, outfile=trimmedfile, minoccupancy=1, refsequence = NULL, trailgap = FALSE)
  
  degendata <- mclapply(drange, function(x) {
    #get per-run data
    outputfile    <- tempfiles[[x]]
    maxdegeneracy <- maxdegeneracies[[x]]
    #calculate degeneracies
    degePrime(alignmentfile=trimmedfile, outfile=outputfile,
              oligolength=oligolengths, maxdegeneracy = maxdegeneracy,
              minimumdepth=minimumdepth, skiplength=skiplength, number_iterations=number_iterations)
    #load the file and return a dataframe
    df <- read.table(outputfile,header = TRUE, stringsAsFactors = FALSE)
    df$degeneracy <- maxdegeneracy
    #remove intermediate file
    rm(outputfile)
    return(df)
  }, mc.cores=ncpus)
  
  #aggregate the data
  aggdata <- Reduce(rbind, degendata)
  
  #add coverage information
  aggdata$coverage <- aggdata$PrimerMatching/aggdata$TotalSeq
  
  # add the primer data to the object
  dgprimer@primerdata <- new("primerdata", aggdata)
  
  return(dgprimer)
}