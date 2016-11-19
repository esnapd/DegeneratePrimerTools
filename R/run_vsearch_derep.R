#' run vsearch global
#'
#' Convenience function for running the vsearch global search. Can accept as inputs either filenames or DNAStringSets and
#' can optionally return a file or a \link{[data.table] data.table} of the results in blastformat.
#'
#' @param query. Required. Either a file location  (string) or a DNAStringSet
#' @param target. Required. Either a file location (string) or a DNAStringSet
#' @param outfile. Optional. Default \code{NULL} Either NUll or a filelocation
#' @param id. Optional. Numeric value between 0 and 1. Default \code{0.95} Percent Identity of targets to subjects
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom data.table fread
#' @export
run_vsearch_global <- function(query, target, outfile=NULL,id=0.95) {
  #intermediate files will be run in temp if no output is sent, andl lets cleaup the tempdir when we are done.
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  # handle the query. if a character, check the file exists
  # if its an DNAString set, write the fastafile out
  if (is.character(query)) {
    query <- file.path(normalizePath(dirname(query)), basename(query))
    if(length(Sys.glob(paste(query, "*", sep="")))<1) stop("query sequence does not exit!")
    queryfile <- query
    
  } else if (is(query, "DNAStringSet")) {
    if(is.null(names(target))) stop("DNAStringSets must have names")
    queryfile <- tempfile(tmpdir = wd, fileext="blastquery")
    writeXStringSet(query, queryfile)
    
  } else {
    stop("query must be a file location or a DNAStringSet")
  }
  
  # handle the target. if a file, make sure that blastcmddb has been run.
  # if its an DNAString set, write a blastdb in the temp directory
  if (is.character(target)) {
    targetdb <-  file.path(normalizePath(dirname(target)), basename(target))
    if(length(Sys.glob(paste(targetdb, "*", sep="")))<1) stop("target sequence file does not exit!")
    
  } else if (is(target, "DNAStringSet")) {
    if(is.null(names(target))) stop("DNAStringSets must have names")
    
    targetdb <- tempfile(tmpdir = wd, fileext = "targetdb")
    writeXStringSet(x=target, filepath = targetdb)
    
  } else {
    stop("target must be a file location or a DNAStringSet")
  }
  
  # create the blast outputfile here
  blastout <- tempfile(tmpdir=wd, fileext = "blastout")
  if (!is.null(outfile)) {
    blastout <- outfile
  }
  
  #run the blast job
  vsearchcmd <- paste("vsearch --usearch_global", queryfile, "-db", targetdb, "-id", id, '--blast6out', blastout)
  print(vsearchcmd)
  system(vsearchcmd)
  
  
  # return. If the fileoutput has been specified, the blast value has been saved there. Otherwise
  # read in the blast table.
  if (is.null(outfile)) {
    try(blastdt <- load_blast(blastout), silent=TRUE)
    if (!exists("blastdt")) stop("BLAST did not return a match!")
    
    return(blastdt)
  }
}
