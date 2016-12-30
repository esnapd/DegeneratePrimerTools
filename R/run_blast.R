#' run blast
#'
#' Convenience function for running blastn. Can accept as inputs either files or DNAStringSets and
#' can optionally return a file or a \link{[data.table] data.table} of the blast results.
#'
#' @param query Required. Either a file location  (string) or a DNAStringSet
#' @param target Required. Either a file location (string) or a DNAStringSet
#' @param outfile Optional. Default \code{NULL} Either NUll or a filelocation
#' @param parallel Boolean. Default \code{FALSE}. Whether or not to use GNU parallel for speed.
#' @param blast_args Optional. Default \code{""}. String arguments to pass to the blast commandline.
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom data.table fread
#' @export
run_blast <- function(query, target, outfile=NULL, parallel=FALSE, blast_args="") {
  #intermediate files will be run in temp if no output is sent.
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  
  # handle the query. if a file, make sure that blastcmddb has been run.
  # if its an DNAString set, write a blastdb in the temp directory
  if (is.character(query)) {
    query <- file.path(normalizePath(dirname(query)), basename(query))
    if(length(Sys.glob(paste(query, "*", sep="")))<1) stop("query sequence does not exit!")
    queryfile <- query
    
  } else if (is(query, "DNAStringSet")) {
    if(is.null(names(query))) stop("DNAStringSets must have names")
    queryfile <- tempfile(tmpdir = wd, fileext="blastquery")
    writeXStringSet(query, queryfile)
    
  } else {
    stop("query must be a file location or a DNAStringSet")
  }
  
  # handle the target. if a file, make sure that blastdbcmd has been run.
  # if its an DNAString set, write a blastdb in the temp directory
  if (is.character(target)) {
    targetdb <-  file.path(normalizePath(dirname(target)), basename(target))
    if(length(Sys.glob(paste(targetdb, "*", sep="")))<1) stop("target sequence file does not exit!")
    
    blastdbinfo <- system2("blastdbcmd", paste("-db", targetdb, "-info"), stdout=TRUE, stderr = TRUE) #system(paste("blastdbcmd -db", targetdb, "-info"), intern=TRUE)
    if (grepl("No alias or index file found", blastdbinfo)) {
      warning("No blast info found for this file. Creating blastdb")
      system(paste("makeblastdb -dbtype nucl -in", targetdb," -parse_seqids"))
    } else{
      cat(paste(blastdbinfo, collapse="\n"))
      cat("\n")
    }
    
  } else if (is(target, "DNAStringSet")) {
    if(is.null(names(target))) stop("DNAStringSets must have names")
    
    targetdb <- tempfile(tmpdir = wd, fileext = "targetdb")
    writeXStringSet(x=target, filepath = targetdb)
    system(paste("makeblastdb -dbtype nucl -in", targetdb))
    
  } else {
    stop("target must be a file location or a DNAStringSet")
  }
  
  # create the blast outputfile here
  blastout <- tempfile(tmpdir=wd, fileext = "blastout")
  if (!is.null(outfile)) {
    blastout <- outfile
  }
  
  #run the blast job
  if (parallel) {
    #run the query using GNU parallel. See https://www.biostars.org/p/63816/
    blastcmd <- paste("cat", queryfile, " | parallel --block 100k --recstart '>' --pipe",
                      "blastn", "-query - -db", targetdb, "-outfmt 6", blast_args, " > ", blastout)
  } else {
    blastcmd <- paste("blastn", "-db", targetdb, "-query", queryfile, "-out", blastout, '-outfmt 6', blast_args)
  }
  message(blastcmd)
  system(blastcmd)
  
  
  # return. If the fileoutput has been specified, the blast value has been saved there. Otherwise
  # read in the blast table.
  if (is.null(outfile)) {
    try(blastdt <- load_blast(blastout), silent=TRUE)
    if (!exists("blastdt")) stop("BLAST did not return a match!")
    
    return(blastdt)
  }
}
