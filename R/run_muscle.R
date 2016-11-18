#' run muscle
#'
#' Convenience function for running \href{http://drive5.com/muscle/}{muscle}. Can accept as inputs either file names or DNAStringSets and
#' can optionally return a file or a \link{[Biostrings] DNAStringSet} of the blast results.
#'
#' @param seqs Required. Either a file location  (string) or a DNAStringSet
#' @param outfile. Optional. Default \code{NULL} Either NUll or a filelocation
#' @param muscle_args. Optional. Default \code{"--maxiters=2"}. Commandline arguments to pass to FastTree.
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readDNAStringSet
#' @export
run_muscle <- function(seqs,  outfile=NULL, muscle_args="-maxiters 2") {
  #intermediate files will be run in temp if no output is sent.
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  # handle the seqs if a file, make sure that FastTree has been run.
  # if its a DNAString set, write it to a file.
  if (is.character(seqs)) {
    seqs <- file.path(normalizePath(dirname(seqs)), basename(seqs))
    if(length(Sys.glob(paste(seqs, "*", sep="")))<1) stop("query sequence does not exit!")
    queryfile <- seqs
    
  } else if (is(seqs, "DNAStringSet")) {
    if(is.null(names(seqs))) stop("DNAStringSets must have names")
    queryfile <- tempfile(tmpdir = wd, fileext = "inputseqs")
    writeXStringSet(seqs, queryfile)
    
  } else {
    stop("seqs must be a file location or a DNAStringSet")
  }
  
  # create the blast outputfile here
  alnout <- tempfile(tmpdir = wd, fileext = "alnout")
  if (!is.null(outfile)) {
    alnout <- outfile
  }
  
  # run FastTree
  musclecmd <- paste("muscle", "-in", queryfile, "-out", alnout, muscle_args )
  system(musclecmd)
  
  
  # return. If the fileoutput has been specified, the blast value has been saved there. Otherwise
  # read in the blast table.
  if (is.null(outfile)) {
    try(aln <- readDNAStringSet(alnout), silent=TRUE)
    if (!exists("aln")) stop("Unable to open the generated alignment file!")
    
    return(aln)
  }
}
