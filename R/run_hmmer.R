#' run hmmer
#'
#' Convenience function for running \href{http://hmmer.org/}{hmmer3}. Can accept as inputs either file names or DNAStringSets and
#' can optionally return a file or a \link{[data.table] data.table} of the results.
#'
#' @param seqs Required. Either a file or an \link{[Biostrings] AAStringSet} to be probed by the PFAMDB
#' @param hmmerdb Required. Location of the hmmer3 binary file.
#' @param outfile. Optional. Default \code{NULL} Either NUll or a filelocation.
#' @param hmmer_args. Optional. Default \code{"--maxiters=2"}. Commandline arguments to pass to FastTree.
#'
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readAAStringSet
#' @export
run_hmmer <- function(seqs, hmmerdb, outfile=NULL, hmmer_args="--domE 1e-30 --cpu 1") {
  #intermediate files will be run in temp if no output is sent.
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  # handle the seqs if a file, make sure that FastTree has been run.
  # if its a DNAString set, write it to a file.
  if (is.character(seqs)) {
    seqs <- file.path(normalizePath(dirname(seqs)), basename(seqs))
    if(length(Sys.glob(paste(seqs, "*", sep="")))<1) stop("query sequence does not exit!")
    queryfile <- seqs
    
  } else if (is(seqs, "AAStringSet")) {
    if(is.null(names(seqs))) stop("AAStringSets must have names")
    queryfile <- tempfile(tmpdir = wd, fileext = "inputseqs")
    writeXStringSet(seqs, queryfile)
    
  } else {
    stop("seqs must be a file location or a AAStringSet")
  }
  
  # create the hmmm outputfile here
  hmmout <- tempfile(tmpdir = wd, fileext = "hmmout")
  if (!is.null(outfile)) {
    hmmout <- outfile
  }
  
  # run Hmmer3
  hmmercmd <- paste("hmmsearch  --domtblout", hmmout, hmmer_args, "-out", hmmerdb, queryfile)
  print(hmmercmd)
  system(hmmercmd)
  
  
  # return. If the fileoutput has been specified, the blast value has been saved there. Otherwise
  # read in the blast table.
  if (is.null(outfile)) {
    try(hmm <- load_hmmmerdomtbl(hmmout), silent=TRUE)
    if (!exists("hmm")) stop("Unable to open the generated alignment file!")
    
    return(hmm)
  }
}
