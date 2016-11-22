#' run FastTree
#'
#' Convenience function for running \href{http://www.microbesonline.org/fasttree/}{FastTree}. Can accept as inputs either file names or DNAStringSets and
#' can optionally return a file or a \link{[data.table] data.table} of the blast results.
#'
#' @param seqs Required. Either a file location  (string) or a DNAStringSet
#' @param outfile. Optional. Default \code{NULL} Either NUll or a filelocation
#' @param fasttree_args. Optional. Default \code{"--maxiters=2"}. Commandline arguments to pass to FastTree.
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom ape read.tree
#' @export
run_fasttree <- function(seqs,  outfile=NULL, fasttree_args="") {
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
  treeout <- tempfile(tmpdir = wd, fileext = "treeout")
  if (!is.null(outfile)) {
    treeout <- outfile
  }
  
  # run FastTree
  fasttreecmd <- paste("FastTree", "-nt", queryfile, fasttree_args, ">", treeout )
  system(fasttreecmd)
  
  
  # return. If the fileoutput has been specified, the blast value has been saved there. Otherwise
  # read in the blast table.
  if (is.null(outfile)) {
    try(tree <- read.tree(treeout), silent=TRUE)
    if (!exists("tree")) stop("Unable to open the generated Tree file!")
    
    return(tree)
  }
}
