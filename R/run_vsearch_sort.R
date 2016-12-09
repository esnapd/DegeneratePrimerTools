#' run vsearch sort
#'
#' Function for running the vsearch dereplication command. Accepts as inputs filename and logfile
#' can optionally return a file or a \link{[data.table] data.table} of the results in blastformat.
#'
#' @param filename. Required. A file location (string)
#' @param SortedFile. Required. A file location (string)
#' @param MinSize. Required. Minimum size of clusters to keep (integer)
#' @param logfile. Required. A file location (string)
#'
#' @importFrom Biostrings writeXStringSet
#' 
#' @export
run_vsearch_sort <- function(filename, SortedFile, MinSize, logfile){
  
  
  #run the vsearch dereplication command.
  sort_args <- paste0("--sortbysize ", filename, " -output ", SortedFile, " -minsize ", MinSize)
  CmdOut <- system2(command = "vsearch", args = sort_args, stdout = TRUE, stderr = TRUE) # Execute sort command
  
  sampleName <- paste(unlist(strsplit(basename(filename), "_"))[1:2], collapse="_")
  cat(paste(Sys.time(), "- Sort sequences with Vsearch minimum size ", MinSize, " for sample", sampleName, "\n"), file=logfile, sep="", append=TRUE)
  cat(paste(Sys.time(), "vsearch", sort_args, "\n"), file=logfile, sep="", append=TRUE)
  
  
  # return path of outputfile of vsearch sort command. Will be deleted after the completion of the function up.
  
}
