#' run vsearch clustering
#'
#' Function for running the vsearch dereplication command. Accepts as inputs filename and logfile
#' can optionally return a file or a \link{[data.table] data.table} of the results in blastformat.
#'
#' @param cmd. Required. Name/Location of the vsearch command (string)
#' @param filename. Required. A file location (string)
#' @param id Required. A number [0-1] (float)
#' @param OutFasta. Required. A file location (string)
#' @param UCFile. Required. A file location (string)
#' @param logfile. Required. A file location (string)
#'
#' @importFrom Biostrings writeXStringSet
#' 
#' @export
run_vsearch_cluster <- function(filename, id=0.95, OutFasta, UCFile, logfile){
  

  cluster_args <- paste0("--cluster_fast ", filename, " -id ", id, " -centroids ", OutFasta, " -uc ", UCFile, " -sizein -sizeout")
  CmdOut <- system2(command = "vsearch", args = cluster_args, stdout = TRUE, stderr = TRUE) # Execute cluster command
  cat(paste(Sys.time(), "- Cluster sequences with vsearch at ", id, "% identity.\n"), file=logfile, sep="", append=TRUE)
  cat(paste(Sys.time(), "vsearch", cluster_args, "\n"), file=logfile, sep="", append=TRUE)

  
  # return path of outputfile of vsearch cluster command. Will be deleted after the completion of the function up.
  
}
