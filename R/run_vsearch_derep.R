#' run vsearch dereplication
#'
#' Function for running the vsearch dereplication command. Accepts as inputs filename and logfile
#' can optionally return a file or a \link{[data.table] data.table} of the results in blastformat.
#'
#' @param cmd. Required. Name/Location of the vsearch command (string)
#' @param filename. Required. A file location (string)
#' @param DerepFile Required. A file location (string)
#' @param logfile. Required. A file location (string)
#'
#' @importFrom Biostrings writeXStringSet
#' 
#' @export
run_vsearch_derep <- function(cmd, filename, DerepFile, logfile){

  
  #run the vsearch dereplication command.
  derep_args <- paste0("--derep_fulllength ", filename, " -output ", DerepFile, " -sizeout")
  CmdOut <- system2(command = cmd, args = derep_args, stdout = TRUE, stderr = TRUE) # Execute dereplicate command
  
  cat(paste(Sys.time(), "- Dereplication with vsearch\n"), file=logfile, sep="", append=T)
  cat(paste(Sys.time(), "vsearch", derep_args, "\n"), file=logfile, sep="", append=T)
  
  
  ###################################################################################################################################
  # Collect information from vsearch derep command output for number of sequences, number of unique sequences and vsearch version.
  
  sequ <- as.numeric(sub(".*in (.*) seqs.*", "\\1", CmdOut[grep("nt in", CmdOut)]))
  derep <- as.numeric(sub("(.*)unique sequences.*", "\\1", CmdOut[grep("unique sequences", CmdOut)]))
  version <- sub("(.*)unique sequences.*", "\\1", CmdOut[grep("vsearch", CmdOut)][1])
  
  message(paste0(sequ, " input sequences -> ", derep, " dereplicated sequences"))
  cat(paste(Sys.time(), sampleName, ": ", sequ, "input sequences ->", derep, " dereplicated sequences\n"), file=logfile, sep="", append=T)
  
  # return path of outputfile of vsearch derep command. Will be deleted after the completion of the function up.

}
