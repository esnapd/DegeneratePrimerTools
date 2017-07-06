#'
#' Read HMMMer DomTbl File
#'
#' HMMER outputs several filetypes this one will read its Domain scanning program
#'
#' @param filename. Required. HMMDom table output file
#' @importFrom data.table fread
#' @export
#'
load_hmmmerdomtbl <- function(filename){
  inputstring = paste("grep -v '^#'", filename)
  hmmtbl <- fread(input = inputstring,
                  col.names  = c('target', 't_accession','tlen', 'queryname','q_accession',
                                 'qlen', 'full_e-val','full_score','full_bias','dom_number',
                                 'dom_of','dom_c-Evalue','dom_i-Evalue','dom_score','dom_bias',
                                 'hmmcoord_from','hmmcoord_to','alncoord_from','alncoord_to',
                                 'envcoord_from', 'envcoord_to','acc','target_decsription'))
}
#' Load Blast Tabular Output Files.
#'
#' Loading blast data into R and returning a \code{\link{[data.table]{data.table}}} and
#' creating an index on the index columns; keys on the QueryI and SubjectID
#'
#' @param blastfile. Required. Path location of a blastfile.
#' @param indexcols. Optional. List of columnvalues to set the index on.
#' @importFrom data.table fread
#' @importFrom data.table setkeyv
#'
#' @export
load_blast <- function(blastfile, indexcols = c("queryID")) {
  column_names <- c("queryID",  "subjectID", "percent.identity",
                    "alignment.length", "mismatches", "gap.openings", "qstart", "qend",
                    "sstart", "send", "evalue", "bitscore")
  
  # check index columns
  for (ival in indexcols) {
    if (!ival %in% column_names) {
      stop(paste("bad values in the indexcols. only valid column names can be used:", paste(column_names, collapse = " ")))
    }
  }
  
  dt <- fread(input=blastfile, header=FALSE, col.names = column_names)
  setkeyv(dt, cols = indexcols)
  
  return(dt)
}