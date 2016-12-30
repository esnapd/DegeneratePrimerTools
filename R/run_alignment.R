#' Run a Multiple Sequence Alignment
#'
#' @param dgprimer Required. A degeprimer object.
#' @param program Optional. Default \code{"msa"}. This option allows oyu to specify which preogrma ot use. By default
#' the \code{\link[msa]{msa}} program is used. Alternatively, "muscle" can be specified to use the CLI version of 
#' muscle
#' @param method Optional. Default \code{"ClustalW"}. The choice of alignment software to pass to \code{\link[msa]{msa}}.
#'            Can be  "ClustalW", "ClustalOmega", or  "Muscle".
#' @param maxiters Optional. Default \code{"default"}
#' @param force Optional. Default is \code{FALSE}. If an msa already exists, you must set FORCE to TRUE to overwrite it.
#' @importFrom msa msa
#' @export
run_alignment <- function(dgprimer, program="msa", method="ClustalW", maxiters="default", force=FALSE,...) {
  
  if (!program %in% c("msa", "muscle")) {
    stop("The program option must be either 'msa' or 'muscle' ")
  }
  
  if (force==FALSE  & !is.null(dgprimer@msa)) {
    stop("Your degeprime object already has an MSA. To overwrite use force=TRUE")
  }
  
  seqs <- dgprimer@refseq
  
  if (program == "msa") {
    aln <- msa::msa(seqs, method=method, maxiters=maxiters, type="dna", ...)
    aln <- as(aln, "DNAMultipleAlignment")
  }
  if (program == "muscle") {
    aln <- run_muscle(seqs) 
    aln <- as(aln, "DNAMultipleAlignment")
  }
    
  dgprimer@msa <- aln
  
  return(dgprimer)
}