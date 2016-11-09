#' Run a Multiple Sequence Alignment
#'
#' @param dgprimer. Required. A degeprimer object.
#' @param method. Optional. Default \code{"ClustalW"}. The choice of alignment software to pass to \code{\link[msa] msa}.
#'            Can be  "ClustalW", "ClustalOmega", or  "Muscle".
#' @param force. Optional. Default is \code{FALSE}. If an msa already exists, you must set FORCE to TRUE to overwrite it.
#' @importFrom msa msa
#' @export
run_alignment <- function(dgprimer, method="ClustalW", maxiters="default", force=FALSE,...) {
  
  if (force==FALSE  & !is.null(dgprimer@msa)) {
    stop("Your degeprime object already has an MSA. To overwrite use force=TRUE")
  }
  
  seqs <- dgprimer@refseq
  aln  <- msa::msa(seqs, method=method, maxiters=maxiters, type="dna", ...)
  
  dgprimer@msa <- as(aln, "DNAMultipleAlignment")
  return(dgprimer)
}