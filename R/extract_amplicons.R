
#' Extract Amplicons from Sequences
#'
#' Find matches to degnerate primers against a nucleotide sequence and
#' return the subsetted sequence or NULL. This matching algorith is based on
#' \code{\link[Biostrings]{matchPattern}} to find matches for the degenerate
#' primer seqeunces in your reference sequecne. This role of for this function
#' is primarily within the context of extractign expected amplicon
#' sequences and we therefore expect only a single amplicon per target and
#' explicitly drop any seqeucence with greater than one match to a target sequence.
#' Therefore this function is not suitable to finding primer targets on genomic
#' DNA, and should be used for finding primer matches within the context of defined
#' reference sequences like single genes, or the DNA seqeunces corresponding to conserved
#' protein domains.
#'
#' @param fp (Required).  Default \code{NULL}. Character string of the forward primer sequence.
#' @param rp (Required).  Default \code{NULL}. Character string of the reverse primer sequence.
#' @param drop.multiple (optional). Default \code{TRUE}. Logical.
#'   If there is more than one match to the forward or reverse sequence, should the result be dropped.
#' @param max.mismatch   (optional). Default\code{2}. Maximum allowed mismatch betwen primer sequences and the target.
#' @export
setGeneric("extract_amplicons", function(object, fp, rp, drop.multiple=TRUE, max.mismatch = 2) standardGeneric("extract_amplicons"))
#'
#' @importFrom  Biostrings DNAStringSet DNAString matchPattern
setMethod("extract_amplicons", "DNAString", function(object, fp, rp, drop.multiple=TRUE, max.mismatch = 2){
  fmatches <- matchPattern(fp,     object, fixed=FALSE, max.mismatch = max.mismatch)
  rmatches <- matchPattern(rc(rp), object, fixed=FALSE, max.mismatch = max.mismatch)
  
  # Return NULL if there are No Matches
  if (length(fmatches) == 0 || length(rmatches) == 0) return(NULL)
  
  # drop the match or provide a warning if there are multiple matches
  if (length(fmatches) > 1 & drop.multiple == TRUE) return(NULL)
  if (length(fmatches) > 1)  warning("More than one match for the forward primer, using the first")
  if (length(rmatches) > 1 & drop.multiple == TRUE) return(NULL)
  if (length(rmatches) > 1)  warning("More than one match for the reverse primer, using the first")
  
  # check for negative indicies
  ampliconlength <- end(rmatches) - start(fmatches)
  if (ampliconlength <= 0) return(NULL)
  
  return(object[start(fmatches):end(rmatches)])
})
#' @importFrom  Biostrings DNAStringSet DNAString matchPattern
#' @importFrom purrr compact
setMethod("extract_amplicons", "DNAStringSet", function(object, fp, rp, drop.multiple=TRUE, max.mismatch = 2){
  amplicons <- lapply(object, function(x) {extract_amplicons(x, fp=fp, rp=rp)})
  amplicons <- compact(amplicons)
  DNAStringSet(amplicons)
})