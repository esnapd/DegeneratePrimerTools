#' Create Paired-End Reference Sequences from DNAStringSets
#'
#' When working with Amplicons that are too long to merge together it is sometimes desireable to concatenate
#' the forward and reverse reads together. In order to classify these sequences using Blast, you would need
#' to generate a reference database of the same "shape" as the amplicons. This function takes a \code{\link[Biostrings]{DNAStringSet}},
#' and two lengths. For amplicon-reference workflows, this can be used downstream of \code{\link{extract_amplicons}} to
#' the ends of your amplicon.
#'
#' @param dnastringset (Required).  A \code{\link[Biostrings]{DNAStringSet}} containing your target sequences.
#' @param trimf (optional). Default \code{240}. Amount to be trimmed from the forward primer match
#' @param trimr (optional). Default \code{175}. Amount to be trimmed from the reverse primer match
#'
#'
#' @importFrom Biostrings xscat subseq
#' @export
#' @examples
#' extract_ends(dnatest, trimf=5, trimr=3)
extract_ends <- function(dnastringset, trimf=240, trimr=175) {
  #get forward and reverse and concatenate together
  amF <- subseq(dnastringset, end=trimf)
  amR <- subseq(dnastringset, start=-trimr)
  ampliconends <- xscat(amF, amR)
  names(ampliconends) <- dnastringset
  
  #filter minimum size
  ampliconends[width(ampliconends) >= trimf+trimr]
}