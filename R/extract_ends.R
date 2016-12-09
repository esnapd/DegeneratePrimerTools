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
#' @param sep (optional).   Default \code{"N"}. A sequence used as a delimiter between the extracted sequences.
#'
#' @importFrom Biostrings subseq
#' @importFrom Biostrings width
#' @importFrom Biostrings xscat 
#' @importFrom Biostrings DNAString
#' @export
#' @examples
#' dnafile <- system.file("sequences","AHBA.fna",package="DegeneratePrimerTools")
#' dnatest <- Biostrings::readDNAStringSet(dnafile)
#' extract_ends(dnatest, trimf=5, trimr=3, sep="N")
extract_ends <- function(dnastringset, trimf=240, trimr=175, sep="N") {
  sep <- DNAString(sep)
  
  #get forward and reverse and concatenate together
  amF <- subseq(dnastringset, end=trimf)
  amR <- subseq(dnastringset, start=-trimr)
<<<<<<< f9628d14de2cbd4f3bf03c59f349033c195ea310
  ampliconends <- xscat(amF,sep, amR)
=======
  ampliconends <- xscat(amF, amR)
>>>>>>> addition of Minimum size parameter instead of hard code in the vsearch sort function.
  names(ampliconends) <- names(dnastringset)
  
  #filter minimum size
  ampliconends[width(ampliconends) >= trimf+trimr]
}
