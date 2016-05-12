#' Reverse Complement
#' lifted from https://github.com/benjjneb/dada2/blob/master/R/misc.R
#'
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @export
rc <- function(sqs) {
  if(length(sqs) < 1) {
    return(character(0))
  } else if(length(sqs) == 1) {
    as(reverseComplement(DNAString(sqs)), "character")
  } else {
    as(reverseComplement(DNAStringSet(sqs)), "character")
  }
}
