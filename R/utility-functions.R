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

#' Take An MSA Object and trim it
#'
#' use trimaligment.pl to trim an MSA object
#'
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @export
trimMSA <- function(msa) {
  temp1  <- tempfile()
  temp2 <- tempfile()

  dnaSS  <- as(msa, "DNAStringSet")
  writeXStringSet(dnaSS, file=temp1)
  trimAlignment(temp1,temp2)

  trimmed <- readDNAStringSet(temp2)
  return(trimmed)
}
#' run the Degen script on a character vector
#'
#' @param infile
#' @param method
#' @export
#'
makedegenerate <- function(infile, method="S") {
  if (!method %in% c("S","Z","SZ")) stop("Degeneracy method is only one one of 'S','Z'. or 'SZ'")
  table        <- paste0("-t=","\'",method,"\'")
  degescript   <- system.file("Degen/Degen_v1_4.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl", degescript, infile, table)
  print(cli)
  system(cli)
}
#' make codons degenerate
#'
#' use the Degen synonymous and nonsynonymous mutations
#'
#' @param sequences
#' @param method
#' @importFrom purrr map_chr
#' @export
#'
makedegenerateseqs <- function(sequences, method="S") {
  if (!method %in% c("S","Z","SZ")) stop("Degeneracy method is only one one of 'S','Z'. or 'SZ'")
  if (!is.character(sequences)) stop("makedegenerate requires characters as input")

  data("degenS")
  data("degenZ")
  data("degenSZ")

  if (method == "S") {
    deglist <- degenS
  } else if (method == "Z") {
    deglist <- degenZ
  } else if  (method == "SZ") {
    deglist <- degenSZ
  }

  if (length(sequences) == 1) {
    # handle single string
    if (nchar(sequences) == 3) {
      return( as.character( deglist[sequences] ))
    } else if (nchar(sequences) < 3 ){
      return(sequences)
    } else if (nchar(sequences) > 3 ) {
      #get triplets, get deenerate codons, return a single string
      codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", sequences), " ")[[1]]
      degcodons <- map_chr(codons, function(x) {makedegenerateseqs(x, method = method)})
      return( paste(degcodons, collapse = ""))
    }
  } else {
    # handle a series of strings
    return(map_chr(sequences, function(x){makedegenerateseqs(x, method = method)}))
  }
}
