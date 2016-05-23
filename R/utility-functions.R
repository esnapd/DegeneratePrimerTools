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
#' Find Primers matching A Set of Reference Sequences
#'
#' Return the best hit from each
#' @importFrom Biostrings matchPattern reverseComplement
#' @importFrom purrr compact
#' @export
setGeneric("find_primers", function(refseq, fp, rp ) standardGeneric("find_primers"))
#'
setMethod("find_primers", c("DNAString", "DNAString", "DNAString"), function(refseq, fp, rp ) {
  # determine if the FP/RP match the sequences.
  p1    <- matchPattern(pattern=fp, subject=refseq, fixed=FALSE, max.mismatch=1)
  p2    <- matchPattern(pattern=reverseComplement(rp), subject=refseq, fixed=FALSE,max.mismatch=5)

  p1loc <- start(p1)[1]
  p2loc <- start(p2)[1]
  if (length(p1) == 0) warning("No Matches for the Forward Primer")
  if (length(p2) == 0) warning("No Matches for the Reverse Primer")
  if (length(p1) >  1) warning("Multiple matches for the forward primer, using the first.")
  if (length(p2) >  1) warning("Multiple matches for the forward primer, using the first.")

  df <- data.frame(start=p1loc,end=p2loc)
  return(df)
})
#'
#'
setMethod("find_primers", c("DNAStringSet", "DNAString", "DNAString"), function(refseq, fp, rp ) {
  # determine if the FP/RP match the sequences.
  primermatches <- lapply(refseq, function(x,fp.=fp,rp.=rp){
    return(find_primers(x, fp., rp.))
  })
  # combine and calculate stats
  df <- Reduce(rbind, primermatches)
  df$expected <- df$end - df$start
  df$sequence <- names(primermatches)
  df <- na.omit(df)
  return(df)
})


