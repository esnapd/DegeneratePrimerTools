#' Add Primer Pairs to the MSA to enable visualization.
#'
#' @param degeprime Required. A \code{degeprimer-class} object.
#' @param max.mismatch Optional. Default \code{3}.
#' @return a \code{\link[Biostrings]{DNAMultipleAlignment}}
#' @importFrom purrr map
#' @importFrom purrr map_df
#' @importFrom purrr map_chr
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAMultipleAlignment
#' @importFrom Biostrings union
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @export
#' 
add_primerpairs_to_MSA <- function(degeprime, max.mismatch=3) {
  # find the matches between the primer and several sequences in the MSA
  primerdata <- purrr::map_df(degeprime@primerpairs, function(ppair){
    # make primermatrix from one refseq against one primer
    pname        <- ppair@name
    hitdf        <- find_primers(DNAStringSet(degeprime@msa), fp=ppair@forwardprimer,rp=ppair@reverseprimer, max.mismatch=max.mismatch)
    hitdf$Match  <- mapply(function(start,end) {ifelse(is.na(start) || is.na(end), "No Match","Match")}, hitdf$start, hitdf$end)
    hitdf$Primer <- pname
    hitdf
  })
  
  #convert primer info to tabular format
  primertable <- make_primer_table(degeprime, wide=TRUE) 
  primertable$reverse_RC <- purrr::map_chr(primertable$reverseprimer, DegeneratePrimerTools:::rc) 
  
  # merge location data and sequence data
  primerdata <- primerdata %>% 
    left_join(primertable, by=c("Primer"="primername")) %>% 
    filter(!is.na(start), !is.na(end))
  
  # add sequences to the MSA
  msa1     <- degeprime@msa
  msawidth <- ncol(msa1)
  
  # the new sequence to be added to the MSA will be
  # -----Forwardprimer-------ReversePrimer---------
  primersaligned <- purrr::map(
    split(primerdata, primerdata$Primer), #split the primerdata df by primer
    function(pdata){
      
      pname <- pdata$Primer[[1]]
      fp    <- pdata$forwardprimer[[1]]
      rp    <- pdata$reverseprimer[[1]]
      rprc  <- pdata$reverse_RC[[1]]
      
      fp_length <- nchar(fp)
      rp_length <- nchar(rp)
      
      start <- pdata$start[[1]]
      end   <- pdata$end[[1]]
      
      # calculate intervals.
      # middledash is total interval minus primers
      startlen   <- start - 1
      endlen     <- msawidth - end
      middlelen  <- (end - startlen - fp_length - rp_length)
      startdash  <- paste(rep("-", startlen), collapse="")
      enddash    <- paste(rep("-", endlen), collapse="")
      middledash <- paste(rep("-", middlelen), collapse="" )
      
      dna <- DNAStringSet(paste0(startdash, fp, middledash, rprc, enddash))
      names(dna) <- pname
      
      dna
    })
  
  primersaligned <- purrr::reduce(primersaligned, Biostrings::union)
  
  dnacombined <- Biostrings::union(DNAStringSet(msa1), primersaligned)
  
  DNAMultipleAlignment(dnacombined)
}