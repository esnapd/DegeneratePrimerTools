#' Add Primers to the MSA to enable visualization.
#' 
#' To add primer pairs to an MSA, see \code{add_primerpairs}.
#'
#' @param degeprime Required. A \code{degeprimer-class} object.
#' @param positions Required. A position (or vector of positions) in the 
#'    DEGEPRIMER output from which degenerate primers will be picked. 
#' @param max.mismatch Optional. Default \code{3}. Maxmimum mismatch between the 
#'    primer and a DNA target.
#' @param windowsize Optional. Default \code{30}. Windowsize of MSA to return. 
#'   A setting of '0' will return the full length alignment. Note: if multiple
#'  positions are specified, it is no longer possible to specify a window.
#' @param mode Optional. Default \code{'consensus'}. If there are multiple 
#'   matches of a degenerate primer against a target sequence, using strict
#'    ill cause the functon to fail. If 'strict' it will only allow single
#'    matches. If "first", it will allow the first match, if 'consensus', 
#'    it will choose the most abundant location.
#' @return a \code{\link[Biostrings]{DNAMultipleAlignment}}
#' @importFrom purrr map
#' @importFrom purrr discard
#' @importFrom purrr reduce
#' @importFrom purrrlyr by_row
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAMultipleAlignment
#' @importFrom Biostrings union
#' @importFrom Biostrings subseq
#' @importFrom Biostrings matchPattern
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @export
#' 
add_primers_to_MSA <- function(degeprime, positions, max.mismatch=3, windowsize=30, mode="consensus") {
  if (is.null(degeprime@primerdata)) stop("There is not primer information associated with this object")
  
  stopifnot(exists("positions"))
  print(class(positions))
  # obtain the primer sequences: find the matches between the primer and several sequences in the MSA
  primerdata <- data.frame(degeprime@primerdata) %>%
    filter(Pos %in% positions) %>%
    arrange(Pos, degeneracy)
  
  if (nrow(primerdata) == 0) stop("There are no degenerate primers at the specified positions.")
  
  msa1       <- degeprime@msa
  msawidth   <- ncol(msa1)
  
  
  # get match location between primer and MSA sequences using matchPattern
  primerdata$msaloc <- lapply(primerdata$PrimerSeq, function(prim){
    
    # get the locations of your primer against each sequence in the MSA
    primermatches <- lapply(DNAStringSet(msa1), function(x){
      p1    <- matchPattern(pattern=prim, subject=x, fixed=FALSE, max.mismatch=max.mismatch)
      p1loc <- start(p1)[1]
      if (length(p1) == 0) warning("No Matches for the Forward Primer")
      if (length(p1) > 1) warning("Multiple matches for the forward primer, using the first.")
      return(p1loc)
    })
    
    # primer location
    primerlocation <- purrr::discard(unlist(primermatches), is.na)
    if (length(unique(primerlocation))  > 1) {
      if (mode=="strict") {
        locs <- paste(unique(primerlocation), collapse=",")
        msg <- paste("There can only be a single location chosen for primers matching an MSA when using strict.",
                     "The primer sequence of interest is", prim , "and it is matching at the following locations:", locs)
        stop(msg)
      } else if (mode=="first"){
        warning("There are multiple matches of your degenerate sequence against one or more target sequences. Since strict is set to FALSE, this will return the first match.")
        primerlocation <- unique(primerlocation)[[1]]
      } else if (mode=="consensus") {
        locfreqs       <- sort( table(primerlocation), decreasing = TRUE)
        consensus      <- as.numeric(names(locfreqs[1]))
        primerlocation <- consensus
        
        locs <- paste(unique(primerlocation), collapse=",")
        msg  <- paste("You are using the consensus mode to add primer location to an MSA. One of your primers with sequence",
                     prim , "and it is matching at the following locations: ", locs, " We have chosen the consensus loction,",
                     consensus,".")
        message(msg)
        
      } else {
        stop("Invalid choice for mode. Can be 'strict', 'first', or 'consensus'")
      }
    }
    
    return(unique(primerlocation))
  })
  
  
  # calculate where on the MSA the sequences should be placed and create a sequecne with
  # dashes in the non-matching parts.
  # the new sequence to be added to the MSA will be
  # -----Primer-------
  primerdata <- mapply(
      function(pos, deg, primerdeg, primerseq, msaloc){
    
      pname     <- paste("Pos", pos[[1]],"Deg",deg[[1]], "Calcdeg", primerdeg[[1]], sep="_")
      fp        <- primerseq[[1]]
      fp_length <- nchar(fp)
      start     <- msaloc[[1]]
      end       <- start + fp_length
      
      # calculate intervals: middledash is total interval minus primers
      startlen   <- start - 1
      endlen     <- msawidth - end + 1
      startdash  <- paste(rep("-", startlen), collapse="")
      enddash    <- paste(rep("-", endlen), collapse="")
      
      dna <- DNAStringSet(paste0(startdash, fp, enddash))
      names(dna) <- pname
      
      dna
    },
    primerdata$Pos,
    primerdata$degeneracy,
    primerdata$PrimerDeg,
    primerdata$PrimerSeq,
    primerdata$msaloc
  ) 
    
  
  # combine the sequences together
  primersaligned <- purrr::reduce(primerdata, Biostrings::union)
  
  # add the new sequences ot the MSA
  dnacombined <- Biostrings::union(DNAStringSet(msa1), primersaligned)
  
  # There is not window slicing if multiple positions are bieng used
  if (length(positions > 1)) return(DNAMultipleAlignment(dnacombined))
  
  # if a single posiiton is used, allow for windowing
  if (windowsize == 0) {
    DNAMultipleAlignment(dnacombined)
  } else {
    if (windowsize < 20) warning("You are slicing the MSA using a small number of basepairs. Consider increasing your windowsize.")
    
    # center the return window on the center of the primer sequence
    pos    <- primerdata$Pos[[1]]
    center <- round(nchar(primerdata$PrimerSeq[[1]]) / 2)
    start  <- pos + center
    left   <- start - round(windowsize /2 )
    right  <- start + round(windowsize /2 )
    windowed_aln <- subseq(dnacombined, left, right)
    return( DNAMultipleAlignment(windowed_aln))
  }  
}