#' Pick Primer Pairs Automatically
#'
#' @importFrom zoo rollapply
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' 
#' @param degeprime  Required.  A \code{\link{degeprimer}} object.
#' @param keepprimers Optional. Default \code{4}. How many primerpairs to keep.
#' @param minsequences Optional. Default \code{3}. The minimum number of 
#'  sequences to be included in peak-picking. This parameter is included to 
#'  avoid the  problem where primers are chosen on regions with small numbers
#'  of sequences, typically found at the end of an alignment.
#' 
#' @return vector of ints denoting the top peaks
#' @export
autofind_primers <- function(degeprime, keepprimers=4, minsequences=3) {
  if (is.null(degeprime@primerdata)) stop("Autopicking Primers requires primer information.")
  
  # descard data where there are few sequences - i.e the beginning and ends of an alignmnet
  degedf <- data.frame(degeprime@primerdata) %>% filter(PrimerMatching >= minsequences)
  
  cutoff <- mean(degedf$coverage) + 2*sd(degedf$coverage)
  
  # calculate local maxima
  degedf$localmaxima <- rollapply(degedf$coverage, 9, function(x) which.max(x)==5, fill=NA)
  
  #return a plot highlighting the high peaks (if there are any)
  # sort to get top coverage, maximum degeneracy, and mamimum primer matching
  degedf <- degedf %>% 
    filter(localmaxima  == TRUE) %>%
    arrange(-coverage, -degeneracy, -PrimerMatching)
  
  fullcoveragecount <- table(degedf$coverage == 1)['TRUE'] 
  
  # return the positions where the best matches occur
  toppositions <- degedf[1:keepprimers,]$Pos
  
  return(toppositions)
  
  # if (fullcoveragecount > keepprimers) {
  #   return(degedf %>% filter(coverage==1))
  # } else {
  #   return(degedf[1:keepprimers,])
  # }
}