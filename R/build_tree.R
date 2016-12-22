#' Create a Tree from An MSA
#'
#' @param dgprimer Required. A degeprimer object.
#' @param force. Optional. Default is \code{FALSE}. If an msa already exists, you must set FORCE to TRUE to overwrite it.
#'
#' @export
build_tree <- function(dgprimer,force=FALSE,...) {
  
  if (is.null(dgprimer@msa)) {
    stop("Your degeprime object does not have a MultipleSequenceAlignment. Try using run-alignment")
  }
  if (force==FALSE  & !is.null(dgprimer@phy_tree)) {
    stop("Your degeprime object already has a phylogenetic tree. To overwrite use force=TRUE")
  }
  
  tree    <-  run_fasttree(dgprimer@msa)
  dgprimer@phy_tree <- tree
  
  return(dgprimer)
}