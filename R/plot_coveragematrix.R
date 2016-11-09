#' Plot the Coverage of Degenerate primers
#'
#' @param degprim
#' @param primerpairlist
#' @param max.mismatch
#' @param tiplabelsize
#' @param align
#' @importFrom ggtree gheatmap 
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tiplab
#' @export
plot_coveragematrix <- function(degprim, primerpairlist=NULL, max.mismatch=1, tiplabelsize=3, align=TRUE, ...) {
  if (!class(degprim) == "degeprimer") {
    stop("The first argument must be of class 'degeprimer'")
  }
  if (is.null(primerpairlist)){
    primerpairlist <- degprim@primerpairs
  }
  
  if (length(primerpairlist) == 0) stop("This function requires primerpairs")
  
  refseq <- degprim@refseq
  
  # Create Matrix of Primer-Sequence Matching
  primerdata <- lapply(primerpairlist, function(ppair){
    # make primermatrix from one refseq against one primer
    pname        <- ppair@name
    hitdf        <- find_primers(refseq, fp=ppair@forwardprimer,rp=ppair@reverseprimer, max.mismatch=max.mismatch)
    hitdf[pname] <- mapply(function(start,end) {ifelse(is.na(start) || is.na(end), "No Match","Match")},
                           hitdf$start, hitdf$end)
    row.names(hitdf) <- hitdf$sequence
    hitdf[pname]
  })
  
  df <- do.call("cbind", primerdata)
  
  # pass the created matrix to ggtree's matrix mapping function
  p  <- ggtree(degprim@phy_tree, ladderize = TRUE)
  p  <- p + geom_tiplab(size = tiplabelsize, align = align)
  gheatmap(p, df, ...)
}