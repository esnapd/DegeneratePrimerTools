#' plot_coveragematrix
#'
#' Plot the Coverage of Degenerate primers
#'
#' @param degprim Required.
#' @param primerpairlist Optional. Default \code{NULL}. Otherwise we expect a list of \code{primerpair}s.
#' @param max.mismatch Optional. Default \code{3}. The maximum mismatch allowed between a primer and target sequence to
#' be considered as a match.
#' @param tiplabelsize Optional. Default \code{3}
#' @param align Optional. Default \code{TRUE}
#' @importFrom ggtree gheatmap
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tiplab
#' @importFrom Biostrings DNAStringSet
#' @export
plot_coveragematrix <- function(degprim, primerpairlist=NULL, max.mismatch=3, tiplabelsize=3, align=TRUE, ...) {
  if (!class(degprim) %in% c("degeprimer", "phyloseq", "DNAStringSet")) {
    stop("The first argument must be of class 'degeprimer', 'phyloseq', or 'DNAStringSet'")
  }
  if (is.null(primerpairlist)){
    message("no primerpairs specified. attmepting to find primer pairs.")
    try(primerpairlist <- degprim@primerpairs)
  }

  if (length(primerpairlist) == 0) stop("Primer pairs are not detected. This function requires primerpairs")

  try(refseq <- degprim@refseq)
  if (class(degprim)== "DNAStringSet") refseq <- degprim

  if (!exists("refseq")) {stop("There is no associated reference sequence with your input object.")}

  # Create Matrix of Primer-Sequence Matching
  primerdata <- lapply(primerpairlist, function(ppair){does_primer_match(dna=refseq, primerpair=ppair,  max.mismatch=max.mismatch)})
  df <- do.call("cbind", primerdata)

  # get the phylogenetic tree
  if (class(degprim) == "DNAStringSet") {
    aln  <- run_muscle(degprim)
    tree <- run_fasttree(aln)

  } else {
    tree <- degprim@phy_tree
  }

  # pass the created matrix to ggtree's matrix mapping function
  p  <- ggtree(tree, ladderize = TRUE)
  p  <- p + geom_tiplab(size = tiplabelsize, align = align)
  gheatmap(p, df, ...)
}
#' does_primer_match
#'
#' match a primerpair against a \code{\link[Biostrings]{DNAStringSet}} and
#' return a single column matrix of values that determine whether a primer
#' pair matches that sequence or not.
#'
#' @param dna Required. A \code{\link[Biostrings]{DNAStringSet}}
#' @param primerpair Required. A \code{primerpair}.
#' @param max.mismatch Optional. Default \code{3}.
#'
#' @return a single column matrix of whether
#'
#' @export
does_primer_match <- function(dna, primerpair, max.mismatch=3) {
  if (!class(dna) == "DNAStringSet") stop("dna must be a DNAStringSet")
  if (!class(primerpair) == "primerpair") stop("primerpair must be a primerpair")

  pname        <- primerpair@name
  hitdf        <- find_primers(dna, fp=primerpair@forwardprimer,rp=primerpair@reverseprimer, max.mismatch=max.mismatch)
  hitdf[pname] <- mapply(function(start,end) {ifelse(is.na(start) || is.na(end), "No Match","Match")},
                         hitdf$start, hitdf$end)
  row.names(hitdf) <- hitdf$sequence
  hitdf[pname]
}
