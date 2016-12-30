#' plot_GC
#' 
#' Plot the GC Content of the reference sequences 
#' of a degeprimer object.
#' 
#' @param degprim Required. A degeprime object or a a phyloseq object.
#' @param lowGC Optional. Default \code{0}. The low GC value for the purpose of plotting.
#' @param highGC Optional. Default \code{1}  The high GC value for the purpose of plotting.
#' @param midGC Optional. Default \code{0.5}  The midpoint GC value for the purpose of
#'  plotting. This is where the color will fade to white.
#' @importFrom ggtree gheatmap 
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tiplab
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 ggtitle
#' 
#' @export
plot_GC <- function(degprim, lowGC=0, highGC=1, midGC=0.5, ...) {
  if (!class(degprim) %in% c("degeprimer", "phyloseq")) {
    stop("The first argument must be of class 'degeprimer' or 'phyloseq'")
  }
  
  if (class(degprim)== "phyloseq") {
    if (is.null(degprim@phy_tree))
      stop("your phyloseq object must contain a phylogenetic tree")
  }
  
  gc <- get_GC(degprim@refseq)
 
  # pass the created matrix to ggtree's matrix mapping function
  p  <- ggtree(degprim@phy_tree,ladderize = T)
  p  <- p + geom_tiplab(size=3, align=FALSE)
  gheatmap(p, gc, ...) + 
    scale_fill_gradient2(limits=c(lowGC, highGC),  midpoint = midGC) +
    ggtitle("GC Content")
}



