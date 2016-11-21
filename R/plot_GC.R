#' Plot the GC Content 
#' 
#'
#' @param degprim. Required. A degeprime object or a a phyloseq object.
#' @importFrom ggtree gheatmap 
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tiplab
#' @importFrom Biostrings alphabetFrequency
#' @importFrom ggplot2 scale_fill_gradient2
#' @export
plot_GC <- function(degprim,  ...) {
  if (!class(degprim) %in% c("degeprimer", "phyloseq")) {
    stop("The first argument must be of class 'degeprimer' or 'phyloseq'")
  }
  
  if (class(degeprim)== "phyloseq") {
    if (is.null(degeprim@phy_tree))
      stop("your phylsoeq object must contain a phylogenetic tree")
  }
  
  refseq <- degprim@refseq
  
  alph <- Biostrings::alphabetFrequency(refseq)
  gc   <-   mapply(FUN = function(a,c,g,t) { 
    (g + c)/ (a + c + t + g)
  },
  alph[,1], alph[,2], alph[,3], alph[,4] )
  
  gc <- data.frame(gc)
  row.names(gc) <- names(refseq)
  
  # pass the created matrix to ggtree's matrix mapping function
  p  <- ggtree(degprim@phy_tree,ladderize = T)
  p  <- p + geom_tiplab(size=3, align=FALSE)
  gheatmap(p, gc, ...) + 
    scale_fill_gradient2(limits=c(0, 1), 
                         midpoint = 0.5) +
    ggtitle("GC Content")
}

