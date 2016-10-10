#' Plot the output from rundegeprime
#'
#' @import ggplot2
#' @importFrom zoo rollapply
#' @importFrom ggrepel geom_text_repel
#' @export
plot_degeprimer <- function(degedf) {
  cutoff <- mean(degedf$coverage) + 2*sd(degedf$coverage)
  degedf$localmaxima <- rollapply(degedf$coverage, 9, function(x) which.max(x)==5, fill=NA)

  #return a plot highlighting the high peaks (if there are any)
  gg <- ggplot(degedf, aes(x=Pos,y=coverage, color=PrimerDeg)) + geom_point()
  gg <- gg + facet_grid(degeneracy~.) + theme_bw()
  gg + geom_text_repel(
    data= subset(degedf, coverage > cutoff & localmaxima == TRUE),
    aes(label = Pos)
  )
}
#' Plot the Coverage of Degenerate primers
#'
#' @param degprim
#' @param primerpairlist
#' @importFrom ggtree gheatmap ggtree geom_tiplab
#' @export
plot_coveragematrix <- function(degprim, primerpairlist=NULL, max.mismatch=1, ...) {
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

  if (ncol(df) == 1) {
    df <- cbind(df,df)
    names(df)[[2]] <- "dummy"
  }

  # pass the created matrix to ggtree's matrix mapping function
  p  <- ggtree(degprim@phy_tree,ladderize = T)
  p  <- p + geom_tiplab(size=3, align=FALSE)
  gheatmap(p, df, ...)
}
#' Plot the GC Content
#'
#' @param degprim
#' @importFrom ggtree gheatmap ggtree geom_tiplab
#' @importFrom Biostrings alphabetFrequency
#' @importFrom ggplot2 scale_fill_gradient2
#' @export
plot_GC <- function(degprim,  ...) {
  if (!class(degprim) == "degeprimer") {
    stop("The first argument must be of class 'degeprimer'")
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
  gheatmap(p, gc, ...) + scale_fill_gradient2(limits=c(0, 1), midpoint = 0.5)
}
