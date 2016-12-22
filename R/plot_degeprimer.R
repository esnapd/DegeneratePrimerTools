#' plot degeprimer
#' 
#' Plot the output from rundegeprime. THe x-axis is the location along the multiple
#' sequence alignment and the y-axis is the coverage of the primer at that location
#' for a given maximum degeneracy value.
#'
#' @param deg Required. 
#' 
#' @import ggplot2
#' @importFrom zoo rollapply
#' @importFrom ggrepel geom_text_repel
#' @export
plot_degeprimer <- function(deg) {
  
  if (class(deg) ==  "degeprimer") deg <- deg@primerdata
  
  cutoff <- mean(deg$coverage) + 2*sd(deg$coverage)
  deg$localmaxima <- rollapply(deg$coverage, 9, function(x) which.max(x)==5, fill=NA)
  
  #return a plot highlighting the high peaks (if there are any)
  gg <- ggplot(deg, aes(x=Pos,y=coverage, color=PrimerDeg)) + geom_point()
  gg <- gg + facet_grid(degeneracy~.) + theme_bw()
  gg + geom_text_repel(
    data = subset(deg, coverage > cutoff & localmaxima == TRUE),
    aes(label = Pos))
}
