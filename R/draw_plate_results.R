#' Function to draw results on 96 well plate layout.
#'
#' @param nbWells Required. Number of wells for the plate, default=96 (integer)
#' @param plateData Required. Data to be displayed in the plate (dataframe)
#' 
#' @import ggplot2
#' @export
draw_plate_results <- function(nbWells = 96, plateData){
  
  ggplot(plateData, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
    geom_point(aes(colour = Molecule), size =10)  + theme_bw() +
    labs(x=NULL, y = NULL) + ggtitle("96 well plate layout") + 
    theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", size=22, hjust=.5))
  
}
