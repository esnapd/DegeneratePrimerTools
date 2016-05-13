#' peakfinder
#'
#' obtain the peaks from a DEGEPRIME dataframe
#'
#' @import miniUI
#' @import ggplot2
#' @importFrom  shiny plotOutput runGadget dialogViewer
#' @export
peakfinder <- function(data, xvar, yvar, width=700, height=400){
  ui <- miniPage(
    gadgetTitleBar("Drag to Select", right = miniTitleBarButton("done", "Done", primary = TRUE)),

    miniContentPanel(
      # The brush="brush" argument means we can listen for
      # brush events on the plot using input$brush.
      plotOutput("plot", height = "100%", brush = "brush")
    ))

  server <- function(input, output, session) {
    #observeEvent(input$done, { stopApp(input$moleculedata) })

    # Render the plot
    output$plot <- renderPlot({
      # Plot the data with x/y vars indicated by the caller.
      ggplot(data, aes_string(xvar, yvar)) + geom_point()
    })

    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(data, input$brush))
    })

  }

  runGadget(ui,
            server,
            viewer =  dialogViewer("Pick Points",
                                   width = width,
                                   height = height))

}
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
#' draw a radial phylogenetic tree
#'
#' @import htmlwidgets
#' @export
plot_phylotree <- function(tree,
                           colordomain = c("Bacteria", "Eukaryota", "Archaea"),
                           outerradius = 480,
                           innerradius = 310, width = NULL, height = NULL) {

  if (class(tree)=="phylo"){
    tree <- write.tree(tree)
  }

  if (!is.character(colordomain)) {
    colordomain <- c("")
  }

  # forward options using x
  params = list(
    tree = tree,
    outerradius = outerradius,
    innerradius = innerradius,
    colordomain = colordomain
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'phylogenetictree',
    x = params,
    width = width,
    height = height,
    package = 'DegeneratePrimerTools'
  )
}
#' Widget output function for use in Shiny
#'
#' @export
phylogenetictreeOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'phylogenetictree', width, height, package = 'DegeneratePrimerTools')
}
#' Widget render function for use in Shiny
#'
#' @export
renderPhylogenetictree <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, phylogenetictreeOutput, env, quoted = TRUE)
}
