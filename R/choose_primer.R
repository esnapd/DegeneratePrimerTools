#' choose a primer from the DEGEPRIMER output
#'
#' a shiny widget to visualize the output form DEGEPRIMER. This interactive plot puts
#' the position of the MSA on the x-axis and coverage on the y-axis. Points can be selected
#' by dragging the mouse. The point/primer information is returned to the R console.
#'
#' @import miniUI
#' @import ggplot2
#' @importFrom  shiny plotOutput runGadget dialogViewer fillRow checkboxGroupInput
#' @export
setGeneric("choose_primer", function(object) standardGeneric("choose_primer"))
#'
#' @export
setMethod("choose_primer", "degeprimer", function(object){choose_primer(object@primerdata)})
#'
#' @export
setMethod("choose_primer", "primerdata", function(object){
  df <- as.data.frame(object)
  
  ui <- miniPage(
    gadgetTitleBar("Drag to Select Primer", right = miniTitleBarButton("done", "Done", primary = TRUE)),
    
    miniContentPanel(
      fillRow(flex = c(1, 3),
              
              #degeneracy values
              checkboxGroupInput('show_degen', 'Primer degeneracy values to show:',
                                 unique(df$degeneracy), selected = max(unique(df$degeneracy))),
              
              # The brush="brush" argument means we can listen for
              # brush events on the plot using input$brush.
              plotOutput("plot", height = "100%", brush = "brush")
      )))
  
  
  server <- function(input, output, session) {
    # Render the plot
    output$plot <- renderPlot({
      # Plot the data with x/y vars indicated by the caller.
      ggplot(df[df$degeneracy %in% input$show_degen,] , aes(Pos, coverage, color=PrimerDeg)) + geom_point()
    })
    
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(data, input$brush))
    })
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(df, input$brush))
    })
  }
  runGadget(ui,
            server,
            viewer =  dialogViewer("Pick Points",
                                   width = 1000,
                                   height = 700))
})