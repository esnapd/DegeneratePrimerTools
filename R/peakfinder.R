#' peakfinder
#'
#' obtain the peask from a DEGEPRIME dataframe
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

