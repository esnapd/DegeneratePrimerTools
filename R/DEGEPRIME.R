#' trim a multiplesequence alignmentfile of gaps
#'
#' @param infile
#' @param outfile
#' @param minoccupancy
#' @pram refsequence
#' @param trailgap
#' @export
#'
trimAlignment <- function(infile, outfile, minoccupancy=0, refsequence = NULL, trailgap = FALSE) {
  minoccupancy <- ifelse(minoccupancy > 0, paste0("-min ", minoccupancy), " ")
  refsequence  <- ifelse(is.null(refsequence), "", paste0("-ref ", refsequence))
  trailgap     <- ifelse(trailgap==FALSE, "","-trailgap")
  degescript   <- system.file("DEGEPRIME/TrimAlignment.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl", degescript, "-i", infile, "-o", outfile, minoccupancy, refsequence, trailgap)
  print(cli)
  system(cli)
}

#' run the DegePrime.pl script
#'
#' @param alignmentfile
#' @param outfile
#' @param oligolength
#' @pram maxdegeneracy
#' @param minimumdepth
#' @param skiplength
#' @param number_iterations
#' @export
degePrime <- function(alignmentfile, outfile,  oligolength, maxdegeneracy, minimumdepth=1, skiplength=20, number_iterations=100) {
  if (!is.numeric(oligolength)) stop("Oligolength must be a number")
  if (!is.numeric(maxdegeneracy)) stop("maxdegenderacy must be a number")

  degescript   <- system.file("DEGEPRIME/DegePrime.pl", package="DegeneratePrimerTools")
  cli          <- paste("perl",   degescript,
                        "-i",     alignmentfile,
                        "-o",     outfile,
                        "-l",     oligolength,
                        "-d",     maxdegeneracy,
                        "-depth", minimumdepth,
                        "-skip",  skiplength,
                        "-iter",  number_iterations)
  print(cli)
  system(cli)
}

#' loadDegePrimerDF
#'
#' load and verify the output of a DEGEPRIME table
#'
#' @param infile
#' @export
loadDegePrimeDF <- function(infile){
  names = c("")
  df = read.table(infile)
  return(df)
}

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

