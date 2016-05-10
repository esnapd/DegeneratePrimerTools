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
  #names = c("Pos","TotalSeq","UniqueMers","Entropy", "PrimerDeg", "PrimerMatching", "PrimerSeq")
  df = read.table(infile, header=TRUE)
  return(df)
}
#' Run and load degeprime and return the dataframe of results
#'
#' @importFrom parallel mclapply
#' @export
rundegeprime <- function(alignmentfile, oligolength, degeneracyrange=c(1,4,100,400,1000),
                         minimumdepth=1, skiplength=20, number_iterations=100, ncpus=1) {

  # use degeneracy range to determine the nubmer of jobs
  drange <- seq(degeneracyrange)
  tempfiles <- lapply(drange, function(x) {tempfile()})

  degendata <- mclapply(drange, function(x) {
    #get per-run data
    outputfile    <- tempfiles[[x]]
    maxdegeneracy <-  degeneracyrange[[x]]
    #calculate degeneracies
    degePrime(alignmentfile=alignmentfile, outfile=outputfile,
              oligolength=oligolength, maxdegeneracy = maxdegeneracy,
              minimumdepth=minimumdepth, skiplength=skiplength, number_iterations=number_iterations)
    #load the file and return a dataframe
    df <- loadDegePrimeDF(outputfile)
    df$degeneracy <- maxdegeneracy
    return(df)
  }, mc.cores=ncpus)

  #aggregate the data
  aggdata <- Reduce(rbind, degendata)

  #add coverage information
  aggdata$coverage <- aggdata$PrimerMatching/aggdata$TotalSeq

  return(aggdata)
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
#' Plot the output from rundegeprime
#'
#' @import ggplot2
#' @importFrom zoo rollapply
#' @importFrom ggrepel geom_text_repel
#' @export
plot_deg <- function(degedf) {
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
