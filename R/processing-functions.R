#' Run a Multiple Sequence Alignment
#'
#' @importFrom msa msa
#' @export
run_alignment <- function(dgprimer, method="ClustalW", maxiters="default", force=FALSE,...) {

  if (force==FALSE  & !is.null(dgprimer@msa)) {
    stop("Your degeprime object already has an MSA. To overwrite use force=TRUE")
  }

  seqs <- dgprimer@refseq
  aln  <- msa::msa(seqs, method=method, maxiters=maxiters, type="dna", ...)

  dgprimer@msa <- as(aln, "DNAMultipleAlignment")
  return(dgprimer)
}
#' Create a Tree from An MSA
#'
#' @importFrom ape bionjs as.DNAbin dist.dna
#' @export
build_tree <- function(dgprimer, method="ClustalW", maxiters="default", force=FALSE,...) {

  if (is.null(dgprimer@msa)) {
    stop("Your degeprime object does not have a MultipleSequenceAlignment. Try using run-alignment")
  }
  if (force==FALSE  & !is.null(dgprimer@phy_tree)) {
    stop("Your degeprime object already has an MSA. To overwrite use force=TRUE")
  }

  aln     <- as.DNAbin(dgprimer@msa)
  nucdist <- dist.dna(aln)
  tree    <- bionjs(nucdist)

  dgprimer@phy_tree <- tree

  return(dgprimer)
}
#' Split FNA by Tree Distance And run Degenera On Each Subset
#'
#' @import ape
#' @importFrom msa msaMuscle msaClustalW
#' @importFrom seqinr dist.alignment
#' @export
split_fna_tree <- function(tree,fnas, splits=2, degeneracyrange=c(1,4,10,50,100,1000), oligolength=21, ncpus=1) {

  distances   <- cophenetic.phylo(tree)
  kmeans_out  <- kmeans(distances, centers = splits)
  clusterdata <- lapply(1:splits, function(x){names(which(kmeans_out$cluster==x))})
  subfnas     <- lapply(seq(clusterdata), function(x){fnas[names(fnas) %in% clusterdata[[x]]]})
  #alns        <- lapply(1:splits, function(x){msaMuscle(subfnas[[x]])})
  alns        <- lapply(1:splits, function(x){msaClustalW(subfnas[[x]])})
  #alndists    <- lapply(alns, function(x){dist.alignment(x)})

  trees       <- lapply(1:splits, function(x){msaClustalW(subfnas[[x]])})

  #run degeprime on each subset
  alnfiles       <- lapply(alns, function(x) tempfile())
  trimfiles      <- lapply(alns, function(x) tempfile())
  degeprimefiles <- lapply(alns, function(x) tempfile())

  primerdata <- lapply(1:splits, function(i) {
    dnaSS  <- as(alns[[i]], "DNAStringSet")
    writeXStringSet(dnaSS, file=alnfiles[[i]])
    trimAlignment(alnfiles[[i]],trimfiles[[i]])
    primdata <- rundegeprime(alignmentfile = trimfiles[[i]],
                 ncpus = ncpus,
                 degeneracyrange = degeneracyrange,
                 oligolength = oligolength)
    primdata$subtree <- paste0("Subtree_", i)
    return(primdata)})

   return(list(fnas=subfnas, primerdata=primerdata))

}
#' Run and load degeprime and return the dataframe of results
#'
#' @importFrom parallel mclapply
#' @export
run_degeprime <- function(alignmentfile, oligolength, degeneracyrange=c(1,4,100,400,1000),
                          minimumdepth=1, skiplength=20, number_iterations=100, ncpus=1) {

  # use degeneracy range to determine the nubmer of jobs
  drange    <- seq(degeneracyrange)
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
#' choose a primer from the DEGEPRIMER output
#'
#'
#'
#' @import miniUI
#' @import ggplot2
#' @importFrom  shiny plotOutput runGadget dialogViewer fillRow checkboxGroupInput
#' @export
setGeneric("choose_primer", function(obj) standardGeneric("choose_primer"))
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
