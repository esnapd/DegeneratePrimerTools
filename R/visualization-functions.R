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

#' Plot how well a primerset covers a tree.
#'
#' note this function relies on the consistent ordering of nodes
#' in phylo objects where tips come first, then nodes in the edgelist
#'https://sites.google.com/site/phylogeneticworkshop2/labs/plotting-phylogenies
#'
#' @importFrom purrr map_chr
#' @export
plot_coverage <- function(dgprimer, primpair) {

  matches    <- find_primers(dgprimer@refseq, primpair@forwardprimer, primpair@reverseprimer)
  matchnames <- matches$sequence
  tree       <- dgprimer@phy_tree

  #names of primer pairs having a match then get the  tip indices
  matchidx <- which(tree$tip.label %in% matchnames)

  #get nodes attached to matched tips
  connectednodes    <- tree$edge[ tree$edge[,2] %in% matchidx, ]
  connectednodes    <- connectednodes[,1] #mat->vec
  connectednodelabs <- connectednodes - length(tree$tip.label)

  print(connectednodelabs)
  tree$node.label[c(2:length(tree$node.label))] <- "No_Match"
  tree$node.label[c(connectednodelabs)]         <- "Primer_Match"

  tree$tip.label <- map_chr(tree$tip.label, function(x){strsplit(x, "_")[[1]][[1]]})

  return(plot_phylotree(tree,colordomain=c("Primer_Match", "No_Match")))

}
#' draw a radial phylogenetic tree
#'
#' @import htmlwidgets
#' @export
plot_phylotree <- function(tree,
                           colordomain = c(""),
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


#' #' Update Tree With Primer Information
#' #'
#' #' Annotate a tree with Primer Information
#' #' @export
#' updateTreePrimerMatch <- function(matchdf ,tree) {
#'   # note this function relies on the consistent ordering of nodes
#'   # in phylo objects where tips come first, then nodes in the edgelist
#'
#'   #names of primer pairs having a match then get the  tip indices
#'   matches  <- matchdf[!is.na(matchdf$expected),]$sequence
#'   matchidx <- which(tree$tip.label %in% matches)
#'
#'   #get nodes attached to matched tips
#'   connectednodes    <- tree$edge[ tree$edge[,2] %in% matchidx, ]
#'   connectednodes    <- connectednodes[,1] #mat->vec
#'   connectednodelabs <- connectednodes - length(tree$tip.label)
#'
#'   #set node names and plot
#'   tree$node.label[c(2:length(tree$node.label))] <- "No_Match"
#'   tree$node.label[c(connectednodelabs)] <- "Primer_Match"
#'
#'   return(tree)
#' }
