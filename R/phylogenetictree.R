#' draw a radial phylogenetic tree
#'
#' @import htmlwidgets
#'
#' @export
phylogenetictree <- function(tree, colordomain, outerradius = 480,
                             innerradius = 310, width = NULL, height = NULL) {

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
    package = 'phylogenetictree'
  )
}

#' Widget output function for use in Shiny
#'
#' @export
phylogenetictreeOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'phylogenetictree', width, height, package = 'phylogenetictree')
}

#' Widget render function for use in Shiny
#'
#' @export
renderPhylogenetictree <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, phylogenetictreeOutput, env, quoted = TRUE)
}
