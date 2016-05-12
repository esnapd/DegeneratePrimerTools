################################################################################
#' Build degeprimer-class objects from their components.
#'
#' \code{degeprimer()} is a constructor method, This is the main method
#' suggested for constructing an experiment-level (\code{\link{degeprimer-class}})
#' object from its component data
#' (component data classes: \code{\link{refseq-class}}, \code{\link{msa-class}},
#'  \code{\link{phylo-class}}, \code{\link{primerdata-class}}).
#'
#' @usage degeprimer(...)
#'
#' @param ... One or more component objects including the \code{phylo}-class
#'  (defined by the \code{\link{ape-package}}), \code{MsaDNAMultipleAlignment}-class
#'  (defined by the \code{\link{msa-package}}), \code{XStringSet}-class
#'  (defined by the \code{\link{Biostrings-package}}) classes as well as the
#'  degeprimer-specific class, \code{\link{primerdata-class}}.
#'
#' @return The class of the returned object depends on the argument
#'  class(es). For an experiment-level object, two or more component data objects
#'  must be provided.
#'  Otherwise, if a single component-class
#'  is provided, it is simply returned as-is.
#'  The order of arguments does not matter.
#'
#' @export
#' @examples
#' data(woodmouse)
#' x1 = degeprimer(phy_tree(woodmouse), phy_tree(esophagus))
#' identical(x1, esophagus)
#' # # data(GlobalPatterns)
#' # # GP <- GlobalPatterns
#' # # phyloseq(sample_data(GP), otu_table(GP))
#' # # phyloseq(otu_table(GP), phy_tree(GP))
#' # # phyloseq(tax_table(GP), otu_table(GP))
#' # # phyloseq(phy_tree(GP), otu_table(GP), sample_data(GP))
#' # # phyloseq(otu_table(GP), tax_table(GP), sample_data(GP))
#' # # phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))
degeprimer <- function(...){

  arglist <- list(...)

  # Remove names from arglist. Will replace them based on their class
  names(arglist) <- NULL

  # ignore all but component data classes.
  arglist  <- arglist[sapply(arglist, is.component.class)]

  # Make the name-replaced, splatted list
  splatlist <- sapply(arglist, splat.degeprimer.objects)

  # rm any forbidden chars in index names (e.g. quotes - phylogenetic tree).
  # Right now, only extra quotes are forbidden.
  # splatlist = lapply(splatlist, function(x){
  #   taxa_names(x) <- gsub("\"", "", taxa_names(x), fixed=TRUE)
  #   taxa_names(x) <- gsub("\'", "", taxa_names(x), fixed=TRUE)
  #   return(x)
  # })

  ####################
  ## Need to determine whether to
  # (A) instantiate a new raw/uncleaned degeprimer object, or
  # (B) return a single component, or
  # (C) to stop with an error because of incorrect argument types.
  if( length(splatlist) > length(get.component.classes()) ){
    stop("Too many components provided\n")
  } else if( length(names(splatlist)) > length(unique(names(splatlist))) ){
    stop("Only one of each component type allowed.\n",
         "For merging multiple objects of the same type/class, try merge_phyloseq(...)\n")
  }
  # else if( length(splatlist) == 1){
  #   return(arglist[[1]])
  # }
  else {
    # Instantiate the degeprimer-class object, dp.
    dp <- do.call("new", c(list(Class="degeprimer"), splatlist) )
  }

  ####################
  ## Reconcile the names of the components in the
  ## refseq, msa, and tree.
  ## in the newly-minted phyloseq object
  #shared_taxa    = intersect_taxa(dp)
  #shared_samples = intersect_samples(dg)


  # if( length(shared_taxa) < 1 ){
  #   stop("Problem with OTU/taxa indices among those you provided.\n",
  #        "Check using intersect() and taxa_names()\n"
  #   )
  # }
  # if( length(shared_samples) < 1 ){
  #   stop("Problem with sample indices among those you provided.\n",
  #        "Check using intersect() and taxa_names()\n"
  #   )
  # }
  #
  # # Start with OTU indices
  # ps = prune_taxa(shared_taxa, ps)
  #
  # # Verify there is more than one component
  # # that describes samples before attempting to reconcile.
  # ps = prune_samples(shared_samples, ps)
  #
  # # Force both samples and taxa indices to be in the same order.
  # ps = index_reorder(ps, "both")
  #
  # # Replace any NA branch-length values in the tree with zero.
  # if( !is.null(phy_tree(ps, FALSE)) ){
  #   ps@phy_tree <- fix_phylo(ps@phy_tree)
  # }
  warning("No checks implemented to ensure data consistency!")


  return(dp)
}
################################################################################
#' Show the component objects classes and slot names.
#'
#' There are no arguments to this function. It returns a named character
#' when called, which can then be used for tests of component data types, etc.
#'
#' @usage get.component.classes()
#'
#' @return a character vector of the component objects classes, where each
#' element is named by the corresponding slot name in the phyloseq-class.
#'
#' @keywords internal
#'
#' @examples #
#' #get.component.classes()
get.component.classes <- function(){
  # define classes vector
  component.classes <- c("XStringSet","MsaDNAMultipleAlignment", "phylo", "primerdata")
  # the names of component.classes needs to be the slot names to match getSlots / splat
  names(component.classes) <- c("refseq", "msa", "phy_tree", "primerdata")
  return(component.classes)
}
# Returns TRUE if x is a component class, FALSE otherwise.
# This shows up over and over again in data infrastructure
#' @keywords internal
is.component.class = function(x){
  inherits(x, get.component.classes())
}
################################################################################
#' Convert \code{\link{degeprimer-class}} into a named list of its non-empty components.
#'
#' This is used in internal handling functions, and one of its key features
#' is that the names in the returned-list match the slot-names, which is useful
#' for constructing calls with language-computing functions like \code{\link{do.call}}.
#' Another useful aspect is that it only returns the contents of non-empty slots.
#' In general, this should only be used by phyloseq-package developers. Standard
#' users should not need or use this function, and should use the accessors and
#' other tools that leave the multi-component object in one piece.
#'
#' @usage splat.degeprimer.objects(x)
#'
#' @param x A \code{\link{degeprimer-class}} object. Alternatively, a component
#'  data object will work, resulting in named list of length 1.
#'
#' @return A named list, where each element is a component object that was contained
#' in the argument, \code{x}. Each element is named according to its slot-name in
#' the phyloseq-object from which it is derived.
#' If \code{x} is already a component data object,
#' then a list of length (1) is returned, also named.
#'
#' @keywords internal
#' @examples #
splat.degeprimer.objects <- function(x){
  if( is.component.class(x) ){
    # Check if class of x is among the component classes already (not phyloseq-class)
    splatx <- list(x)
    names(splatx) <- names(which(sapply(get.component.classes(), function(cclass, x) inherits(x, cclass), x)))
  } else if( inherits(x, "degeprimer") ){
    # Else, check if it inherits from degeprimer, and if-so splat
    slotnames = names(getSlots("degeprimer"))
    allslots  = sapply(slotnames, function(i, x){access(x, i, FALSE)}, x)
    splatx    = allslots[!sapply(allslots, is.null)]
  } else {
    # Otherwise, who knows what it is, silently return NULL.
    return(NULL)
  }
  return(splatx)
}
