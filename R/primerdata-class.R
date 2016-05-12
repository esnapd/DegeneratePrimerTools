################################################################################
#' Build or access the primerdata table.
#'
#' This is the suggested method for both constructing and accessing
#' \code{\link{primerdata-class}}) objects.
#'
#' @usage primerdata(object,errorIfNULL=TRUE)
#'
#' @param object (Required). An data.frame, \code{\link{primerdata-class}},
#'  or \code{\link{degeprimer-class}}.
#'
#' @return An \code{\link{primerdata-class}} object.
#'
#' @docType methods
#' @rdname primerdata-methods
#' @export
setGeneric("primerdata", function(object){standardGeneric("primerdata")})
# Access the primerdata slot.
#' @aliases primerdata,degeprimer-method
#' @rdname primerdata-methods
setMethod("primerdata", "degeprimer", function(object){object@primerdata})
# Return the primer data as is.
#' @aliases primerdata,degeprimer-method
#' @rdname primerdata-methods
setMethod("primerdata", "primerdata", function(object){object})
# Instantiate a primerdata table from a dataframe.
#' @aliases primerdata,data.frame-method
#' @rdname primerdata-methods
setMethod("primerdata", "data.frame", function(object){
  # check names to match degeprimer headings
  degnames <- c("Pos", "TotalSeq","UniqueMers","Entropy", "PrimerDeg",
                "PrimerMatching", "PrimerSeq", "degeneracy", "coverage")
  objnames <- names(object)
  for (dn in degnames) {
    if (!dn %in% objnames) {
      stop("Dataframe names do not match DEGEPRIME output. Is your primer dataframe correct?")
    }
  }

  primdata <- new("primerdata", object)
  return(primdata)
})
