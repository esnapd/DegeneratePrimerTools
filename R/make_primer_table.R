#' make primer table
#'
#' @importFrom purrr map_df
#' @export
#' 
setGeneric("make_primer_table", function(object, wide) standardGeneric("make_primer_table"))
#'
setMethod("make_primer_table", "degeprimer", function(object, wide=TRUE) {
  primerpairs <- object@primerpairs
  
  if (length(primerpairs) == 0) stop("There are no primers defined. Add primers first.")
  
  primerdata <- map_df(primerpairs, function(p) {
    df = data.frame(
      primername = p@name,
      forwardprimer = as.character(p@forwardprimer),
      reverseprimer = as.character(p@reverseprimer),
      expectedlength = p@expectedsize,
      stringsAsFactors = FALSE)
  })
  
  if (wide) return(primerdata)
  
  # if tall then stack forward and reverse primers on top of one another
  primerdataF <- primerdata
  primerdataR <- primerdata
  
  primerdataF$primername <- paste0(primerdataF$primername, "_F")
  primerdataF <- primerdataF[, c("primername", "forwardprimer")]
  names(primerdataF) <- c("primername", "primer")
  
  primerdataR$primername <- paste0(primerdataR$primername, "_R")
  primerdataR <- primerdataR[, c("primername", "reverseprimer")]
  names(primerdataR) <- c("primername", "primer")
  
  return(rbind(primerdataF, primerdataR))
})