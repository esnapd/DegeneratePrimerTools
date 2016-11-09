#' Add Primer
#'
#' A convenience function to add a primer pair to the primerl list based on the
#' degenracy and location of forward and reverse primers.
#'
#' @importFrom purrr map_chr
#' @export
add_primerpair <-function(dgprimer, name, fpos, fdeg, rpos, rdeg) {
  
  pdf           <- dgprimer@primerdata
  existingpairs <- dgprimer@primerpairs
  
  #check degenracy and position values
  if (!fpos %in% unique(pdf$Pos))        stop("fpos value invalid. must be a Pos value present in the degeprimer table.")
  if (!rpos %in% unique(pdf$Pos))        stop("rpos value invalid. must be a Pos value present in the degeprimer table.")
  if (!fdeg %in% unique(pdf$degeneracy)) stop("fdeg value invalid. must be a degeneracy value present in the degeprimer table.")
  if (!rdeg %in% unique(pdf$degeneracy)) stop("rdeg value invalid. must be a degeneracy value present in the degeprimer table.")
  
  # check name uniqueness
  if(length(existingpairs) == 0) {
    existingnames <- c()
  } else {
    existingnames <- map_chr(existingpairs, function(x) {x@name})
  }
  
  if (name %in% existingnames) {
    stop("Primer-pair names must be unique.")
  }
  
  #get sequences
  fseq <-     pdf[pdf$Pos==fpos & pdf$degeneracy==fdeg,]$PrimerSeq
  rseq <-  rc(pdf[pdf$Pos==rpos & pdf$degeneracy==rdeg,]$PrimerSeq)
  
  newprimer <- new("primerpair",
                   name=name,
                   forwardprimer  = DNAString(fseq),
                   reverseprimer  = DNAString(rseq))
  
  # Add the primer pair to the primerlist
  num_existing_primers <- length(existingpairs)
  
  if (num_existing_primers == 0) {
    plist <- list(newprimer)
  } else {
    plist <- as.list( c(existingpairs, newprimer))
  }
  
  # add the new list to the primer list slot
  dgprimer@primerpairs <- plist
  
  return(dgprimer)
}