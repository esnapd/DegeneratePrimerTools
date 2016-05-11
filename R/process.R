#' Visualize the coverage of a primer on an Interactive Pylogenetic Tree
#'
#' @importFrom Biostrings DNAString countPattern
visualize_coverage <- function(forwardseq, reverseseq, referenceseq, tree) {
  #load the data
  fp  <- DNAString(forwardseq)
  rp  <- DNAString(reverseseq)
  ref <- Map(function(x){DNAString(x)}, referenceseq)

  pattern = DNAString("TSGGTGTNT")
  matchPattern(pattern, DNAString("TGGGTGTCTTGGGTGTATTTGGTGTAT"), fixed = FALSE)
}
#'
#'
find_primers <- function(fp, rp, refseq) {
  if (!class(ks1) == "DNAString") {
    stop("Expecting a DNAString as the reference object")
  }

  #load the data
  fp  <- DNAString(forwardseq)
  rp  <- DNAString(reverseseq)


  #calculate the positions
  p1    <- matchPattern(pattern=fp, subject=ref, fixed=FALSE)
  p2    <- matchPattern(pattern=fp, subject=ref, fixed=FALSE)
  p1loc <- start(p1)[1]
  p2loc <- start(p2)[1]

  if (length(p1) == 0) stop("No Matches for the Forward Primer")
  if (length(p2) == 0) stop("No Matches for the Reverse Primer")
  if (length(p1) >  1) warning("Multiple matches for the forward primer, using the first.")
  if (length(p2) >  1) warning("Multiple matches for the forward primer, using the first.")

  return(list(start=p1loc,end=p2loc))

}
#' Split FNA by Tree Distance And ReCreate
#'
#' @import ape
#' @importFrom msa msaMuscle msaClustalW
#' @export
split_fna_tree <- function(tree,fnas, splits=2, degeneracyrange=c(1,4,10,50,100,1000), oligolength=21, ncpus=1) {

  distances   <- cophenetic.phylo(tree)
  kmeans_out  <- kmeans(distances, centers = splits)
  clusterdata <- lapply(1:splits, function(x){names(which(kmeans_out$cluster==x))})
  subfnas     <- lapply(seq(clusterdata), function(x){fna[names(fna) %in% clusterdata[[x]]]})
  #alns        <- lapply(1:splits, function(x){msaMuscle(subfnas[[x]])})
  alns        <- lapply(1:splits, function(x){msaClustalW(subfnas[[x]])})

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
