#' Visualize the coverage of a primer on an Interactive Pylogenetic Tree
#'
visualize_coverage <- function(matchdf, tree) {


  pattern = DNAString("TSGGTGTNT")
  matchPattern(pattern, DNAString("TGGGTGTCTTGGGTGTATTTGGTGTAT"), fixed = FALSE)
}
#' Find Primers matching A Set of Reference Sequences
#'
#' @importFrom Biostrings matchPattern
find_primers <- function(fp, rp, refseq) {
  if (!class(fp) == "DNAString") {
    stop("Expecting a DNAString as the reference object")
  }
  if (!class(rp) == "DNAString") {
    stop("Expecting a DNAString as the reference object")
  }
  if (!class(refseq) == "DNAStringSet") {
    stop("Expecting a DNAString as the reference object")
  }

  # determine if the FP/RP match the sequences.
  primermatches <- lapply(refseq, function(x){
    p1    <- matchPattern(pattern=fp, subject=x, fixed=FALSE)
    p2    <- matchPattern(pattern=rc(rp), subject=x, fixed=FALSE)
    p1loc <- start(p1)[1]
    p2loc <- start(p2)[1]
    if (length(p1) == 0) warning("No Matches for the Forward Primer")
    if (length(p2) == 0) warning("No Matches for the Reverse Primer")
    if (length(p1) >  1) warning("Multiple matches for the forward primer, using the first.")
    if (length(p2) >  1) warning("Multiple matches for the forward primer, using the first.")
    return(data.frame(start=p1loc,end=p2loc))
  })

  # combine and calculate stats
  df          <- Reduce(rbind, primermatches)
  df$expected <- df$end - df$start
  df$sequence <- names(primermatches)
  return(df)
}
#' Split FNA by Tree Distance And run Degenera On Each Subset
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
