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
#' Find Primers matching A Set of Reference Sequences
#'
#' Return the best hit from each
#' @importFrom Biostrings matchPattern reverseComplement
#' @importFrom purrr compact
#' @export
setGeneric("find_primers", function(refseq, fp, rp ) standardGeneric("find_primers"))
#'
setMethod("find_primers", c("DNAString", "DNAString", "DNAString"), function(refseq, fp, rp ) {
  # determine if the FP/RP match the sequences.
  p1    <- matchPattern(pattern=fp, subject=refseq, fixed=FALSE)
  p2    <- matchPattern(pattern=reverseComplement(rp), subject=refseq, fixed=FALSE)

  p1loc <- start(p1)[1]
  p2loc <- start(p2)[1]
  if (length(p1) == 0) warning("No Matches for the Forward Primer")
  if (length(p2) == 0) warning("No Matches for the Reverse Primer")
  if (length(p1) >  1) warning("Multiple matches for the forward primer, using the first.")
  if (length(p2) >  1) warning("Multiple matches for the forward primer, using the first.")

  df <- data.frame(start=p1loc,end=p2loc)
  return(df)
})
#'
#'
setMethod("find_primers", c("DNAStringSet", "DNAString", "DNAString"), function(refseq, fp, rp ) {
  # determine if the FP/RP match the sequences.
  primermatches <- lapply(refseq, function(x,fp.=fp,rp.=rp){
    return(find_primers(x, fp., rp.))
  })
  # combine and calculate stats
  df <- Reduce(rbind, primermatches)
  df$expected <- df$end - df$start
  df$sequence <- names(primermatches)
  df <- na.omit(df)
  return(df)
})
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
