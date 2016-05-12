#' Visualize the coverage of a primer on an Interactive Pylogenetic Tree
#' @export
visualize_coverage <- function(matchdf, tree) {
  # note this funciton relies on the consistent ordering of nodes
  # in phylo objects where tips come first, then nodes in the edgelist

  #names of primer pairs having a match then get the  tip indices
  matches  <- matchdf[!is.na(matchdf$expected),]$sequence
  matchidx <- which(tree$tip.label %in% matches)

  #get nodes attached to matched tips
  connectednodes <- tree$edge[ tree$edge[,2] %in% matchidx, ]
  connectednodes <- connectednodes[,1] #mat->vec
  connectednodelabs <- connectednodes - length(tree$tip.label)

  print(connectednodelabs)
  tree$node.label[c(2:length(tree$node.label))] <- "No_Match"
  tree$node.label[c(connectednodelabs)] <- "Primer_Match"

  tree$tip.label <- purrr::map_chr(tree$tip.label, function(x){strsplit(x, "_")[[1]][[1]]})

  phylogenetictree(tree,colordomain=c("Primer_Match", "No_Match"))
  return(tree)
}
#' Update Tree With Primer Information
#'
#' Annotate a tree with Primer Information
#' @export
updateTreePrimerMatch <- function(matchdf ,tree) {
  # note this function relies on the consistent ordering of nodes
  # in phylo objects where tips come first, then nodes in the edgelist

  #names of primer pairs having a match then get the  tip indices
  matches  <- matchdf[!is.na(matchdf$expected),]$sequence
  matchidx <- which(tree$tip.label %in% matches)

  #get nodes attached to matched tips
  connectednodes    <- tree$edge[ tree$edge[,2] %in% matchidx, ]
  connectednodes    <- connectednodes[,1] #mat->vec
  connectednodelabs <- connectednodes - length(tree$tip.label)

  #set node names and plot
  tree$node.label[c(2:length(tree$node.label))] <- "No_Match"
  tree$node.label[c(connectednodelabs)] <- "Primer_Match"

  return(tree)
}
#' Find Primers matching A Set of Reference Sequences
#'
#' @importFrom Biostrings matchPattern
#' @export
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
