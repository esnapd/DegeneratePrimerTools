#' Split FNA by Tree Distance And run DEGEPRIME On Each Subset
#'
#' @import ape
#' @importFrom msa msaMuscle msaClustalW
#' @importFrom seqinr dist.alignment
#' @export
split_fna_tree <- function(tree, fnas, splits=2, degeneracyrange=c(1,4,10,50,100,1000), oligolength=21, ncpus=1) {
  
  distances   <- cophenetic.phylo(tree)
  kmeans_out  <- kmeans(distances, centers = splits)
  clusterdata <- lapply(1:splits, function(x){names(which(kmeans_out$cluster==x))})
  subfnas     <- lapply(seq(clusterdata), function(x){fnas[names(fnas) %in% clusterdata[[x]]]})
  alns        <- lapply(1:splits, function(x){msaClustalW(subfnas[[x]])})
  
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