################################################################################
#' The S4 class for storing DEGEPRIME output data.
#'
#' DEGEPRIME is used to find degenerate primers from Multiple Sequence Alignments.
#' The output dataframe is stored as a dataframe where one or more additional columns
#' can be used to store metadata.
#'
#' \describe{
#'    \item{.Data}{data-frame data, inherited from the data.frame class.}
#'  }
#' @name primerdata-class
#' @rdname primerdata-class
#' @exportClass primerdata
setClass("primerdata", contains = "data.frame")
################################################################################
#' S3 class placeholder definition (list) for phylogenetic trees.
#'
#' The ape package does not export a version of its \code{\link[ape]{phylo}}-class,
#' partly because it is not really defined formally anywhere.
#' Instead, it is an S3 class extended from the base class, \code{\link{list}} --
#' this is a very common and easy approach --
#' and proper behavior of any method taking an instance of this class
#' requires exact naming conventions for element names of the components.
#' The phyloseq package does not provide any validity checks that a given phylo
#' instance is valid (conforms to the conventions in the ape package). Yet.
#' If problems arise, this might be considered, and they could be defined
#' judiciously and within phyloseq.
#' Similarly, if a formal definition for the the phylo-class is ever exported
#' by ape, the current philosophy of phyloseq would be to remove this
#' internal definition and import the former. Note that there is still some
#' work going on for the phylobase package, which is addressing these same
#' exact issues for S4 phylogenetic tree interaction.
#' A very large number of packages (around 60 at my last count), depend on ape,
#' making it easily the de facto standard for representing phylogenetic trees in R;
#' and the phyloseq team would prefer to use any exported definitions from
#' the ape package if possible and available.
#'
#' @seealso
#' \code{\link[ape]{phylo}}
#'
#' @keywords internal
phylo <- structure(list(), class = "phylo")
################################################################################
# If this ever works
# @importClassesFrom ape phylo
################################################################################
#' An S4 placeholder of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree.
#' It is used throughout the ape package.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
#' @name phylo-class
#' @rdname phylo-class
#' @exportClass phylo
setOldClass("phylo")
################################################################################
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.
################################################################################
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
#' @importClassesFrom Biostrings BStringSet
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom Biostrings RNAStringSet
#' @importClassesFrom Biostrings AAStringSet
#' @importClassesFrom Biostrings QualityScaledXStringSet
#' @importClassesFrom Biostrings XStringQuality
#' @importClassesFrom Biostrings PhredQuality
#' @importClassesFrom Biostrings SolexaQuality
#' @importClassesFrom Biostrings IlluminaQuality
#' @importClassesFrom Biostrings QualityScaledBStringSet
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importClassesFrom Biostrings QualityScaledRNAStringSet
#' @importClassesFrom Biostrings QualityScaledAAStringSet
#' @importClassesFrom Biostrings XStringSet
#' @keywords internal
#' @keywords internal
setClassUnion("XStringSetOrNULL", c("XStringSet", "NULL"))
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom Biostrings XStringSet
#' @keywords internal
setClassUnion("DNAMultipleAlignmentOrNULL", c("DNAMultipleAlignment", "NULL"))
#' @keywords internal
setClassUnion("primerdataOrNULL", c("primerdata", "NULL"))
################################################################################
#' The main experiment-level class for degprimer data
#'
#' Contains all currently-supported component data classes:
#' \code{\link[Biostrings]{XStringSet-class}} (\code{"refseq"} slot).
#' \code{\link[Biostrings]{DNAMultipleAlignment}} (\code{"msa"} slot).
#' \code{\link[ape]{phylo}}-class (\code{"phy_tree"} slot),
#' \code{{primerdata-class}} (\code{"primerdata"} slot),
#'
#' There are several advantages to storing your primer-design data as an
#'instance of the
#' phyloseq class, not the least of which is that it is easy to return to the
#' data later and feel confident that the different data types ``belong'' to
#' one another. Furthermore, the \code{\link{phyloseq}} constructor ensures that
#' the different data components have compatible indices (e.g. OTUs and samples),
#' and performs the necessary trimming automatically when you create your
#' ``experiment-level'' object. Downstream analyses are aware of which data
#' classes they require -- and where to find them -- often making your
#' \code{phyloseq-class} object the only data argument required for analysis and plotting
#' functions (although there are many options and parameter arguments available
#' to you).
#'
#' In the case of missing component data, the slots are set to \code{NULL}. As
#' soon as a \code{phyloseq-class} object is to be updated with new component
#' data (previously missing/\code{NULL} or not), the indices of all components
#' are re-checked for compatibility and trimmed if necessary. This is to ensure
#' by design that components describe the same taxa/samples, and also that these
#' trimming/validity checks do not need to be repeated in downstream analyses.
#'
#' slots:
#' \describe{
#'    \item{refseq}{ a biological sequence set object of a class that
#'         inherits from the \code{\link[Biostrings]{XStringSet-class}}, from the Biostrings package.}
#'    \item{msa}{ a single object of the \code{\link[Biostrings]{DNAMultipleAlignment-class}}, from the Biostrings package.}
#'    \item{phy_tree}{ a single object of the \code{\link[ape]{phylo}}-class, from the ape package.}
#'    \item{phy_tree}{ a single object of the \code{primerdata-class}}.}
#' }
#' @seealso
#'  The constructor, \code{\link{degeprimer}} and  the component
#'  constructor/accessors \code{\link{otu_table}}, \code{\link{sample_data}},
#'  \code{\link{tax_table}}, \code{\link{phy_tree}}, and \code{\link{refseq}}.
#'
#' @import BiocGenerics
#' @importClassesFrom Biostrings XStringSet
#' @name degeprimer-class
#' @rdname degeprimer-class
#' @exportClass degeprimer
setClass(Class="degeprimer",
         representation=representation(
           refseq     = "XStringSetOrNULL",
           msa        = "DNAMultipleAlignmentOrNULL",
           phy_tree   = "phyloOrNULL",
           primerdata = "primerdataOrNULL"),
         prototype=prototype(refseq=NULL, msa=NULL,phy_tree=NULL, primerdata=NULL)
)
################################################################################
