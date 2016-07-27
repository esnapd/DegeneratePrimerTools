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
#' Split FNA by Tree Distance And run Degenera On Each Subset
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
#' Design Primers
#'
#' A convenience function to find degenerate primers usign the
#' seqeuences of the degenerateprimers object.  Return teh resutls to
#' the primerdata slot.
#' @importFrom parallel mclapply
#' @export
design_primers <- function(dgprimer, oligolength=21, degeneracyrange=c(1,4,100,400,1000),
                          minimumdepth=1, skiplength=20, number_iterations=100, ncpus=1,
                          force=FALSE) {

  # check if degepreime data ia associated with the class. if so, require force=TRUE
  if (!is.null(dgprimer@primerdata) & force==FALSE) {
    stop("Degenerate primer data alread exists. To overwrite please use force=TRUE")
  }

  # check if degepreime data ia associated with the class. if so, require force=TRUE
  if (is.null(dgprimer@msa)) {
    stop("Primer Calculation requires a multiple sequence alignment")
  }

  # use degeneracy range to determine the nubmer of jobs
  if (length(degeneracyrange) == 1 & is.numeric(degeneracyrange)) {
    degeneracyrange <- c(degeneracyrange)
  }

  drange    <- seq(degeneracyrange)
  tempfiles  <- lapply(drange, function(x) {tempfile()})

  #create the trimmed file for DEGEPRIMER input
  #write alignments to disk and trim the alignment
  alignfile   <- tempfile()
  trimmedfile <- tempfile()
  writeXStringSet( as(dgprimer@msa,"DNAStringSet"), alignfile)
  trimAlignment(infile=alignfile, outfile=trimmedfile, minoccupancy=0, refsequence = NULL, trailgap = FALSE)

  degendata <- mclapply(drange, function(x) {
    #get per-run data
    outputfile    <- tempfiles[[x]]
    maxdegeneracy <- degeneracyrange[[x]]
    #calculate degeneracies
    degePrime(alignmentfile=trimmedfile, outfile=outputfile,
              oligolength=oligolength, maxdegeneracy = maxdegeneracy,
              minimumdepth=minimumdepth, skiplength=skiplength, number_iterations=number_iterations)
    #load the file and return a dataframe
    df <- read.table(outputfile,header = TRUE, stringsAsFactors = FALSE)
    df$degeneracy <- maxdegeneracy
    return(df)
  }, mc.cores=ncpus)

  #aggregate the data
  aggdata <- Reduce(rbind, degendata)

  #add coverage information
  aggdata$coverage <- aggdata$PrimerMatching/aggdata$TotalSeq


  # add the primer data to the object
  dgprimer@primerdata <- new("primerdata", aggdata)

  return(dgprimer)
}
#' Add Primer
#'
#' A convenience function to add a primer pair to the primerl list based on the
#' degenracy and location of forward and reverse primers.
#'
#' @importform purrr map_chr
#' @export
add_primerpair <-function(dgprimer, name, fpos, fdeg, rpos, rdeg) {

  pdf           <- dgprimer@primerdata
  existingpairs <- dgprimer@primerpairs

  #check degenracy and position values
  if (!fpos %in% unique(pdf$Pos)) stop("fpos value invalid. must be a Pos value present in the degeprimer table.")
  if (!rpos %in% unique(pdf$Pos)) stop("rpos value invalid. must be a Pos value present in the degeprimer table.")
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
  print(num_existing_primers)

  if (num_existing_primers == 0) {
    plist <- list(newprimer)
  } else {
    print("Appending Primers")
    plist <- as.list( c(existingpairs, newprimer))
  }

  # add the new list to the primer list slot
  dgprimer@primerpairs <- plist

  return(dgprimer)
}
#' choose a primer from the DEGEPRIMER output
#'
#'
#'
#' @import miniUI
#' @import ggplot2
#' @importFrom  shiny plotOutput runGadget dialogViewer fillRow checkboxGroupInput
#' @export
setGeneric("choose_primer", function(object) standardGeneric("choose_primer"))
#'
#' @export
setMethod("choose_primer", "degeprimer", function(object){choose_primer(object@primerdata)})
#'
#' @export
setMethod("choose_primer", "primerdata", function(object){
  df <- as.data.frame(object)

  ui <- miniPage(
    gadgetTitleBar("Drag to Select Primer", right = miniTitleBarButton("done", "Done", primary = TRUE)),

    miniContentPanel(
      fillRow(flex = c(1, 3),

              #degeneracy values
              checkboxGroupInput('show_degen', 'Primer degeneracy values to show:',
                                 unique(df$degeneracy), selected = max(unique(df$degeneracy))),

              # The brush="brush" argument means we can listen for
              # brush events on the plot using input$brush.
              plotOutput("plot", height = "100%", brush = "brush")
      )))


  server <- function(input, output, session) {
    # Render the plot
    output$plot <- renderPlot({
      # Plot the data with x/y vars indicated by the caller.
      ggplot(df[df$degeneracy %in% input$show_degen,] , aes(Pos, coverage, color=PrimerDeg)) + geom_point()
    })

    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(data, input$brush))
    })
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(df, input$brush))
    })
  }
  runGadget(ui,
            server,
            viewer =  dialogViewer("Pick Points",
                                   width = 1000,
                                   height = 700))

})
#' Extract Amplicons from Sequences
#'
#' Find matches to degnerate primers against a nucleotide sequence and
#' return the subsetted sequence or NULL
#'
#' @importFrom  Biostrings DNAStringSet DNAString matchPattern
#' @export
setGeneric("extract_amplicons", function(object, fp, rp, drop.multiple=TRUE, max.mismatch = 2) standardGeneric("extract_amplicons"))
#'
#' @importFrom  Biostrings DNAStringSet DNAString matchPattern
#' @export
setMethod("extract_amplicons", "DNAString", function(object, fp, rp, drop.multiple=TRUE, max.mismatch = 2){
  fmatches <- matchPattern(fp,     object, fixed=FALSE, max.mismatch = max.mismatch)
  rmatches <- matchPattern(rc(rp), object, fixed=FALSE, max.mismatch = max.mismatch)

  # Return NULL if there are No Matches
  if (length(fmatches) == 0 || length(rmatches) == 0) return(NULL)

  # drop the match or provide a warning if there are multiple matches
  if (length(fmatches) > 1 & drop.multiple == TRUE) return(NULL)
  if (length(fmatches) > 1)  warning("More than one match for the forward primer, using the first")
  if (length(rmatches) > 1 & drop.multiple == TRUE) return(NULL)
  if (length(rmatches) > 1)  warning("More than one match for the reverse primer, using the first")

  # check for negative indicies
  ampliconlength <- end(rmatches) - start(fmatches)
  if (ampliconlength <= 0) return(NULL)

  return(object[start(fmatches):end(rmatches)])
})
#'
#' @export
#' @importFrom  Biostrings DNAStringSet DNAString matchPattern
#' @importFrom purrr compact
setMethod("extract_amplicons", "DNAStringSet", function(object, fp, rp, drop.multiple=TRUE,max.mismatch = 2){
  amplicons <- lapply(object, function(x) {extract_amplicons(x, fp=fp, rp=rp)})
  amplicons <- compact(amplicons)
  DNAStringSet(amplicons)
})
#' Create Paired-End Reference Sequences from DNAStringSets
#'
#' When working with Amplicons that are too long to merge together it is sometiems desireable to concatenate
#' the forward and reverse reads together. In order to calssify these sequences using Blast, you would need
#' to generate referecnce sequences of the same "shape" as the amplicons. This function takes a DNAStringSet
#' and two primer sequences and extracts the amplicon but concatenated the ends together, mimickign how a paired
#' end read looks.
#'
#' @param dnastringset (Required).  A \code{\link[Biostrings]{DNAStringSet}} containing your target sequences.
#' @param fp (Required).  Default \code{NULL}. Character string of the forward primer sequence.
#' @param rp (Required).  Default \code{NULL}. Character string of the reverse primer sequence.
#' @param trimf (optional). Default \code{240}. Amount to be trimmed from the forward primer match
#' @param trimr (optional). Default \code{175}. Amount to be trimmed from the reverse primer match
#' @param drop.multiple (optional). Default \code{TRUE}. Logical.
#'   If there is more than one match to the forward or reverse sequence, should the result be dropped.
#' @param max.mismatch   (optional). Default\code{2}. Maximum allowed mismatch betwen primer sequences and the target.
#'
#'
#' @importFrom Biostrings xscat subseq
#' @export
extract_trimmed_amplicons <- function(dnastringset, fp, rp, trimf=240, trimr=175, drop.multiple=TRUE, max.mismatch = 2) {
  #get the amplicons
  amplicons     <- extract_amplicons(dnastringset, fp, rp)
  ampliconnames <- names(amplicons)

  #get forward and reverse and concatenate together
  amF <- subseq(amplicons, end=trimf)
  amR <- subseq(amplicons, start=-trimr)
  ampliconends <- xscat(amF, amR)
  names(ampliconends) <- ampliconnames

  #filter minimum size
  ampliconends[width(ampliconends) >= trimf+trimr]
}
