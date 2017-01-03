#' Primer validation and plot functions
#'
#' Library of functions to run Primer validation analysis.
#' 
#' @import rBLAST
#' @import phyloseq
#' @import dplyr
#' @import ggrepel
#' @import ggplot2
#' @import purrr
#' @import ape
#' @import Biobase
#' @import gridExtra
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' 
#' @examples
#' PFAM <- "PF13714"
#' 
#' @param file. Required. File location (string)
#' @param Ref. Required. Reference sequences for the target investigated. Either file location (string) or DNAstringSet.
#' @param logfile. Required. A file location (string)
#' 
#' @export
primer_filtering <- function(file, Ref, logfile, setwd=NULL){

  wd <- tempdir()
  on.exit(unlink(list.files(wd)))

  ###################################################################################################################################
  # From the regex and the file name, extract sample name:
  sampleName <- paste(unlist(strsplit(basename(file), "_"))[1:3], collapse="_")
  
  ###################################################################################################################################
  # Add new step in log file
  cat(paste(Sys.time(), "- Initiate Clustering with Vsearch for sample", sampleName, "\n"), file=logfile, sep="", append=TRUE)
  message(paste0("Reading ", file))

  ###################################################################################################################################
  # Sort dereplicated sequences
  SortedFile <- paste0(wd, "/", sampleName, "_sorted.fasta")
  #run_vsearch_sort(DerepFile, SortedFile, 3, logfile)
  
  ###################################################################################################################################
  # Cluster sequences at 97% identity
  #OutFasta97 <- paste0(wd, "/", sampleName, "_clustered97.fasta")
  #UCFile <- paste0(wd, "/", sampleName, "_clustered97.uc")
  #run_vsearch_cluster(SortedFile, id=0.97, OutFasta97, UCFile, logfile)
  
  ###################################################################################################################################
  #Read the sorted sequences and the clustered sequences at 97% into a DNAStringSet
  #SortedSeqs <- readDNAStringSet(SortedFile)
  SortedSeqs <- readDNAStringSet(file)
  #Cluster97Seqs <- readDNAStringSet(OutFasta97)

  ###################################################################################################################################
  # BLAST sequences to filter out unrelated amplicons to the PFAM of interest
  cat(paste(Sys.time(), "Running BLAST to filter out unrelated sequences to the Reference set for the target", "\n"), file=logfile, sep="", append=TRUE)
  cat(paste(Sys.time(), "Blast clusters on Reference blast database\n"), file=logfile, sep="", append=TRUE)
  
  FilteredNames <- sample(unique(sort(run_blast(SortedSeqs, Ref, blast_args = "-evalue 1e-10")$queryID)),size = 30000) # Save the query sequence name that blast hit at e-10 eValue the Target sequences.
  
  ###################################################################################################################################
  # Filter out the sequences from the original sorted file.
  cat(paste(Sys.time(), "Filter blast results to remove unrelated sequences\n"), file=logfile, sep="", append=TRUE)
  Filteredseqs <- SortedSeqs[(names(SortedSeqs)) %in% FilteredNames]
  FilteredFile <- paste0(wd, "/", sampleName, "_filtered.fasta")
  Biostrings::writeXStringSet(Filteredseqs, FilteredFile, format = "fasta")
  
  
  ###################################################################################################################################
  # Dereplicate with vsearch
  DerepFile <- paste0(wd, "/", sampleName, "_filtered_derep.fasta")
  #run_vsearch_derep(FilteredFile, DerepFile, logfile)
  
  ###################################################################################################################################
  # Sort the filtered sequences
  SortedFile97 <- paste0(wd, "/", sampleName, "_sorted97.fasta")
  run_vsearch_sort(FilteredFile, SortedFile97, MinSize =1, logfile)
  
  return(SortedFile97)
}







#' Primer analysis function to validate a set of primers for a target and plot their rarefaction curve and annotated tree.
#' 
#' @param Target. Required. Target name (string)
#' @param fileList. Required. List of file locations (vector of strings)
#' @param Ref. Required. Reference sequences for the target investigated. Either file location (string) or DNAstringSet.
#' @param Primers. Required. Primer sequences. Either file location (string) or DNAstringSet.
#' @param OTUsizeFilter. Required. Size of the OTUs to filter out (integer).
#' 
#' @export
primer_analysis <- function(Target, fileList, Ref, Primers, OTUsizeFilter = 3){
  
  logpath <-getwd()
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  ###################################################################################################################################
  # Check log file.
  
  if(!is.null(setwd)){logfile <- paste(logpath, "/log.txt", sep="")} else {logfile <- "log.txt"}
  
  df <- NULL
  Cluster95AllSeq <- NULL

  for (i in 1:length(fileList)) {
      
      # Read the filename and extract the sample name, assuming a field separator of "_"
      # Need to add a check on that field.
      sampleName <- paste(unlist(strsplit(basename(fileList[i]), "_"))[1:3], collapse="_")
  
      FilteredSeqsFasta <- primer_filtering(file = fileList[i], Ref, logfile)
      
      
      ###################################################################################################################################
      # Cluster filtered sequences at 95% identity
      
      OutFasta95 <- paste0(wd, "/", sampleName, "_clustered95.fasta")
      UCFile95 <- paste0(wd, "/", sampleName, "_clustered95.uc")
      
      run_vsearch_cluster(FilteredSeqsFasta, id=0.95, OutFasta95, UCFile95, logfile)
      Cluster95Seqs <- readDNAStringSet(OutFasta95)
      UC95phylo <- filter_taxa(import_usearch_uc(UCFile95), function(x) x>= OTUsizeFilter, TRUE)
      
      
      ###################################################################################################################################
      # phyloseq-tools rarefy
      
      rarefaction_curve_data <- (calculate_rarefaction_curves(UC95phylo, c('Observed'), seq(1, max(sample_sums(UC95phylo)), length.out=100), parallel = FALSE))[c(1,4)]
      #rarefaction_curve_data$Alpha_diversity_mean <- lowess(rarefaction_curve_data$Depth, rarefaction_curve_data$Alpha_diversity_mean, iter = 10000)$y
      rarefaction_curve_data <- cbind(rarefaction_curve_data, rep(as.character(sampleName),100))
      rarefaction_curve_data <- cbind(rarefaction_curve_data, rep(as.character(Target),100))
      rarefaction_curve_data <- cbind(seq(1:100), rarefaction_curve_data)
      
      colnames(rarefaction_curve_data) <- c("idx","Depth","Diversity", "Sample", "Target")
      
      #summary(rarefaction_curve_data)
      
      df <- rbind (df, rarefaction_curve_data)
      Cluster95AllSeq <- append(Cluster95AllSeq, Cluster95Seqs, after=length(Cluster95AllSeq))

  }
  
  
  
  
  ### Now For each Target in the entire set we want to plot their respective rarefaction curve and phylo tree 
  #targets <- unique(df$Target)

  #for (i in 1:length(targets)){
  
  #Target <- targets[i]
  Samples <- levels(df$Sample[which(df$Target %in% Target)])
  
  # Prepare the colours for each individual sample within the Target.
  #numbertiplabs<-length(Tree$tip.label)
  #colourtips <- rep("black",numbertiplabs)
  #for (sam in 1:length(Samples)){
  #  colourtips[grep(unlist(strsplit(Samples[sam], "_"))[2], Tree$tip.label)] <- ColourPalette[sam]
  #}
  #colourtips[grep("PEP_Control", Tree$tip.label)] <- ColourPalette[which(Samples=="PEP_Control")]
  #plot(TREE.laz,tip.color=colourtips,adj=1,cex=0.3,use.edge.length=F)
  
  #legend('topleft', Samples, lty=1, col=ColourPalette, bty='n', cex=0.8)
      
      grouping <- rep(Samples, each = 100)

      RarefPlot <- ggplot(df, aes(x=Depth,y=Diversity, group=grouping, colour=grouping)) + 
        geom_point(shape = 1, size = 0.6) + 
        geom_smooth(span = 1) +
        #geom_line() +
        #scale_color_manual(values=ColourPalette) +
        #geom_text_repel(data = maxperprimer,aes(label=Samples)) +
        guides(fill=guide_legend(title=NULL)) +
        facet_wrap(~Target, scales="free")


      ###################################################################################################################################
      # Multiple sequence alignment with Muscle to generate the tree
      # CONCATENATE WITH THE REFERENCE...
      if (is.character(Ref)){
        Ref <- readDNAStringSet(Ref)
      }
      #if (is.character(Primers)){
      #  Primers <- readDNAStringSet(Primers)
      #}
      
      # Get the predicted amplicons from the primer pair sequences for the Reference sequences with the extract_amplicons function
      # Uses matchpattern from the BioString package.
      
      #for (i in 1:length(RefObj@primerdata[,7])){
      #  
      #  PrimerMatch <- vmatchPattern(RefObj@primerdata[i,7], subject = Ref)
      #  if ((!is.null(PrimerMatch@ends)) & (start > (min(unlist(PrimerMatch@ends))+PrimerMatch@width0))){
      #    start <- min(unlist(PrimerMatch@ends))
      #  }
      #  if ((!is.null(PrimerMatch@ends)) & (end < max(unlist(PrimerMatch@ends)))){
      #    end <- max(unlist(PrimerMatch@ends))
      #  }
     # }
      
      # Trim the set of reference sequences with the extreme start and end positions
      RefTrim <- extract_ends(extract_amplicons(Ref, "STGCGGGTGCTGCCSGACGAC", "SGCGTASAGGTACTGCAGC"))
      
      # CONCATENATE WITH THE REFERENCE...
      ConcatMSA <- c(extract_ends(Cluster95AllSeq), RefTrim)
      # Then run the multiple sequence alignment with muscle, calling the run_muscle function.
      MSAClustered95 <- NULL
      MSAClustered95 <- run_muscle(ConcatMSA)
      
      ###################################################################################################################################
      # Generation of the tree with FastTree
      
      # check if the names of the sequences contain ":" or ";" and switch to "_" if any:
      names(MSAClustered95) <- gsub("\\:", "_", names(MSAClustered95))
      names(MSAClustered95) <- gsub("\\;", "_", names(MSAClustered95))
      
      # Then call the FastTree function on the multiple sequence alignment results from Muscle
      Tree <- run_fasttree(MSAClustered95)
      
      ###################################################################################################################################
      # Now plot tree and rarefaction curve
      
      # Prepare the labels 
      RefIDs <- colsplit(string=names(Ref), pattern=" ", names=c("1","2"))[1]
      
      LABS <- strsplit(Tree$tip.label, "_M038")
      LABELS <- do.call(rbind, LABS)[,1]
      
      SampleMatrix <- NULL
      SampleMatrix <- do.call("cbind", list(LABELS))
      SampleMatrix <- cbind(Tree$tip.label, SampleMatrix, SampleMatrix)
      
      SampleMatrix[which(SampleMatrix[,2] %in% unlist(RefIDs)),2] <- ""
      SampleMatrix[which(!SampleMatrix[,3] %in% unlist(RefIDs)),3] <- ""
      rownames(SampleMatrix) <- Tree$tip.label
      colnames(SampleMatrix) <- c("TipLab","Samples", "Reference")
      dd <- data.frame(SampleMatrix)
      SampleMatrix[which(SampleMatrix[,3] %in% unlist(RefIDs)),3] <- "Reference"
      
      
      
      
      #Tree$tip.label <- TipLabels
      TREE.laz<-ladderize(Tree,TRUE)
      
      
      # Plot tree
      T <- ggtree(TREE.laz) #, branch.length="none") +
            T <- T %<+% dd + geom_tiplab(aes(label=Reference, subset=isTip), size =2, color="darkblue")
      # Create the heatmap displaying the different 
      H <- gheatmap(T, SampleMatrix[,-1], offset = 0.1, width=0.2, colnames = FALSE) # + scale_fill_manual(values=ColourPalette)
      
      # Organize both Rarefaction curve and tree with the sample annotation Heatmap side by side with grid.arrange.
      grid.arrange(RarefPlot, H, ncol = 2)
      
  #}
}




  
  
  