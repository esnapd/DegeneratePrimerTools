#' Primer validation and plot functions
#'
#' Library of functions to run Primer validation analysis.
#' 
#' @import purrr
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#' 
#' @param file Required. File location (string)
#' @param PFAM Required. PFAM id of the target for the file being filtered (string)
#' @param Ref. Required. Reference sequences for the target investigated. Either file location (string) or DNAstringSet.
#' @param logfile. Required. A file location (string)
#' @param SampleSize. Optional. Size to sample, default is 100000. (integer)
#' 
#' @export
primer_filtering <- function(file, Ref, logfile, SampleSize=100000, setwd=NULL){

  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  ###################################################################################################################################
  # From the regex and the file name, extract sample name:
  sampleName <-unlist(strsplit(basename(file), "\\."))[1]
  sampleName <- paste(unlist(strsplit(sampleName, "_"))[1:3], collapse="_")
  
  
  ###################################################################################################################################
  # Add new step in log file
  cat(paste(Sys.time(), "- Reading sequences for sample", sampleName, "\n"), file=logfile, sep="", append=TRUE)
  message(paste0("Reading ", file))
  
  
  ###################################################################################################################################
  #Read the sequences into a DNAStringSet
  Seqs <- readDNAStringSet(file)
  names(Seqs) <- gsub(" ", "_", names(Seqs))
  
  
  ###################################################################################################################################
  # Check that there are enough sequences to perform analysis (i.e. >100000)
  if (length(Seqs) < SampleSize){
    stop("Number of sequences must be over 100,000 to run the analysis.")
  }
  
  ###################################################################################################################################
  # BLAST sequences to filter out unrelated amplicons to the PFAM of interest
  cat(paste(Sys.time(), "Running BLAST to filter out unrelated sequences to the Reference set for the target", "\n"), file=logfile, sep="", append=TRUE)
  cat(paste(Sys.time(), "Blast clusters on Reference blast database\n"), file=logfile, sep="", append=TRUE)
  
  FilteredNames <- unique(sort(run_blast(Seqs, Ref, blast_args = "-evalue 1e-10")$queryID)) # Save the query sequence name that blast hit at e-10 eValue the Target sequences.
  SampledNames <- sample(FilteredNames, size = SampleSize)
  
  
  ###################################################################################################################################
  # Filter out the sequences from the original sorted file.
  cat(paste(Sys.time(), "Filter blast results to remove unrelated sequences\n"), file=logfile, sep="", append=TRUE)
  Filteredseqs <- Seqs[(names(Seqs)) %in% FilteredNames]
  FilteredFile <- paste0(wd, "/", sampleName, "_filtered.fasta")
  writeXStringSet(Filteredseqs, FilteredFile, format = "fasta")
  
  
  ###################################################################################################################################
  # Filter out the sequences from the original sorted file and subsample at size required.
  cat(paste(Sys.time(), "Filter blast results to remove unrelated sequences\n"), file=logfile, sep="", append=TRUE)
  Sampledseqs <- Seqs[(names(Seqs)) %in% SampledNames]
  SampledFile <- paste0(wd, "/", sampleName, "_sampled.fasta")
  Biostrings::writeXStringSet(Sampledseqs, SampledFile, format = "fasta")
  
  
  ###################################################################################################################################
  # Sort the filtered sequences
  SortedFile <- paste0(wd, "/", sampleName, "_sorted.fasta")
  run_vsearch_sort(SampledFile, SortedFile, MinSize=1, logfile)
  
  return(SortedFile)
}


#' Primer analysis function to validate a set of primers for a target and plot their rarefaction curve and annotated tree.
#' 
#' @param target Required. target name (string)
#' @param fileList Required. List of file locations (vector of strings)
#' @param Ref Required. Reference sequences for the target investigated. Either file location (string) or DNAstringSet or Degeprimer RDS object.
#' @param Primers. Required. Primer sequences. Either file location (string) or DNAstringSet.
#' @param OTUsizeFilter. Required. Size of the OTUs to filter out (integer).
#' 
#' @import ggplot2
#' @importFrom Biostrings writeXStringSet
#' @importFrom ggtree ggtree
#' @importFrom ggtree gheatmap
#' @importFrom gridExtra grid.arrange
#' @export
primer_analysis <- function(Target, fileList, Ref, Primers, OTUsizeFilter = 3){
  
  logpath <-getwd()
  wd <- tempdir()
  on.exit(unlink(list.files(wd)))
  
  ###################################################################################################################################
  # Check log file.
  
  if(!is.null(setwd)){logfile <- paste(logpath, "/log.txt", sep="")} else {logfile <- "log.txt"}
  
  df85 <- NULL
  df90 <- NULL
  Cluster90AllSeq <- NULL
  Cluster85AllSeq <- NULL
  SummaryTable <- NULL

  for (i in 1:length(fileList)) {
      
      # Read the filename and extract the sample name, assuming a field separator of "_"
      # Need to add a check on that field.
      sampleName <-unlist(strsplit(basename(fileList[i]), "\\."))[1]
      sampleName <- paste(unlist(strsplit(sampleName, "_"))[1:3], collapse="_")
      
      FilteredSeqsFasta <- primer_filtering(file = fileList[i], Ref, SampleSize=100000, logfile)
      
      ###################################################################################################################################
      # Cluster filtered sequences at 90% identity
      
      OutFasta90 <- paste0(wd, "/", sampleName, "_clustered90.fasta")
      UCFile90 <- paste0(wd, "/", sampleName, "_clustered90.uc")
      
      run_vsearch_cluster(FilteredSeqsFasta, id=0.90, OutFasta90, UCFile90, logfile)
      Cluster90Seqs <- readDNAStringSet(OutFasta90)
      UC90phylo <- filter_taxa(import_usearch_uc(UCFile90), function(x) x>= OTUsizeFilter, TRUE)
      
      
      ###################################################################################################################################
      # Cluster filtered sequences at 85% identity
      
      OutFasta85 <- paste0(wd, "/", sampleName, "_clustered85.fasta")
      UCFile85 <- paste0(wd, "/", sampleName, "_clustered85.uc")
      
      run_vsearch_cluster(FilteredSeqsFasta, id=0.85, OutFasta85, UCFile85, logfile)
      Cluster85Seqs <- readDNAStringSet(OutFasta85)
      UC85phylo <- filter_taxa(import_usearch_uc(UCFile85), function(x) x>= OTUsizeFilter, TRUE)
      
      
      ###################################################################################################################################
      # phyloseq-tools rarefy for 90% cluster identity
      
      rarefaction_curve_data90 <- (calculate_rarefaction_curves(UC90phylo, c('Observed'), seq(1, max(sample_sums(UC90phylo)), length.out=100), parallel = FALSE))[c(1,4)]
      rarefaction_curve_data90 <- cbind(rarefaction_curve_data90, rep(as.character(sampleName),100))
      rarefaction_curve_data90 <- cbind(rarefaction_curve_data90, rep(as.character(Target),100))
      rarefaction_curve_data90 <- cbind(seq(1:100), rarefaction_curve_data90)
      
      colnames(rarefaction_curve_data90) <- c("idx","Depth","Diversity", "Sample", "Target")
      
      df90 <- rbind (df90, rarefaction_curve_data90)
      Cluster90AllSeq <- append(Cluster90AllSeq, Cluster90Seqs, after=length(Cluster90AllSeq))
      
      
      ###################################################################################################################################
      # phyloseq-tools rarefy for 85% cluster identity
      
      rarefaction_curve_data85 <- (calculate_rarefaction_curves(UC85phylo, c('Observed'), seq(1, max(sample_sums(UC85phylo)), length.out=100), parallel = FALSE))[c(1,4)]
      rarefaction_curve_data85 <- cbind(rarefaction_curve_data85, rep(as.character(sampleName),100))
      rarefaction_curve_data85 <- cbind(rarefaction_curve_data85, rep(as.character(Target),100))
      rarefaction_curve_data85 <- cbind(seq(1:100), rarefaction_curve_data85)
      
      colnames(rarefaction_curve_data85) <- c("idx","Depth","Diversity", "Sample", "Target")
      
      df85 <- rbind (df85, rarefaction_curve_data85)
      Cluster85AllSeq <- append(Cluster85AllSeq, Cluster85Seqs, after=length(Cluster85AllSeq))
      
      
      ###################################################################################################################################
      # Create summary table:
      
      NbTotalSeqs <- length(readDNAStringSet(fileList[i]))
      NbFilteredSeqs <- length(readDNAStringSet(paste0(wd,"/",sampleName,"_filtered.fasta")))
      NbSampledSeqs <- length(readDNAStringSet(paste0(wd,"/",sampleName,"_sampled.fasta")))
      Nb85OTUs <- length(which(as.numeric(gsub('.*=(.*);', '\\1', names(Cluster85Seqs))) >= OTUsizeFilter))
      Nb90OTUs <- length(which(as.numeric(gsub('.*=(.*);', '\\1', names(Cluster90Seqs))) >= OTUsizeFilter))
        
      SummaryTable = rbind(SummaryTable, data.frame(
                                            sampleName,
                                            NbTotalSeqs,
                                            paste(NbFilteredSeqs," (~", round(NbFilteredSeqs/NbTotalSeqs*100, 2), "%)", sep = ""),
                                            length(readDNAStringSet(paste0(wd,"/",sampleName,"_sampled.fasta"))),
                                            Nb85OTUs,
                                            Nb90OTUs))
      
  }
  
  
  Samples85 <- levels(df85$Sample[which(df85$Target %in% Target)])
  Samples90 <- levels(df90$Sample[which(df90$Target %in% Target)])

  
  grouping85 <- rep(Samples85, each = 100)
  grouping90 <- rep(Samples90, each = 100)

      RarefPlot85 <- ggplot(df85, aes(x=Depth,y=Diversity, group=grouping85, colour=grouping85)) + 
        geom_point(shape = 1, size = 0.6) + 
        geom_smooth(span = 1) +
        guides(fill=guide_legend(title=NULL)) +
        facet_wrap(~Target, scales="free")
      RarefPlot90 <- ggplot(df90, aes(x=Depth,y=Diversity, group=grouping90, colour=grouping90)) + 
        geom_point(shape = 1, size = 0.6) + 
        geom_smooth(span = 1) +
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
      
      # Get the predicted amplicons from the primer pair sequences for the Reference sequences with the extract_amplicons function
    
      # Uses matchpattern from the BioString package.
      #
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
      RefTrim <- extract_ends(extract_amplicons(Ref, "TTCCAGAACGGMAARCTGATG", "GGAACATSGCCATGTAGTGSG"))
      if (length(RefTrim) == 0) {
        warning("No reference sequence can be extracted with the primer data.")
      }
      
      # CONCATENATE WITH THE REFERENCE...
      ConcatMSA85 <- c(extract_ends(Cluster85AllSeq), RefTrim)
      ConcatMSA90 <- c(extract_ends(Cluster90AllSeq), RefTrim)
      
      # Then run the multiple sequence alignment with muscle, calling the run_muscle function.
      MSAClustered85 <- NULL
      MSAClustered85 <- run_muscle(ConcatMSA85)
      MSAClustered90 <- NULL
      MSAClustered90 <- run_muscle(ConcatMSA90)
      
      ###################################################################################################################################
      # Generation of the tree with FastTree
      
      # check if the names of the sequences contain ":" or ";" and switch to "_" if any:
      names(MSAClustered85) <- gsub("\\:", "_", names(MSAClustered85))
      names(MSAClustered85) <- gsub("\\;", "_", names(MSAClustered85))
      names(MSAClustered90) <- gsub("\\:", "_", names(MSAClustered90))
      names(MSAClustered90) <- gsub("\\;", "_", names(MSAClustered90))
      
      # Then call the FastTree function on the multiple sequence alignment results from Muscle
      Tree85 <- run_fasttree(MSAClustered85)
      Tree90 <- run_fasttree(MSAClustered90)
      
      ###################################################################################################################################
      # Now plot tree and rarefaction curve
      
      # Prepare the labels 
      RefIDs <- reshape2:::colsplit(string=names(Ref), pattern=" ", names=c("1","2"))[1]
      
      LABS85 <- strsplit(Tree85$tip.label, "_M038") # How to cut this?
      LABELS85 <- do.call(rbind, LABS85)[,1]
      LABS90 <- strsplit(Tree90$tip.label, "_M038") # How to cut this?
      LABELS90 <- do.call(rbind, LABS90)[,1]
      
      SampleMatrix85 <- NULL
      SampleMatrix85 <- do.call("cbind", list(LABELS85))
      SampleMatrix85 <- cbind(Tree85$tip.label, SampleMatrix85, SampleMatrix85)
      SampleMatrix90 <- NULL
      SampleMatrix90 <- do.call("cbind", list(LABELS90))
      SampleMatrix90 <- cbind(Tree90$tip.label, SampleMatrix90, SampleMatrix90)
      
      SampleMatrix85[which(SampleMatrix85[,2] %in% unlist(RefIDs)),2] <- ""
      SampleMatrix85[which(!SampleMatrix85[,3] %in% unlist(RefIDs)),3] <- ""
      SampleMatrix90[which(SampleMatrix90[,2] %in% unlist(RefIDs)),2] <- ""
      SampleMatrix90[which(!SampleMatrix90[,3] %in% unlist(RefIDs)),3] <- ""
      rownames(SampleMatrix85) <- Tree85$tip.label
      colnames(SampleMatrix85) <- c("TipLab","Samples", "Reference")
      rownames(SampleMatrix90) <- Tree90$tip.label
      colnames(SampleMatrix90) <- c("TipLab","Samples", "Reference")
      dd85 <- data.frame(SampleMatrix85)
      SampleMatrix85[which(SampleMatrix85[,3] %in% unlist(RefIDs)),3] <- "Reference"
      dd90 <- data.frame(SampleMatrix90)
      SampleMatrix90[which(SampleMatrix90[,3] %in% unlist(RefIDs)),3] <- "Reference"
      
      
      
      
      #Tree$tip.label <- TipLabels
      TREE85.laz<-ladderize(Tree85,TRUE)
      TREE90.laz<-ladderize(Tree90,TRUE)
      
      
      # Plot tree
      T85 <- ggtree(TREE85.laz) #, branch.length="none") +
        T85 <- T85 %<+% dd85 + # geom_tiplab(aes(label=Reference, subset=isTip), size =2, color="darkblue") +
          geom_text_repel(aes(label = Reference),size = 3, force = 1, nudge_x = 0.4, nudge_y=1, color="red")
          
      # Create the heatmap displaying the different 
      H85 <- gheatmap(T85, SampleMatrix85[,-1], offset = 0.5, width=0.2, colnames = FALSE) # + scale_fill_manual(values=ColourPalette)
      # Plot tree
      T90 <- ggtree(TREE90.laz) #, branch.length="none") +
        T90 <- T90 %<+% dd90 + #geom_tiplab(aes(label=Reference, subset=isTip), size =2, color="darkblue")
        geom_text_repel(aes(label = Reference),size = 3, force = 1, nudge_x = 0.4, nudge_y=1, color="red")
        
      # Create the heatmap displaying the different 
      H90 <- gheatmap(T90, SampleMatrix90[,-1], offset = 0.5, width=0.2, colnames = FALSE) # + scale_fill_manual(values=ColourPalette)
      
      # Organize both Rarefaction curve and tree with the sample annotation Heatmap side by side with grid.arrange.
      #grid.arrange(RarefPlot85, H85, ncol = 2)
      #grid.arrange(RarefPlot90, H90, ncol = 2)
      
      # Create a table plot
      names(SummaryTable) <- c("Sample Names", "Total paired reads", "Filtered", "Sampling", paste("85% clusters at size>=",OTUsizeFilter), paste("90% clusters at size>=",OTUsizeFilter))
      
      # Set theme to allow for plotmath expressions
      tbl85 <- tableGrob(SummaryTable[,-6], rows=NULL)
      tbl90 <- tableGrob(SummaryTable[,-5], rows=NULL)
      
      # Plot chart and table into one object
      grid.arrange(tbl85, RarefPlot85, H85,
                   nrow=2,
                   as.table=TRUE,
                   heights=c(1,3), layout_matrix = rbind(c(1,1), c(2,3)))
      # Plot chart and table into one object
      grid.arrange(tbl90, RarefPlot90, H90,
                   nrow=2,
                   as.table=TRUE,
                   heights=c(1,3), layout_matrix = rbind(c(1,1), c(2,3)))
      # Plot chart and table into one object
      grid.arrange(tbl85, RarefPlot85, H85,tbl90, RarefPlot90, H90,
                   nrow=4,
                   as.table=TRUE,
                   heights=c(1,3,1,3), layout_matrix = rbind(c(1,1), c(2,3),c(4,4), c(5,6)))
      
}
  