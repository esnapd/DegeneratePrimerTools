<<<<<<< 4179cbd6e2e43077fe998d8cd427cf6f71ec1c24
#' primer_validation
=======
#' Primer validation and plot functions
>>>>>>> Large improvements of the code for the primer validation functions
#'
#' Library of functions to run Primer validation analysis.
#'
<<<<<<< 4179cbd6e2e43077fe998d8cd427cf6f71ec1c24
#' @import phyloseq
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' #library(rBLAST)
primer_validation <- function(file, vsearchpath="integrated", PFAM, cmd="", setwd=NULL){

  logpath <-getwd()
  wd <- tempdir()
  #setwd("/Users/christophelemetre/Documents/Work/eSNaPD3/")
  #on.exit(unlink(list.files(wd)))
=======
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
 

logpath <-getwd()
wd <- tempdir()
on.exit(unlink(list.files(wd)))


primer_filtering <- function(file, vsearchpath, PFAM, cmd="", setwd=NULL){
  
>>>>>>> Large improvements of the code for the primer validation functions
  
  ###### Preparing the PFAM sequences to blast against for filtering
  #PFAM <- "PF13714" # for PEP
  #PFAMnucSeqs <- retrieve_PFAM_nucleotide_sequences(PFAM, alignmenttype = "full")
  
  #PFAMseqs <- DNAStringSet(PFAMnucSeqs$domainsequence[which(!is.na(PFAMnucSeqs$domainsequence))],use.names = TRUE)
  #names(PFAMseqs) <- paste(sep="_",PFAMnucSeqs$PFAM_ID[which(!is.na(PFAMnucSeqs$domainsequence))],PFAMnucSeqs$Accession[which(!is.na(PFAMnucSeqs$domainsequence))],PFAMnucSeqs$EMBL_ID[which(!is.na(PFAMnucSeqs$domainsequence))])
  #message("Done with PFAM file")
  
  
  vsearchCmd <- paste0(vsearchpath,"vsearch")
  
  ###################################################################################################################################
  # Check log file.
  
  if(!is.null(setwd)){logfile <- paste(logpath, "/log.txt", sep="")} else {logfile <- "log.txt"}
  if(!file.exists(file)){meep <- paste(logpath, "/", sep="")} else{meep <- ""}
  
  
  ###################################################################################################################################
  # From the regex and the file name, extract target and sample name:
  
  filename <- gsub("(.*).fasta$", "\\1", basename(file))
  sampleName <- paste(unlist(strsplit(basename(file), "_"))[1:2], collapse="_")
  target <- unlist(strsplit(basename(file), "_"))[1]
  
  
  ###################################################################################################################################
  # Add new step in log file
  
  cat(paste(Sys.time(), "- Initiate Clustering with Vsearch for sample", sampleName, "\n"), file=logfile, sep="", append=TRUE)
  message(paste0("Reading ", file))
  
  
  ###################################################################################################################################
  # Dereplicate with vsearch
  
  DerepFile <- paste0(wd, "/", sampleName, "_combined_derep.fasta")
  run_vsearch_derep(vsearchCmd, file, DerepFile, logfile)

  
  ###################################################################################################################################
  # Sort dereplicated sequences
  
  SortedFile <- paste0(wd, "/", sampleName, "_sorted.fasta")
  run_vsearch_sort(vsearchCmd, DerepFile, SortedFile, logfile)
  
  
  ###################################################################################################################################
  # Cluster sequences at 97% identity
  
  OutFasta97 <- paste0(wd, "/", sampleName, "_clustered97.fasta")
  UCFile <- paste0(wd, "/", sampleName, "_clustered97.uc")
  run_vsearch_cluster(vsearchCmd, SortedFile, id=0.97, OutFasta97, UCFile, logfile)
  
  #Read the sorted sequences and the clustered sequences at 97% into a DNAStringSet
  SortedSeqs <- readDNAStringSet(SortedFile)
  Cluster97Seqs <- readDNAStringSet(OutFasta97)

  
  ###################################################################################################################################
  # BLAST sequences to filter out unrelated amplicons to the PFAM of interest
  
  cat(paste(Sys.time(), "Running BLAST to filter out unrelated sequences to the PFAM of interest", "\n"), file=logfile, sep="", append=TRUE)
  cat(paste(Sys.time(), "Blast clusters on PFAM blast database\n"), file=logfile, sep="", append=TRUE)
  
  FilteredNames <- unique(sort(run_blast(SortedSeqs, PFAMseqs, blast_args = "-evalue 1e-10")$queryID)) # Save the query sequence name that blast hit at e-10 eValue the Target sequences.
  
  
  ###################################################################################################################################
  # Filter out the sequences from the original sorted file.
  
  cat(paste(Sys.time(), "Filter blast results to remove unrelated sequences\n"), file=logfile, sep="", append=TRUE)
  Filteredseqs <- SortedSeqs[(names(SortedSeqs)) %in% FilteredNames]
  FilteredFile <- paste0(wd, "/", sampleName, "_filtered.fasta")
  Biostrings::writeXStringSet(Filteredseqs, FilteredFile, format = "fasta")
  
  
  ###################################################################################################################################
  # Sort the filtered sequences
  
  SortedFile97 <- paste0(wd, "/", sampleName, "_sorted97.fasta")
  run_vsearch_sort(vsearchCmd, FilteredFile, SortedFile97, logfile)
  
  
  ###################################################################################################################################
  # Clean up
<<<<<<< 4179cbd6e2e43077fe998d8cd427cf6f71ec1c24
  # delete all the intermediate files.
  

  
}

  
  
  
  
  
  
  
  
=======
  # delete all the intermediate files and/or variables..
  
  
  
  return(SortedFile97)
}


PFAMlist <- c("PF13714", "PF13714", "PF13714", "PF13714", "PF13714", "PF13714")
fileList <- c("/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_01_combined.fasta",
              "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_02_combined.fasta",
              "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_04_combined.fasta",
              "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_06_combined.fasta",
              "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_Control_combined.fasta")
df <- NULL
Cluster95AllSeq <- NULL

for (i in 1:length(fileList)) {
      
      # Read the filename and extract the sample name and its target, assuming a field separator of "_"
      # Need to add a check on that field.
      filename <- gsub("(.*).fq$", "\\1", basename(fileList[i]))
      sampleName <- paste(unlist(strsplit(basename(fileList[i]), "_"))[1:2], collapse="_")
      target <- unlist(strsplit(basename(fileList[i]), "_"))[1]
  
  
      FilteredSeqsFasta <- primer_filtering(file = fileList[i], PFAM = PFAMlist[i], vsearchpath = "/Users/christophelemetre/miniconda2/bin/vsearch-2.3.0-osx-x86_64/bin/")
      
      
      ###################################################################################################################################
      # Cluster filtered sequences at 95% identity
      
      OutFasta95 <- paste0(wd, "/", sampleName, "_clustered95.fasta")
      UCFile95 <- paste0(wd, "/", sampleName, "_clustered95.uc")
      run_vsearch_cluster(vsearchCmd, FilteredSeqsFasta, id=0.97, OutFasta95, UCFile95, logfile)
      Cluster95Seqs <- readDNAStringSet(OutFasta95)
      UC95phylo <- import_usearch_uc(UCFile95)
      
      
      ###################################################################################################################################
      # phyloseq-tools rarefy
      
      # Change the Path
      #source("/Users/christophelemetre/Documents/Work/Primer_Design/R/phyloseq-tools.R")
      rarefaction_curve_data <- (calculate_rarefaction_curves(UC95phylo, c('Observed'), seq(1, max(sample_sums(UC95phylo)), length.out=100), parallel = FALSE))[c(1,4)]
      #rarefaction_curve_data$Alpha_diversity_mean <- lowess(rarefaction_curve_data$Depth, rarefaction_curve_data$Alpha_diversity_mean, iter = 10000)$y
      rarefaction_curve_data <- cbind(rarefaction_curve_data, rep(as.character(sampleName),100))
      rarefaction_curve_data <- cbind(rarefaction_curve_data, rep(as.character(target),100))
      rarefaction_curve_data <- cbind(seq(1:100), rarefaction_curve_data)
      
      colnames(rarefaction_curve_data) <- c("idx","Depth","Diversity", "Sample", "Target")
      
      #summary(rarefaction_curve_data)
      
      df <- rbind (df, rarefaction_curve_data)
      Cluster95AllSeq <- append(Cluster95AllSeq, Cluster95Seqs, after=length(Cluster95AllSeq))

}

maxperprimer <- df %>% 
  arrange(rank(-Diversity)) %>%
  group_by(Sample) %>%
  #slice(1:1) %>%
  ungroup()


### Now For each Target in the entire set we want to plot their respective rarefaction curve and phylo tree 
targets <- unique(df$Target)

for (i in 1:length(targets)){

      Target <- targets[i]
      Samples <- levels(df$Sample[which(df$Target %in% Target)])
      
      #Defining the colour palette from the number of samples for the Target.
      #ColourPalette <- palette(rainbow(length(Samples)))
      
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
      
      MSAClustered95 <- run_muscle(Cluster95AllSeq)


      ###################################################################################################################################
      # Generation of the tree with FastTree
      
      # check if the names of the sequences contain ":" or ";" and switch to "_" if any:
      names(MSAClustered95) <- gsub("\\:", "_", names(MSAClustered95))
      names(MSAClustered95) <- gsub("\\;", "_", names(MSAClustered95))
      
      # Then call the FastTree function on the multiple sequence alignment results from Muscle
      Tree <- run_fasttree(MSAClustered95)
      
      
      ###################################################################################################################################
      # Now plot tree and rarefaction curve
      
      # Prepare the colours for each individual sample within the Target.
      #numbertiplabs<-length(Tree$tip.label)
      #colourtips <- rep("black",numbertiplabs)
      #for (sam in 1:length(Samples)){
      #  colourtips[grep(unlist(strsplit(Samples[sam], "_"))[2], Tree$tip.label)] <- ColourPalette[sam]
      #}
      #colourtips[grep("PEP_Control", Tree$tip.label)] <- ColourPalette[which(Samples=="PEP_Control")]
      #plot(TREE.laz,tip.color=colourtips,adj=1,cex=0.3,use.edge.length=F)
      #legend('topleft', Samples, lty=1, col=ColourPalette, bty='n', cex=0.8)
      
      # Prepare the labels 
      LABS <- strsplit(Tree$tip.label, "_M038")
      LABELS <- do.call(rbind, LABS)[,1]
      SampleMatrix <- NULL
      SampleMatrix <- do.call("cbind", list(LABELS))
      rownames(SampleMatrix) <- Tree$tip.label
      colnames(SampleMatrix) <- "Samples"
      
      #Tree$tip.label <- LABELS
      TREE.laz<-ladderize(Tree,FALSE)
      
      
      # Plot tree
      T <- ggtree(TREE.laz, branch.length="none")# + geom_tippoint(color=colourtips, shape=20, size=1)
      H <- gheatmap(T, SampleMatrix, offset = 0, width=0.2, colnames = FALSE)# + scale_fill_manual(values=ColourPalette)
      
      # Organize both Rarefaction curve and tree with the sample annotation Heatmap side by side with grid.arrange.
      grid.arrange(RarefPlot, H, ncol = 2)
      
}


>>>>>>> Large improvements of the code for the primer validation functions



