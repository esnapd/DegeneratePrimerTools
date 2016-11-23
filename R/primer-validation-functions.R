#' eSNaPD core functions
#'
#' Library of functions to run eSNaPD.
#'
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' 
#' @examples
#' PFAM <- "PF13714"
 
#devtools::install_github("zachcp/phyloseq-tools")
#source("/Users/christophelemetre/Documents/Work/DegeneratePrimerTools-1/R/phyloseq-tools.R")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
library(rBLAST)



primer_validation <- function(file, vsearchpath="integrated", PFAM, cmd="", setwd=NULL){

  logpath <-getwd()
  wd <- tempdir()
  setwd("/Users/christophelemetre/Documents/Work/eSNaPD3/")
  #on.exit(unlink(list.files(wd)))
  
  ###### Preparing the PFAM sequences to blast against for filtering
  #PFAM <- "PF13714" # for PEP
  PFAMnucSeqs <- retrieve_PFAM_nucleotide_sequences(PFAM, alignmenttype = "full")
  
  PFAMseqs <- DNAStringSet(PFAMnucSeqs$domainsequence[which(!is.na(PFAMnucSeqs$domainsequence))],use.names = TRUE)
  names(PFAMseqs) <- paste(sep="_",PFAMnucSeqs$PFAM_ID[which(!is.na(PFAMnucSeqs$domainsequence))],PFAMnucSeqs$Accession[which(!is.na(PFAMnucSeqs$domainsequence))],PFAMnucSeqs$EMBL_ID[which(!is.na(PFAMnucSeqs$domainsequence))])
  message("Done with PFAM file")
  
  # To change for pre check of pre existence of datasets on the repo
  if (vsearchpath=="integrated"){
    if (operating_system=="MacOSX"){
      vsearchpath <- paste("/Users/christophelemetre/miniconda2/bin", "/vsearch-2.3.0-osx-x86_64/bin/", sep="")
      usearchpath <- paste("/Users/christophelemetre/miniconda2/bin", "/usearch8.1/", sep="")
      Sys.chmod(vsearchpath, mode = "0777", use_umask = TRUE)} else {
        vsearchpath <- paste(system.file(package="Primer_Design"), "/vsearch-2.3.0-osx-x86_64", sep="")
        Sys.chmod(vsearchpath, mode = "0777", use_umask = TRUE)
      }
  }
  
  
  ###################################################################################################################################
  # Check log file.
  
  if(!is.null(setwd)){logfile <- paste(getwd(), "/log.txt", sep="")} else {logfile <- "log.txt"}
  if(!file.exists(file)){meep <- paste(getwd(), "/", sep="")} else{meep <- ""}
  
  
  ###################################################################################################################################
  # Create the new Vsearch folder for all the vsearch cluster files.
  
  dir.create(paste(getwd(), "/Vsearch", sep=""), showWarnings=F)
  filename <- gsub("(.*).fq$", "\\1", basename(file))
  sampleName <- paste(unlist(strsplit(basename(file), "_"))[1:2], collapse="_")
  target <- unlist(strsplit(basename(file), "_"))[1]
  
  
  ###################################################################################################################################
  # Add new step in log file
  
  cat(paste(Sys.time(), "- Initiate Clustering with Vsearch for sample", sampleName, "\n"), file=logfile, sep="", append=T)
  
  
  ###################################################################################################################################
  # Dereplicate with vsearch
  
  vsearch_cmd <- paste0(vsearchpath, "vsearch")
  derep_args <- paste0("--derep_fulllength ", meep, file, " -output ", getwd(), "/Vsearch/", filename, "_derep.fasta -sizeout")
  CmdOut <- system2(command = vsearch_cmd, args = derep_args, stdout = TRUE, stderr = TRUE) # Execute dereplicate command
  cat(paste(Sys.time(), "- Dereplication with Vsearch for sample", sampleName, "\n"), file=logfile, sep="", append=T)
  cat(paste(Sys.time(), vsearch_cmd, derep_args, "\n"), file=logfile, sep="", append=T)
  
  
  ###################################################################################################################################
  # Collect information from vsearch derep command output for number of sequences, number of unique sequences and vsearch version.
  
  sequ <- as.numeric(sub(".*in (.*) seqs.*", "\\1", CmdOut[grep("nt in", CmdOut)]))
  derep <- as.numeric(sub("(.*)unique sequences.*", "\\1", CmdOut[grep("unique sequences", CmdOut)]))
  version <- sub("(.*)unique sequences.*", "\\1", temp[grep("vsearch", temp)][1])
  
  
  message(paste0("Reading ", filename,".fasta"))
  message(paste0(sequ, " input sequences -> ", derep, " dereplicated sequences"))
  cat(paste(Sys.time(), sampleName, ": ", sequ, "input sequences ->", derep, " ereplicated sequences\n"), file=logfile, sep="", append=T)

  
  ###################################################################################################################################
  # Sort dereplicated sequences
  
  sort_args <- paste0("--sortbysize ", getwd(), "/Vsearch/", filename, "_derep.fasta", " -output ", getwd(), "/Vsearch/", filename, "_sorted.fasta -minsize 2")
  CmdOut <- system2(command = vsearch_cmd, args = sort_args, stdout = TRUE, stderr = TRUE) # Execute sort command
  cat(paste(Sys.time(), "- Sort sequences with Vsearch minimum size 2 for sample", sampleName, "\n"), file=logfile, sep="", append=T)
  cat(paste(Sys.time(), vsearch_cmd, sort_args, "\n"), file=logfile, sep="", append=T)
  
  
  ###################################################################################################################################
  # Cluster sequences at 97% identity
  
  cluster_args <- paste0("--cluster_fast ", getwd(), "/Vsearch/", filename, "_sorted.fasta", " -id 0.97", " -centroids ", getwd(), "/Vsearch/", filename, "_cluster97.fasta -uc ", getwd(), "/Vsearch/", filename, "_cluster97.uc -sizein -sizeout")
  CmdOut <- system2(command = vsearch_cmd, args = cluster_args, stdout = TRUE, stderr = TRUE) # Execute cluster at 97% command
  Cluster97Seqs <- readDNAStringSet(paste0(getwd(), "/Vsearch/",filename, "_cluster97.fasta"))
  cat(paste(Sys.time(), "- Sort sequences with Vsearch minimum size 2 for sample", sampleName, "\n"), file=logfile, sep="", append=T)
  cat(paste(Sys.time(), vsearch_cmd, sort_args, "\n"), file=logfile, sep="", append=T)
  
  
  ###################################################################################################################################
  # BLAST sequences to filter
  
  Biostrings::writeXStringSet(PFAMseqs, paste0("PFAM_",PFAM,"_seqsBlastDB.fasta"),format = "fasta")
  cat(paste(Sys.time(), "Write PFAM nucletoide sequences", "\n"), file=logfile, sep="", append=T)
  makeblastdb_args <- paste0("-in ", getwd(), "/PFAMseqsBlastDB.fasta -dbtype nucl")
  system2(command = "makeblastdb", args = makeblastdb_args)
  cat(paste(Sys.time(), "Create PFAM Blast database\n"), file=logfile, sep="", append=T)
  bl <- blast(db = paste0(getwd(), paste0("PFAM_",PFAM,"_seqsBlastDB.fasta")))
  cat(paste(Sys.time(), "Blast clusters on PFAM blast database\n"), file=logfile, sep="", append=T)
  
  print(bl, info=TRUE)
  cat(paste(Sys.time(), "Filter blast results to remove unrelated sequences\n"), file=logfile, sep="", append=T)
  FilteredNames <- predict(bl, Cluster97Seqs, BLAST_args = "-evalue 1e-10")[1]
  Filtered97seqs <- Cluster97Seqs[(names(Cluster97Seqs)) %in% FilteredNames$QueryID]

  
  ###################################################################################################################################
  # Cluster filtered sequences at 95% identity
  
  cluster_args <- paste0("--cluster_fast ", getwd(), "/Vsearch/", filename, "_sorted.fasta -id 0.95 -centroids ", getwd(), "/Vsearch/", filename, "_cluster95.fasta -uc ", getwd(), "/Vsearch/", filename, "_cluster95.uc -sizein -sizeout")
  CmdOut <- system2(command = vsearch_cmd, args = cluster_args, stdout = TRUE, stderr = TRUE) # Execute cluster at 95% command
  Cluster95Seqs <- readDNAStringSet(paste0(getwd(), "/Vsearch/",filename, "_cluster97.fasta"))
  UCfile <- import_usearch_uc(paste0(getwd(), "/Vsearch/", filename, "_cluster95.uc"), readDelimiter = "_M03834")
  
  
  ###################################################################################################################################
  # phyloseq-tools rarefy
  
  # Change the Path!!!
  source("/Users/christophelemetre/Documents/Work/Primer_Design/R/phyloseq-tools.R")
  rarefaction_curve_data <- (calculate_rarefaction_curves(UCfile, c('Observed'), seq(1, max(sample_sums(UCfile)), length.out=100), parallel = FALSE))[c(1,4,2)]
  rarefaction_curve_data <- cbind(seq(1,100), rarefaction_curve_data)
  rarefaction_curve_data <- cbind(rarefaction_curve_data, rep(target,100))
  
  colnames(rarefaction_curve_data) <- c("idx","Depth","Diversity", "Sample", "Target")
  
  #summary(rarefaction_curve_data)
  
  return(rarefaction_curve_data)
  
  
  ###################################################################################################################################
  # Clean up
  # delete all the intermediate files.
  
  
  
  
}

df <- ''
PFAMlist <- c("PF13714", "PF13714", "PF01041")
fileList <- c("/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_01_combined.fq",
               "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/PEP_02_combined.fq",
               "/Users/christophelemetre/Documents/Work/eSNaPD3/Sample_fq/AHBA_02_combined.fq")
for (i in length(fileList)) {
  RAR <- primer_validation(file = fileList[i], PFAM = PFAMlist[i], vsearchpath = "integrated")
  df <- do.call("rbind", as.data.frame(RAR))
  #df2 <- rbind(df2, RAR)
}

fsdfsd
d1 <- lapply(Sys.glob("rarefactiondata/*.txt"), loadrarefaction)
d1 <- lapply(Sys.glob("Vsearch/*rarefied95.txt"), loadrarefaction)
df <- do.call(rbind,d1)

maxperprimer <- df %>% 
  arrange(desc(diversity)) %>%
  group_by(SampleName) %>%
  slice(1:1) %>%
  ungroup()


grouping <- rep(c("PEP_01","PEP_02","PEP_04","PEP_06","PEPcontrol"), each = 100)
RarefPlot <- ggplot(df, aes(x=depth,y=diversity, group=primer))+#, colour=grouping)) + 
  geom_line() + 
  #scale_color_manual(values=c("red","orange","yellow","blue","green")) +
  #geom_text_repel(data = maxperprimer,aes(label=primer)) +
  facet_wrap(~target, scales="free")

RarefPlot

























  
  rarefaction_depths = c(5000,10000)
  
  rarefy_obj <- function(physeq, i,minsize, depth=5000) {
    # prune samples
    ad95rare <- rarefy_even_depth(physeq, sample.size = depth ,trimOTUs = TRUE, rngseed = i)
    saveRDS(ad95rare, file=paste0("Vsearch/rarefaction/min",minsize,"/ad95_",depth,"_",i,".RDS"))
  }
  
  cores = 20
  
  for (depth in rarefaction_depths) {
    print(paste0("Rarefying to Depth of ", depth))
    #for (minsize in c(1,2,10,100,1000)) {
    for (minsize in c(2,5,10,100,1000)) {
      print(paste0("Using Minimum OTU size of ", minsize))
      physeqminsize <- prune_taxa(taxa_sums(PhyloSeqOTA) >= minsize, PhyloSeqOTA)
      rarefn <- function(i) {rarefy_obj(physeq=physeqminsize, i=i, minsize = minsize, depth=depth)}
      mclapply(1:100, rarefn, mc.cores = cores)
    }
  }
  ###############################  ###############################  ###############################  ###############################
  
  PhyloSeqObj <- phyloseq(PhyloSeqOTA, "PEP_TT")
  RarefPrimer <- rarefy_even_depth(PhyloSeqObj, sample.size = 1000 , rngseed = 500)
  plot_richness(PhyloSeqObj)
  ###############################  ###############################  ###############################
  
  # vegan rarefy
  data(BCI,package = "vegan")
  H <- vegan::diversity(BCI)
  simp <- diversity(BCI, "simpson")
  invsimp <- diversity(BCI, "inv")
  r.2 <- vegan::rarefy(BCI, 2)
  alpha <- fisher.alpha(BCI)
  pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
  ## Species richness (S) and Pielou's evenness (J):
  S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
  J <- H/log(S)
  
  
  
  
  
  
  
  # rarefaction plots
  
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(purrr)
  
  loadrarefaction <- function(filename) {
    basename <- strsplit(filename, "/")[[1]][[2]]
    primer   <- strsplit(basename, "\\_")[[1]][[2]]
    target   <- strsplit(basename, "\\_")[[1]][[1]]
    #target   <- strsplit(basename, "\\.")[[1]][[1]]
    #primer   <- "primer"
    
    read.table(filename, col.names = c("idx","depth","diversity"), stringsAsFactors = FALSE) %>%
      mutate(basename = basename,
             primer=primer,
             target=target)
  }
  
  
  rm(df)
  rm(d1)
  setwd("/users/clemetre/Project_Ulysses/AmpliconsBlastAnalysis/AmpliconsData/20160724_UlyssesPrimers_processing/Filtering_Reads/PARAMETER_TESTING/")
  #setwd("PEP/97cluster_prior_filter_then_95cluster/ClusterSize_greater_2/eVal_e-10/")
  ###Filtering_Reads/PEP_filtering_with_new_sequences_from_package/")
  
  d1 <- lapply(Sys.glob("rarefactiondata/*.txt"), loadrarefaction)
  d1 <- lapply(Sys.glob("Vsearch/*rarefied95.txt"), loadrarefaction)
  df <- do.call(rbind,d1)
  
  maxperprimer <- df %>% 
    arrange(desc(diversity)) %>%
    group_by(primer) %>%
    slice(1:1) %>%
    ungroup()
  
  
  grouping <- rep(c("PEP_01","PEP_02","PEP_04","PEP_06","PEPcontrol"), each = 100)
  RarefPlot <- ggplot(df, aes(x=depth,y=diversity, group=primer))+#, colour=grouping)) + 
    geom_line() + 
    #scale_color_manual(values=c("red","orange","yellow","blue","green")) +
    #geom_text_repel(data = maxperprimer,aes(label=primer)) +
    facet_wrap(~target, scales="free")
  
  RarefPlot
  
  
  
  
  
  
  
  
  
  
  
  



