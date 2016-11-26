#' primer_validation
#'
#' Library of functions to run eSNaPD.
#'
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

  
  
  
  
  
  
  
  



