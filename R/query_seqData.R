#' Functions to interrogate sequencing data for a particular reference set.
#'
#' @param Type Required. Type of the analysis, either array (with plate well information) or crude (with geolocalisation information) (string)
#' @param runFolder Required. path location for the folder with the sequencing data (string)
#' @param metaData Required. path location for the file containing metadata for the sequencing data (string)
#' @param Ref Required. Reference sequences for the target investigated. Either file location (string) or DNAstringSet or Degeprimer RDS object.
#' 
#' 
#' @import ggplot2
#' @importFrom Biostrings DNAStringSet
#' @export
query_SeqData <- function(Type, runFolder, metaData, Ref){
  
  runFolder <- "/Volumes/data/AmpliconProcessing/Processing/20160730_DFD_AD2/DFDPlate1_AD2/processed/"
  print(runFolder)
  sampleSeqFiles <- list.files(runFolder, pattern = "*.fasta", full.names = TRUE)
  sampleList <- sapply(1:length(sampleSeqFiles), function(x) unlist(strsplit(basename(sampleSeqFiles[x]), "\\.r00"))[1])
  
  Ref <- DNAStringSet("/Volumes/data/Christophe/Reference_sequences_Primer_validation/AHBA_Ulysse.fna")
  
  MolecResults <- NULL
  sum =0
  for (i in 1:2){#96){#length(sampleSeqFiles)){
    RES <- NULL
    
    RES <- run_blast(sampleSeqFiles[i], Ref, parallel = TRUE, blast_args = "-evalue 1e-30 -max_target_seqs 1")
    if (!is.null(RES)){
      MolecResults <- rbind(MolecResults,
                            data.frame(Sample = unlist(strsplit(RES$queryID, "\\.r00"))[1],
                                       Reference = unlist(strsplit(RES$subjectID[1], "_"))[3],
                                       eValue = RES$evalue[1]))
      }
    
  }
  colnames(MolecResults)[1] <- "Sample"
  PlateLayout <- data.frame(Sample=sampleList[1:96], rown = rep (letters[1:8], 12), coln = rep (1:12, each = 8))
  plateData <- merge(PlateLayout, MolecResults, by="Sample", all.x = TRUE)
  geoData <- read.csv("/Volumes/data/DFDData/CitizenScience/_data/DFD_Samples.csv")
  geoData <- merge(MolecResults, geoData, by.x="Sample", by.y="DFD_UID")
  
  if(Type == "Geo"){
    return(DegeneratePrimerTools::draw_plate_results(96, plateData))}
  if(Type == "Plate"){
    return(DegeneratePrimerTools::draw_geomap_results(geoData))}
  
  
  
}