#' obtain the start and stop locations for a PFAM ID
#'
#' Use a PFAM ID to retrieve the PFAM alignment. PFAM maintains alignments for the seed sequences as well as
#' representative sequences (rp15, rp35, rp55, rp75), as well as sequecnes from NCBI and UNIPPROT.  The alignments are in Stockholm
#' format and contain protein information
#'
#' @param pfamid Required. A string coresponding to a PFAM accession id.
#' @param alignmnettype Optional. Will determine which of PFAM's prebuilt alignments to download. Can choose from "seed"
#'             "full", "rp15", "rp35", "rp55", "rp75", "uniprot", "ncbi", "meta".
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom purrr map_df
#' @seealso \url{http://pfam.xfam.org}
#' @export
#' @examples
#' retrieve_PFAM_ids("PF16997")
retrieve_PFAM_ids <- function(pfamid, alignmenttype = "uniprot") {
  alignmenttypes <- c("seed", "full", "rp15", "rp35", "rp55", "rp75", "uniprot", "ncbi", "meta")
  if (!alignmenttype %in% alignmenttypes) stop(paste("Alignment type must be one of: ", paste0(alignmenttypes, collapse=" ")))
  
  pfamsite <- paste0("http://pfam.xfam.org/family/", pfamid ,"/alignment/", alignmenttype)
  
  try(r <- GET(pfamsite))
  
  if (!grepl("^PF", pfamid)) stop("pfamids are prefixed with 'PF'")
  if (r$status_code == 400)  stop("Invalid HTTP request. Check that your PFAM ID is correct and that the alignment type is available. For example the ncbi and meta are not always avaiable.")
  if (r$status_code == 500)  stop("Invalid HTTP request. Status Code == 500. Error with PFAM server. Try later?")
  
  
  if (is.null(r)) stop("There was a problem downloading the PFAM")
  
  rows <- scan(text=content(r,"text", encoding="UTF8"), what = "", sep = "\n", strip.white = TRUE, quiet = TRUE, blank.lines.skip = FALSE)
  
  if (length(rows) < 3 || !identical(grep("^# STOCKHOLM", rows[1L]), 1L)) stop("invalid Stockholm file")
  
  markupLines <- grep("(^\\s*|^#.*|^//\\s*)$", rows, perl=TRUE)
  seqrows <- rows[!1:length(rows) %in% markupLines]
  seqrows <- strsplit(seqrows, "\\s+")
  
  df <- suppressWarnings(
    map_df(seqrows, function(s) {
      id    <- strsplit(s[[1]], "/")[[1]][[1]]
      rng   <- strsplit(s[[1]], "/")[[1]][[2]]
      start <- strsplit(rng, "-")[[1]][[1]]
      end   <- strsplit(rng, "-")[[1]][[2]]
      data.frame(PFAM_ID = pfamid, Accession = id, start=as.numeric(start),end=as.numeric(end))
    }))
  
  df$Accession <- gsub("\\..*$", "", df$Accession) #remove trailing version number
  df
}
#' Combine the PFAM info for adjacent Protein Domains
#'
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table setkey
#' @export
combine_domains <- function(nterm, cterm, gapsize = 100, allowableoverlap = 15) {
  protnames <- c("PFAM_ID", "Accession", "start", "end")
  stopmesg  <-"This function is meant to work on the four column dataframe returned by
  retrieve_PFAM_ids()"
  
  if (!all.equal(names(nterm),protnames)) stop(stopmessg)
  if (!all.equal(names(cterm),protnames)) stop(stopmessg)
  
  names(nterm) <- c("NPFAM_ID", "Accession", "Nstart", "Nend")
  names(cterm) <- c("CPFAM_ID", "Accession", "Cstart", "Cend")
  nterm <- as.data.table(nterm)
  cterm <- as.data.table(cterm)
  setkey(nterm, "Accession")
  setkey(cterm, "Accession")
  
  # merge
  # no rows where and N and C domains overlap too much or are negative
  # only keep adjacent domains within a reasonable amino acid distance of one another
  together <- nterm[cterm, allow.cartesian = TRUE][Cstart - Nend > - allowableoverlap][ Cstart - Nend < gapsize]
  
  together[, PFAM_ID := mapply(function(x,y) {paste(x,y,sep="_")}, NPFAM_ID, CPFAM_ID)]
  
  #fix some columns and return
  together$start   <- together$Nstart
  together$end     <- together$Cend
  together[, c("PFAM_ID", "Accession", "start", "end"), with=FALSE]
}
#' Get EMBL IDs from Uniprot IDs
#'
#' Use the REST interface at UNIPROT to get ENA mappings for a protein.
#'
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom purrr map_df
#' @seealso \url{http://www.uniprot.org/help/programmatic_access}
#' @export
#' @examples
#' retrieve_UNIPROT_to_EMBL("Q4SMD5")
#' retrieve_UNIPROT_to_EMBL(c("Q4SMD5", "A3CRU7"))
retrieve_UNIPROT_to_EMBL <- function(uniprotids, chunksize=400) {
  
  baseurl <- "http://www.uniprot.org/mapping/"
  
  if (is.character(uniprotids)) uniprotids <- c(uniprotids)
  
  # split into groups no greater in size than the chunksize
  groups <- split(uniprotids, ceiling(seq_along(uniprotids)/chunksize))
  
  # batch request ID mapping using the REST API
  df <- map_df(
    groups,
    function(ids) {
      try(r <- POST(url=baseurl, query=list( from="ACC", to="EMBL", query= paste0(ids, collapse = " "), format="tab")))
      
      if (is.null(r)) stop("There was a problem obtaining the UNIPROT->ENA Mapping")
      
      df <- read.table(text=content(r, "text", encoding="UTF8"), header = TRUE, stringsAsFactors = FALSE)
      names(df) <- c("UNIPROT_ID", "EMBL_ID")
      
      return(df)
    }
  )
}
#' Get UNIPROT IDS from EMLB/Gennbank IDs
#'
#' Use the REST interface at UNIPROT to get UNIPROT mappings for a protein
#' given the EMBL/Genbank ID
#'
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom purrr map_df
#' @seealso \url{http://www.uniprot.org/help/programmatic_access}
#' @export
#' @examples
#' retrieve_EMBL_to_UNIPROT("AEK75490.1")
#' retrieve_EMBL_to_UNIPROT(c("AEK75490.1", "AEK75491.1"))
retrieve_EMBL_to_UNIPROT <- function(uniprotids, chunksize=400) {
  
  baseurl <- "http://www.uniprot.org/mapping/"
  
  if (is.character(uniprotids)) uniprotids <- c(uniprotids)
  
  # split into groups no greater in size than the chunksize
  groups <- split(uniprotids, ceiling(seq_along(uniprotids)/chunksize))
  
  # batch request ID mapping using the REST API
  df <- map_df(
    groups,
    function(ids) {
      try(r <- POST(url=baseurl, query=list( from="EMBL", to="ACC", query= paste0(ids, collapse = " "), format="tab")))
      
      if (is.null(r)) stop("There was a problem obtaining the UNIPROT->ENA Mapping")
      
      df <- read.table(text=content(r, "text"), header = TRUE, stringsAsFactors = FALSE)
      names(df) <- c("EMBL_ID", "UNIPROT_ID")
      
      return(df)
    }
  )
}
#' Retrieve nucleotide sequences from EMBL ENA
#'
#' Use the REST API to retrieve sequences from EMBL
#'
#' @param emblids. Required. A character or character vector of EMBL nucleotide accessions.
#' @param chunksize. Optional. Default \code{200}. The size of the batch request to send to EMBL.
#' @importFrom httr GET
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings union
#' @return a DNAStringSet
#' @export
#' @seealso \url{https://www.ebi.ac.uk/ena/browse/data-retrieval-rest}
#' @examples
#' retrieve_EMBL_sequences(c("A00145","A00146"))
retrieve_EMBL_sequences <- function(emblids, chunksize=400) {
  
  baseurl <- "https://www.ebi.ac.uk/ena/data/view/"
  opts    <- "&display=fasta&download=text"
  
  if (is.character(emblids)) emblids <- c(emblids)
  
  # split into groups no greater in size than the chunksize, retrieve and write those IDs to file
  # then read them in
  groups <- split(emblids, ceiling(seq_along(emblids)/chunksize))
  temps  <- lapply(groups, function(x){return(tempfile())})
  
  seqs <- lapply(seq_along(groups), function(i) {
    
    paste0("Attempting to retrieve: ", paste0(groups[[i]], collapse= " "))
    r <- GET(paste0(baseurl, paste0(groups[[i]], collapse=","), opts))
    
    # remove lines where sequences return error messages
    lines <- scan(text=content(r, "text", encoding="UTF8"), what = "", sep = "\n", strip.white = TRUE, quiet = TRUE, blank.lines.skip = FALSE)
    lines <- lines[!grepl("has been suppressed at the submitter's request on", lines)]
    
    write(x=lines, file=temps[[i]])
  })
  
  
  try(dnas <- lapply(temps, readDNAStringSet))
  try(dna  <- Reduce(append, dnas))
  if (exists("dna")) return(dna)
  
  #if dna doesn't exist, print the offending values
  
  warning("DNA retrieved from EMBL contains non-fasta sequences in the return values.
          Attempting to locate and print the offending values")
  
  for (tmp in temps) {
    tryCatch(dna <- readDNAStringSet(tmp),
             #finally = print(paste(scan(file=tmp, what = "", sep = "\n", strip.white = TRUE, quiet = TRUE, blank.lines.skip = FALSE), tmp))
             finally = file.copy(tmp, "~/Desktop/badfile.txt")
    )
  }
}
#' A one-stop-shop for obtaining nucelotides from PFAM sequences.
#'
#' @param pfamids Required. A string corresponding to a PFAM accession id, a vector of PFAM strings, or a semi-colon delimited string of PFAM ids.
#' @param alignmnettype Optional. Will determine which of PFAM's prebuilt alignments to download. Can choose from "seed"
#'             "full", "rp15", "rp35", "rp55", "rp75", "uniprot", "ncbi", "meta".
#' @param gapsize Optional. Default \code{100}. Maximum distance in aminoacids that two domains can be if you want them to be combined.
#' @param allowableoverlap Optional. Default \code{15}. Maximum ovelap between two domains to still be included.
#' @param retrivaltype Optional. Default \code{"rest"}. PFAM data is retrieved wither by the PFAM REST api or by the prepackaged DB.
#' @param pfamdb Optional. Default \code{NULL}. Location of the PFAM Database. To use this option, th retrievaltyp emust be set to 'db'
#'
#' @importFrom purrr map
#' @importFrom purrr reduce
#' @importFrom Biostrings width
#' @importFrom Biostrings subseq
#' @importFrom parallel mcmapply
#' @export
#' @examples
#' # WAP1 Domain
#' retrieve_PFAM_nucleotide_sequences("PF16997")
#'
#' \dontrun{
#' # SemiColon Format
#' retrieve_PFAM_nucleotide_sequences("PF00330;PF00694"), alignmenttype = "uniprot")
#'
#' # Aconitase Domain
#' retrieve_PFAM_nucleotide_sequences("PF00330", alignmenttype = "seed")
#'
#' # Aconitase Domain plus Aconitase-C  Domain
#' retrieve_PFAM_nucleotide_sequences(c("PF00330", "PF00694"), alignmenttype = "uniprot")
#'
#' # Use the db functionality
#' retrieve_PFAM_nucleotide_sequences("PF16997", pfamdb = "PFAM_Domains.sqlite", retrievaltype = "db")
#' retrieve_PFAM_nucleotide_sequences(c("PF00330"), pfamdb = "PFAM_Domains.sqlite", retrievaltype = "db")
#'}
retrieve_PFAM_nucleotide_sequences <- function(pfamids, alignmenttype = "uniprot", gapsize = 100, allowableoverlap = 15, retrievaltype="db", pfamdb=NULL) {
  
  ################################################################
  ################################################################
  ### Handle pfam IDs
  # check for the presence of semicolons. If there, split and return
  if (length(pfamids)==1 & grepl(";", pfamids)) {
    pfamids <- strsplit(pfamids,";")[[1]]
  }
  
  # handle list/character pfams
  if (!is.vector(pfamids)) { pfamids <- c(pfamids)}
  
  ################################################################
  ################################################################
  ### Obtain PFAM-> UNIPROT mappings
  
  message("Obtaining UNIPROT IDs from PFAM IDs....")
  if (retrievaltype == "db") {
    
    if (is.null(pfamdb)) stop("pfamdb must specify a path to a valid PFAM database if the retrievaltype is set to 'db' ")
    
    pfamdfs <- map(pfamids, ~select_PFAMs(pfamdb=pfamdb, pfamids=.))
    pfamdfs <- map(pfamdfs, function(df) {
      names(df) <- names(df) <- c("PFAM_ID", "Accession", "start", "end")
      df})
    
  } else if (retrievaltype == "rest"){
    pfamdfs    <- map(pfamids, ~retrieve_PFAM_ids(., alignmenttype=alignmenttype))
  } else {
    stop("retrievaltype must be either 'db' or 'rest' ")
  }
  
  
  if (length(unique(pfamdfs$PFAM_ID)) == 1) {
    pfamdf <- pfamdfs[[1]]
  } else {
    pfamdf <- reduce(pfamdfs, ~combine_domains(.x, .y, gapsize=gapsize, allowableoverlap = allowableoverlap))
  }
  
  
  if (nrow(pfamdf) == 0) stop(
    "There are PFAM records returned from your search. There are several possiblities that can cause this.
    First, if using a Database check that your PFAM is there. For example, the PFAM server contians many
    domains, but the defualt DB for eSNaPD only has PFAM-A sequences.  If you entered more than
    one PFAM domain you shoudl also verify that the combination of domaisn you are interested in exists,
    and that you have appropriate 'gapsize' and 'overlapbetween' sequence settings.")
  
  ################################################################
  ################################################################
  ### UNIPROT-> EMBL -> ENA mappings
  message("Obtaining UNIPROT  to EMBL mappings ....")
  uniprotids <- unique(pfamdf$Accession)
  embldf     <- retrieve_UNIPROT_to_EMBL(uniprotids)
  message("Obtaining EMBL sequences ....")
  seqs       <- retrieve_EMBL_sequences(unique(embldf$EMBL_ID))
  
  
  ################################################################
  ################################################################
  ### Process the output Data.
  
  message("Processing EMBL sequences to obtain nucleotide sequences of input PFAM domains....")
  #collate the data
  masterdf<- merge(pfamdf, embldf, by.x="Accession", by.y="UNIPROT_ID")
  #calculate DNA start/stop positions
  masterdf$dnastart <- 3 * (masterdf$start - 1) + 1 #(n-1 * 3)  + 1
  masterdf$dnaend   <- 3 * masterdf$end             #(n * 3)
  # remove the version number form the EMBL IDs
  masterdf$EMBL_ID <- gsub("\\..*$", "", masterdf$EMBL_ID)
  
  # extract all of the domain sequences; save the master data.table in case of failure.
  domainsequences <- NULL
  try(domainsequences <- mcmapply(
    FUN = function(x) {extract_sequences(id,start,end, sequences = seqs)},
    masterdf$EMBL_ID,
    masterdf$dnastart,
    masterdf$dnaend)
  )
  
  # save the downloaded
  if (is.null(domainsequences)) {
    pfamid <- paste(pfamids, collapse="_")
    name   <- paste0(pfamid, "_fulloutput.RDS")
    message(paste("there was an error in extracting your DNA sequences. Saving to ", name, "."))
    saveRDS(masterdt, file = name)
  } else {
    #pull out the dna
    masterdf$domainsequence <- domainsequences
    # this is needed becuase data table has a different semantic for selecting columns thant data.frame does
    # this will be changed in data.table 1.9.8 so it can be changed them.
    # now I am forcing the conversion to data.table and then selecting columns, being sure to use the 'with='
    # clause.
    masterdt <- as.data.table(masterdf)
    
    return(masterdt[,c("PFAM_ID", "Accession", "start", "end", "EMBL_ID", "dnastart", "dnaend", "domainsequence"), with=FALSE])
  }
}

#' Extraction function
#'
#' Helper function to retrieve a sequecne form a list of sequecnes, and extract a particular portion of it.
#'
#' @importFrom Biostrings width
#' @importFrom Biostrings subseq
extract_sequences <-  function(seqname, start, end, sequences) {
  #pull out a sequence form the DNAstring set using the name, stop, and end.
  targetsequence <- sequences[grepl(seqname, names(sequences))]
  
  # sequence not found can be legitimate
  if (length(targetsequence) == 0) {
    warning(paste0(seqname, " is not found in the header of any sequences and may have been removed or suppressed from the ENA"))
    return(NA_character_)
  }
  
  # multiple sequences is a problem
  if (length(targetsequence) > 1)  {
    target_names <- paste(names(targetsequence), collapse = ",")
    warning(paste("There is apossible error in the ENA fasta headers or in the EMBL ID list. Returning the first sequence. The
                  accession values or this target are", target_names))
    
  }
  
  # I encountered this in some archael sequences and we should flag it when this happend.
  if (end > width(targetsequence)) {
    warning(paste0("The width of ", seqname, "is too narrow for the specified value. This can happen in certain sequences due to a frameshift. See http://www.uniprot.org/uniprot/D3E4V4 and its
                   associated nucleotide fasta for an example. This sequence is not being included and should be verified."))
    return(NA_character_)
  }
  
  dnadomain <- subseq(targetsequence, start=start, end=end)
  return(as.character(dnadomain))
}
