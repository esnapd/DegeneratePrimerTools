#' obtain the start and stop locations for a PFAM ID
#'
#' Use a PFAM ID to retrieve the PFAM alignment. PFAM maintains alignments for the seed sequences as well as
#' representative sequences (rp15, rp35, rp55, rp75), as well as sequecnes from NCBI and UNIPPROT.  The alignments are in Stockholm
#' format and contain protein information
#'
#' @param pfamid. Required. A string coresponding to a PFAM accession id.
#' @param alignmnettype. Optional. Will determine which of PFAM's prebuilt alignments to download. Can choose from "seed"
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

    rows <- scan(text=content(r,"text"), what = "", sep = "\n", strip.white = TRUE, quiet = TRUE, blank.lines.skip = FALSE)

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
  # no rows where and N and C domains overlap to omuch or are negative
  # only keep adjacent domains within a reasonable amino acid distance of one another
  together <- nterm[cterm][Cstart - Nend > - allowableoverlap][ Cstart - Nend < gapsize]
  
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
retrieve_UNIPROT_to_EMBL <- function(uniprotids, chunksize=200) {

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

      df <- read.table(text=content(r, "text"), header = TRUE, stringsAsFactors = FALSE)
      names(df) <- c("UNIPROT_ID", "EMBL_ID")

      return(df)
    }
  )
}
#' Retrieve nucleotide sequences from EMBL ENA
#'
#' Use the REST API to retrieve sequences from EMBL
#'
#' @importFrom httr GET
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings union
#' @return a DNAStringSet
#' @export
#' @seealso \url{https://www.ebi.ac.uk/ena/browse/data-retrieval-rest}
#' @examples
#' retrieve_EMBL_sequences(c("A00145","A00146"))
retrieve_EMBL_sequences <- function(emblids, chunksize=100) {

  baseurl <- "https://www.ebi.ac.uk/ena/data/view/"
  opts    <- "&display=fasta&download=text"

  if (is.character(emblids)) emblids <- c(emblids)

  # split into groups no greater in size than the chunksize, retrieve and write those IDs to file
  # then read them in
  groups <- split(emblids, ceiling(seq_along(emblids)/chunksize))
  temps <- lapply(groups, function(x){return(tempfile())})

  seqs <- lapply(seq_along(groups), function(i) {

    paste0("Attempting to retrieve: ", paste0(groups[[i]], collapse= " "))
    r <- GET(paste0(baseurl, paste0(groups[[i]], collapse=","), opts))

    # remove lines where sequences return error messages
    lines <- scan(text=content(r, "text"), what = "", sep = "\n", strip.white = TRUE, quiet = TRUE, blank.lines.skip = FALSE)
    lines <- lines[!grepl("has been suppressed at the submitter's request on", lines)]
    write(x=lines, file=temps[[i]])
  })

  dnas <- lapply(temps, readDNAStringSet)
  return(Reduce(append, dnas))
}
#' A one-stop-shop for obtaining nucelotides from PFAM sequences.
#'
#' @param pfamids. Required. A string corresponding to a PFAM accession id.
#' @param alignmnettype. Optional. Will determine which of PFAM's prebuilt alignments to download. Can choose from "seed"
#'             "full", "rp15", "rp35", "rp55", "rp75", "uniprot", "ncbi", "meta".
#' @param gapsize. Optional. Default \code{100}. Maximum distance in aminoacids that two domains can be if you want them to be combined.
#' @param allowableoverlap. Optional. Default \code{15}. Maximum ovelap between two domains to still be included.
#'              
#' @importFrom purrr map
#' @importFrom purrr reduce
#' @export
#' @examples
#' # WAP1 Domain
#' retrieve_PFAM_nucleotide_sequences("PF16997")
#' 
#' \dontrun{
#' # Aconitase Domain
#' retrieve_PFAM_nucleotide_sequences("PF00330", alignmenttype = "seed")
#'
#' # Aconitase Domain plus Aconitase-C  Domain
#' retrieve_PFAM_nucleotide_sequences(c("PF00330", "PF00694"), alignmenttype = "uniprot")
#' 
#'}
retrieve_PFAM_nucleotide_sequences <- function(pfamids, alignmenttype = "uniprot", gapsize = 100, allowableoverlap = 15) {
  
  if (!is.vector(pfamids)) { pfamids <- c(pfamids)}
  pfamdfs    <- map(pfamids, ~retrieve_PFAM_ids(., alignmenttype=alignmenttype))
  
  if (length(unique(pfamdfs$PFAM_ID)) == 1) {
    pfamdf <- pfamdfs[[1]]
  } else {
    pfamdf <- reduce(pfamdfs, ~combine_domains(.x, .y, gapsize=gapsize, allowableoverlap = allowableoverlap))
  }

  if (nrow(pfamdf) == 0) stop("There are no sequences that have these adjacent domains. You can try changing the alignment type, the
                             gapsize or the allowable overlapbetween sequences.")
  
  
  uniprotids <- unique(pfamdf$Accession)
  embldf     <- retrieve_UNIPROT_to_EMBL(uniprotids)
  seqs       <- retrieve_EMBL_sequences(unique(embldf$EMBL_ID))
  #collate the data
  masterdf<- merge(pfamdf, embldf, by.x="Accession", by.y="UNIPROT_ID")
  #calculate DNA start/stop positions
  masterdf$dnastart <- 3 * (masterdf$start - 1) + 1 #(n-1 * 3)  + 1
  masterdf$dnaend   <- 3 * masterdf$end             #(n * 3)
  # remove the version number form the EMBL IDs
  masterdf$EMBL_ID_clean <- gsub("\\..*$", "", masterdf$EMBL_ID)
  #pull out the dna
  masterdf$domainsequence <- mapply(
    function(seqname, start, end, sequences=seqs) {
      #pull out a sequence form the DNAstinrg set using the name, stop, and end.
      targetsequence <- sequences[grepl(seqname, names(sequences))]

      # sequence not found can be legitimate
      if (length(targetsequence) == 0) {
        warning(paste0(seqname, " is not found in the header of any sequences and may have been removed or suppressed from the ENA"))
        return(NA_character_)
      }

      #multiple sequences is a problem
      if (length(targetsequence) > 1)  stop("There is an error in the ENA fasta headers or in the EMBL ID list")

      # I encoutnered this in some archael sequences and we shoudl flag it whent his happend.
      if (end > width(targetsequence)) {
        warning(paste0("The width of ", seqname, "is too narrow for the specified value. This can happen in certain sequences ldue to a framshift. See http://www.uniprot.org/uniprot/D3E4V4 and its
                       associated nucleotide fasta for an example. This sequence is not being included and should be verified."))
        return(NA_character_)
      }

      dnadomain <- subseq(targetsequence, start=start, end=end)
      return(as.character(dnadomain))
      },
    masterdf$EMBL_ID_clean,
    masterdf$dnastart,
    masterdf$dnaend)

  return(masterdf[c("PFAM_ID", "Accession", "start", "end", "EMBL_ID", "dnastart", "dnaend", "domainsequence")])
}

