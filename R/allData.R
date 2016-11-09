################################################################################
#' small Degeprime Object derived from 11 AHBA sequences
#'
#' 
#' Degeprimer glues together data needed for finding and assessing degeenrate primers.
#' This is a small degeprimer object derived from 11 AHBA sequences using the process 
#' found in the package vignette.
#'
#' 
#' @name data-ahba
#' @aliases AHBA_degeprimer
#' @docType data
#' @keywords data
#' @examples
#' data(AHBA_degeprimer)
#' plot_gc(AHBA_degeprime)
#' msa <- add_primers_to_MSA(AHBA_degeprime)
#' plot_msa(msa)
################################################################################
NA

#' Codon List for "Z" method of Degen.pl
#'
#' A named list of codons corresponding to the "Z"
#' variables of the Degen.pl algorithm.
#'
#'
#' @format A named list where names are codons sequences
#' and the values are the degenerate codon for that sequence
#' using the "Z" method of Degen.pl
#'
"degenZ"

#' Codon List for "S" method of Degen.pl
#'
#' A named list of codons corresponding to the "S"
#' variables of the Degen.pl algorithm.
#'
#'
#' @format A named list where names are codons sequences
#' and the values are the degenerate codon for that sequence
#' using the "S" method of Degen.pl
#'
"degenS"

#' Codon List for "SZ" method of Degen.pl
#'
#' A named list of codons corresponding to the "SZ"
#' variables of the Degen.pl algorithm.
#'
#'
#' @format A named list where names are codons sequences
#' and the values are the degenerate codon for that sequence
#' using the "SZ" method of Degen.pl
#'
"degenSZ"

#' p7 adaptor sequences for paired end illumina barcoding
#'
#' These sequences are the Illumina barcodes used for the indexing read
#' on paired end runs. Intended for the MiSeq instrument.
#'
#' @format A two-column data frame containing the full p7 sequence and
#' the barcode needed for the sample file.
#'
"illumina_p7"
