################################################################################
# Use testthat to test that onlinde sequence database retrieval
################################################################################
library("DegeneratePrimerTools"); packageVersion("DegeneratePrimerTools")
library("Biostrings"); packageVersion("Biostrings")
library("testthat"); packageVersion("testthat")
context('Checking sequence retrieval from online databases')

test_that("PFAM Ids are being correctly downloaded from the PFAM Website", {

  df <- retrieve_PFAM_ids("PF16997", alignmenttype = "uniprot")
  expect_is(df, "data.frame")
  expect_equal(names(df), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(df), c(131,4))

  # check non-default alignmenttypes
  # seed
  dfseed <- retrieve_PFAM_ids("PF16997", alignmenttype = "seed")
  expect_is(dfseed, "data.frame")
  expect_equal(names(dfseed), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfseed), c(14,4))
  # full
  dffull <- retrieve_PFAM_ids("PF16997", alignmenttype = "full")
  expect_is(dffull, "data.frame")
  expect_equal(names(dffull), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dffull), c(18,4))
  # rp15
  dfrp15 <- retrieve_PFAM_ids("PF16997", alignmenttype = "rp15")
  expect_is(dfrp15, "data.frame")
  expect_equal(names(dfrp15), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfrp15), c(2,4))
  # rp35
  dfrp35 <- retrieve_PFAM_ids("PF16997", alignmenttype = "rp35")
  expect_is(dfrp35, "data.frame")
  expect_equal(names(dfrp35), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfrp35), c(5,4))
  # rp55
  dfrp55 <- retrieve_PFAM_ids("PF16997", alignmenttype = "rp55")
  expect_is(dfrp55, "data.frame")
  expect_equal(names(dfrp55), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfrp55), c(14,4))
  # rp75
  dfrp75 <- retrieve_PFAM_ids("PF16997", alignmenttype = "rp75")
  expect_is(dfrp75, "data.frame")
  expect_equal(names(dfrp75), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfrp75), c(17,4))
  # ncbi
  dfncbi <- retrieve_PFAM_ids("PF16997", alignmenttype = "ncbi")
  expect_is(dfncbi, "data.frame")
  expect_equal(names(dfncbi), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfncbi), c(61,4))
  # meta
  dfmeta <- retrieve_PFAM_ids("PF07752", alignmenttype = "meta")
  expect_is(dfmeta, "data.frame")
  expect_equal(names(dfmeta), c("PFAM_ID", "Accession", "start", "end"))
  expect_equal(dim(dfmeta), c(2,4))

  # check incorrect alignment type
  expect_error(retrieve_PFAM_ids("PF16997", alignmenttype = "wrongthing"))
  # check for errors in non-existent REST endpoints
  expect_error(retrieve_PFAM_ids("PF16997", alignmenttype = "meta"), "Invalid HTTP request. Check that your PFAM ID is correct and that the alignment type is available. For example the ncbi and meta are not always avaiable.")
  expect_error(retrieve_PFAM_ids("PF00501", alignmenttype = "ncbi"))
  # check for errors in non-existent REST endpoints
  expect_error(retrieve_PFAM_ids("P00501", alignmenttype = "seed"),  "pfamids are prefixed with 'PF'")
  # check that dots and have been reomved from the return Accesison values
  expect_false( any(grep("\\.", df$Accession)) )
})


test_that("Uniprot IDs can be converted to EMBL ids", {
  pfamdf    <- retrieve_PFAM_ids("PF16997", alignmenttype = "uniprot")
  uniprotdf <- retrieve_UNIPROT_to_EMBL(pfamdf$Accession)

  expect_is(uniprotdf, "data.frame")
  expect_equal(names(uniprotdf), c("UNIPROT_ID", "EMBL_ID"))
  expect_equal(dim(pfamdf),    c(131,4))
  expect_equal(dim(uniprotdf), c(135,2))

  #test a single ID
  expect_identical(retrieve_UNIPROT_to_EMBL("Q4SMD5"), data.frame(UNIPROT_ID=c("Q4SMD5"), EMBL_ID=c("CAF98197.1"),stringsAsFactors = FALSE))
})

test_that("EMBL/ENA nucleotide sequnces are fetched correctly", {
  fnas <- retrieve_EMBL_sequences(c("A00145","A00146"))
  expect_is(fnas, "DNAStringSet")
  expect_equal(length(fnas), 2)
})

test_that("EMBL/ENA Retrieval Code Correctly handles suppressed sequences", {
  #good, good, bad
  seqs <- retrieve_EMBL_sequences(c("KGE72166.1","KJE27958.1","AIE20100.1"))
  expect_equal(length(seqs), 2)
})

test_that("All retrieval functions work as an integrated pipeline", {

  seqdf <-  retrieve_PFAM_nucleotide_sequences("PF16997")
  expect_is(seqdf, "data.frame")
  expect_equal(dim(seqdf), c(135,8))

  #test in frame translation
  for (i in 1:10) {
    dna <- seqdf$domainsequence[[i]]
    expect_is(translate(DNAString(dna)), "AAString")
  }
})

