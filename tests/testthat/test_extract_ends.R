################################################################################
# Use testthat to test extract_ends
################################################################################
library("DegeneratePrimerTools"); packageVersion("DegeneratePrimerTools")
library("Biostrings"); packageVersion("Biostrings")
library("testthat"); packageVersion("testthat")
context('Checking sequence retrieval from online databases')

test_that("Extract-Ends Preserves Names", {
  dnafile <- system.file("sequences","AHBA.fna",package="DegeneratePrimerTools")
  dnatest <- Biostrings::readDNAStringSet(dnafile)
  trimmed <- extract_ends(dnatest, trimf=5, trimr=3, sep="N")
  expect(all.equal(names(dnatest), names(trimmed)))
})

