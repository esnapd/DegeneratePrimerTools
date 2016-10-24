################################################################################
# Use testthat to test that basic pocesisng funcitonality is working
################################################################################
library("DegeneratePrimerTools"); packageVersion("DegeneratePrimerTools")
library("Biostrings"); packageVersion("Biostrings")
library("testthat"); packageVersion("testthat")
context('Checking processing functions')

data("AHBA_degeprimer")

test_that("Primers are being added to the MSA", {
  expect_is(add_primers_to_MSA(AHBA_degeprimer, position = 120), "DNAMultipleAlignment")
  
})
