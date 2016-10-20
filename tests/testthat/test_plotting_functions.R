################################################################################
# Use testthat to test that onlinde sequence database retrieval
################################################################################
library("DegeneratePrimerTools"); packageVersion("DegeneratePrimerTools")
library("testthat"); packageVersion("testthat")
context('Checking Plotting Functions')

data("AHBA_degeprimer")

test_that("Plot_GC Works", {
  expect_is(plot_GC(AHBA_degeprimer), "ggplot")
})

test_that("Plot_MSA Works", {
  msa2 <- add_primers_to_MSA(AHBA_degeprimer)
  expect_is(plot_msa(AHBA_degeprimer@msa), "ggplot")
  expect_is(plot_msa(msa2), "ggplot")
})

test_that("Plot_CoverageMatrix Works", {
  expect_is(plot_coveragematrix(AHBA_degeprimer@msa), "ggplot")
  expect_is(plot_msa(msa2), "ggplot")
})

test_that("plot_degeprimer Works", {
  expect_is(plot_degeprimer(AHBA_degeprimer@primerdata), "ggplot")
})
