library(DegeneratePrimerTools)
context('Checking the Codon Tables')

test_that("the 3 basic degeneracy tables work", {
  expect_equal(makedegenerateseqs("ACT"), "ACN")
  expect_equal(makedegenerateseqs("AGG", method="S"), "MGN")
  expect_equal(makedegenerateseqs("AGG", method="Z"), "CGN")
  expect_equal(makedegenerateseqs("AGG", method="SZ"), "MGN")

  expect_equal(makedegenerateseqs("TCG", method="SZ"), "NNN")
  expect_equal(makedegenerateseqs("TCR", method="SZ"), "NNN")
  expect_equal(makedegenerateseqs("TCY", method="SZ"), "NNN")
  expect_equal(makedegenerateseqs("TCW", method="SZ"), "NNN")

})

test_that("make degenerateseqs can hancdle multiple types of text", {
  expect_equal(makedegenerateseqs("ACT"), "ACN")
  expect_equal(makedegenerateseqs(c("ACT","ACT")), c("ACN","ACN"))
  expect_equal(makedegenerateseqs("ACTACT"), "ACNACN")
  expect_equal(makedegenerateseqs("ACTACTT"), "ACNACNT")
  expect_equal(makedegenerateseqs(c("ACTACTT","ACTACTT")), c("ACNACNT","ACNACNT"))

})
