# Audit finding B1. Removing from one margin only emptied the other margin.

mods_table = function(row, col) {
  data.frame(reason = "filter", row = row, col = col, value = NA)
}

test_that("removing only loci keeps every individual", {
  gl = recalc_metrics(make_gl(matrix(c(0, 1, 2, 0, 1, 2, 2, 1, 0, 0, 0, 1), nrow = 4)))

  out = apply_bioflow_modifications(gl, mods_table(row = NA, col = 2))

  expect_equal(adegenet::nInd(out), 4)
  expect_equal(adegenet::locNames(out), c("snp1", "snp3"))
})

test_that("removing only individuals keeps every locus", {
  gl = recalc_metrics(make_gl(matrix(c(0, 1, 2, 0, 1, 2, 2, 1, 0, 0, 0, 1), nrow = 4)))

  out = apply_bioflow_modifications(gl, mods_table(row = c(1, 4), col = NA))

  expect_equal(adegenet::indNames(out), c("ind2", "ind3"))
  expect_equal(adegenet::nLoc(out), 3)
})

test_that("removing from both margins still works", {
  gl = recalc_metrics(make_gl(matrix(c(0, 1, 2, 0, 1, 2, 2, 1, 0, 0, 0, 1), nrow = 4)))

  out = apply_bioflow_modifications(gl, mods_table(row = c(1, NA), col = c(NA, 3)))

  expect_equal(adegenet::indNames(out), c("ind2", "ind3", "ind4"))
  expect_equal(adegenet::locNames(out), c("snp1", "snp2"))
})

test_that("a table with imputation rows only returns the object unchanged", {
  gl = recalc_metrics(make_gl(matrix(c(0, 1, 2, 0, 1, 2, 2, 1, 0, 0, 0, 1), nrow = 4)))
  mods = data.frame(reason = "imputation", row = 1, col = 1, value = 2)

  expect_warning(out <- apply_bioflow_modifications(gl, mods), "not applied")

  expect_equal(dim(as.matrix(out)), c(4L, 3L))
})

test_that("imputation rows are reported as not applied", {
  gl = recalc_metrics(make_gl(matrix(c(0, 1, 2, 0, 1, 2, 2, 1, 0, 0, 0, 1), nrow = 4)))
  mods = data.frame(reason = c("filter", "imputation"), row = c(NA, 1),
                    col = c(2, 1), value = c(NA, 2))

  expect_warning(out <- apply_bioflow_modifications(gl, mods), "not applied")
  expect_equal(adegenet::locNames(out), c("snp1", "snp3"))
})
