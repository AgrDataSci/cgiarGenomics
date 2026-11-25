# Audit finding M1. get_maf() divided the allele frequency by ploidy twice.

test_that("get_maf returns the true minor allele frequency in a diploid", {
  mt = cbind(c(0, 0, 0, 2),   # alternative frequency 2/8
             c(2, 2, 2, 2),   # fixed for the alternative allele
             c(1, 1, NA, 2),  # 4/6 over the non-missing calls
             c(1, 1, 1, 1))   # 4/8
  gl = make_gl(mt)

  expect_equal(unname(get_maf(gl, ploidy = 2)), c(0.25, 0, 1 / 3, 0.5))
})

test_that("get_maf does not depend on the ploidy argument", {
  mt = cbind(c(0, 4, 4, 4),   # 12/16, so MAF is 0.25
             c(1, 1, 1, 1))   # 4/16
  gl = make_gl(mt, ploidy = 4)

  expect_equal(unname(get_maf(gl, ploidy = 4)), c(0.25, 0.25))
  expect_equal(get_maf(gl, ploidy = 4), get_maf(gl))
})

test_that("recalc_metrics stores the corrected MAF and filters use it", {
  mt = cbind(c(0, 0, 0, 2),
             c(0, 1, 1, 2),
             c(2, 2, 2, 2))
  gl = recalc_metrics(make_gl(mt))

  expect_equal(gl@other$loc.metrics$maf, c(0.25, 0.5, 0))

  # MAF of 0.25 passes a threshold of 0.2. The old value, 0.125, did not.
  out = filter_gl(gl, "maf", 0.2, ">")
  expect_equal(unname(out$index), c(1L, 2L))
})
