# Audit finding M2. Heterozygosity counted missing calls in the denominator.

test_that("heterozygosity is the share of heterozygotes among called genotypes", {
  mt = cbind(c(1, 1, NA, NA),   # two heterozygotes, two missing
             c(0, 1, 2, 2),
             c(0, 0, 0, NA))
  gl = make_gl(mt)

  expect_equal(unname(get_loc_heterozygosity(gl)), c(1, 0.25, 0))
  expect_equal(unname(get_ind_heterozygosity(gl)), c(1 / 3, 2 / 3, 0, 0))

  both = get_heterozygosity_metrics(gl, ploidy = 2)
  expect_equal(unname(both$het_loc), c(1, 0.25, 0))
  expect_equal(unname(both$het_ind), c(1 / 3, 2 / 3, 0, 0))
})

test_that("partial heterozygotes count in a tetraploid", {
  mt = cbind(c(0, 1, 3, 4), c(2, 2, NA, 4))
  gl = make_gl(mt, ploidy = 4)

  expect_equal(unname(get_loc_heterozygosity(gl, ploidy = 4)), c(0.5, 2 / 3))
})
