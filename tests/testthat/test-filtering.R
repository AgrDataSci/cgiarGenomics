# Audit findings F1 and F3.

three_loci = function() {
  recalc_metrics(make_gl(cbind(c(0, 0, 0, 2), c(0, 1, 1, 2), c(2, 2, 1, 2))))
}

test_that("a filter step that nothing passes stops with an error", {
  gl = three_loci()
  steps = list(list(param = "maf", operator = ">", threshold = 0.9))

  expect_error(apply_sequence_filtering(gl, steps),
               class = "cgiarGenomics_empty_filter")
})

test_that("filter_gl lists every name as removed when nothing passes", {
  gl = three_loci()

  out = filter_gl(gl, "maf", 0.9, ">")

  expect_length(out$index, 0)
  expect_equal(out$filter_out, c("snp1", "snp2", "snp3"))
})

test_that("a filter step that everything passes keeps the object", {
  gl = three_loci()
  steps = list(list(param = "loc_miss", operator = "<=", threshold = 0.5))

  out = apply_sequence_filtering(gl, steps)

  expect_equal(adegenet::nLoc(out$gl), 3)
  expect_equal(adegenet::nInd(out$gl), 4)
  expect_length(out$filt_log$loc_miss_1$filter_out, 0)
})

test_that("filters are applied in sequence and logged", {
  gl = three_loci()
  steps = list(list(param = "maf", operator = ">", threshold = 0.2))

  out = apply_sequence_filtering(gl, steps)

  expect_equal(adegenet::locNames(out$gl), c("snp1", "snp2"))
  expect_equal(out$filt_log$maf_1$filter_out, "snp3")
})

test_that("random_select_loci works when two chromosomes tie for the most loci", {
  set.seed(1)
  mt = matrix(sample(0:2, 10 * 8, replace = TRUE), nrow = 10)
  mt[1, ] = 0
  mt[2, ] = 2
  gl = recalc_metrics(make_gl(mt, chromosome = rep(c("1", "2"), each = 4)))

  out = random_select_loci(gl, ind_miss = 1, loc_miss = 1, maf = 0,
                           size = 5, seed = 1)

  expect_equal(adegenet::nLoc(out), 5)
})

test_that("random_select_loci keeps the loci in the order of the object", {
  set.seed(2)
  mt = matrix(sample(0:2, 10 * 12, replace = TRUE), nrow = 10)
  mt[1, ] = 0
  mt[2, ] = 2
  gl = recalc_metrics(make_gl(mt, chromosome = rep(c("1", "2", "3"), c(6, 4, 2))))

  out = random_select_loci(gl, ind_miss = 1, loc_miss = 1, maf = 0,
                           size = 6, seed = 1)
  kept = match(adegenet::locNames(out), adegenet::locNames(gl))

  expect_equal(adegenet::nLoc(out), 6)
  expect_false(is.unsorted(kept))
  # 6 x 6/12 = 3 from chromosome 1, 2 from 2, 1 from 3
  expect_equal(as.vector(table(as.character(adegenet::chromosome(out)))), c(3, 2, 1))
})
