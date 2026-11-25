# Frequency imputation, masking and accuracy. Audit findings I2, I3 and I6.

test_that("a complete panel is returned unchanged with an empty log", {
  mt = matrix(rep(c(0, 1, 2, 1), 6), nrow = 8)
  gl = make_gl(mt)

  out = suppressMessages(impute_gl(gl, ploidity = 2, method = "frequency"))

  expect_equal(unname(as.matrix(out$gl)), mt)
  expect_length(out$log, 0)
})

test_that("frequency imputation fills gaps and keeps observed calls", {
  set.seed(1)
  mt = matrix(sample(0:2, 60, replace = TRUE), nrow = 10)
  masked = mt
  masked[c(3, 17, 42)] = NA
  gl = make_gl(masked)

  out = suppressMessages(impute_gl(gl, ploidity = 2, method = "frequency"))
  imputed = unname(as.matrix(out$gl))

  expect_false(anyNA(imputed))
  expect_equal(imputed[-c(3, 17, 42)], mt[-c(3, 17, 42)])
  expect_true(all(imputed %in% 0:2))
  expect_setequal(unlist(out$log, use.names = FALSE), c(3, 17, 42))
  expect_equal(adegenet::position(out$gl), adegenet::position(gl))
})

test_that("loci with no data are removed with a warning", {
  mt = cbind(c(0, 1, 2, NA), rep(NA_real_, 4), c(2, 2, 1, 0))
  gl = make_gl(mt)

  expect_warning(
    out <- suppressMessages(impute_gl(gl, ploidity = 2, method = "frequency")),
    "all missing"
  )
  expect_equal(adegenet::locNames(out$gl), c("snp1", "snp3"))
  expect_false(anyNA(as.matrix(out$gl)))
})

test_that("a ploidy that disagrees with the object gives a warning", {
  gl = make_gl(cbind(c(0, 1, 2, NA), c(2, 2, 1, 0)))

  expect_warning(
    suppressMessages(impute_gl(gl, ploidity = 4, method = "frequency")),
    "ploidy"
  )
})

test_that("set_fraction_na keeps the object and masks called cells only", {
  mt = matrix(rep(c(0, 1, 2, 1), 10), nrow = 8)
  mt[1:4] = NA
  gl = make_gl(mt, chromosome = rep(c("1", "2"), c(2, 3)))

  set.seed(1)
  masked = set_fraction_na(gl, fraction = 0.25)
  m = as.matrix(masked)

  expect_equal(sum(is.na(m)), 4 + 10)
  expect_length(masked@other$masked_cells, 10)
  expect_false(any(masked@other$masked_cells %in% 1:4))
  expect_true(all(is.na(m[masked@other$masked_cells])))
  expect_equal(adegenet::indNames(masked), adegenet::indNames(gl))
  expect_equal(adegenet::locNames(masked), adegenet::locNames(gl))
  expect_equal(adegenet::position(masked), adegenet::position(gl))
  expect_equal(as.character(adegenet::chromosome(masked)),
               as.character(adegenet::chromosome(gl)))
  expect_equal(adegenet::alleles(masked), adegenet::alleles(gl))
})

test_that("get_accuracy scores the masked cells only", {
  truth = matrix(rep(c(0, 2), 30), nrow = 10)
  ref = make_gl(truth)
  masked_mt = truth
  masked_mt[1, 1] = NA
  masked_mt[2:3, 2] = NA
  masked = make_gl(masked_mt)
  imputed_mt = truth
  imputed_mt[1, 1] = 2        # wrong
  imputed_mt[2, 2] = 0        # wrong, truth is 2
  imputed = make_gl(imputed_mt)

  acc = get_accuracy(ref, imputed, masked_gl = masked)

  expect_equal(as.vector(acc), c(0, 0.5, NA, NA, NA, NA))
  expect_equal(unname(attr(acc, "n_scored")), c(1, 2, 0, 0, 0, 0))
})

test_that("get_accuracy warns when no mask is given", {
  truth = matrix(rep(c(0, 2), 30), nrow = 10)
  imputed_mt = truth
  imputed_mt[1, 1] = 2

  expect_warning(acc <- get_accuracy(make_gl(truth), make_gl(imputed_mt)), "inflated")
  expect_equal(as.vector(acc)[1], 0.9)
})

test_that("random forest imputation works with one flanking marker", {
  skip_if_not_installed("randomForest")
  truth = matrix(rep(rep(c(0, 2), 20), 2), nrow = 40)
  masked = truth
  masked[1:2, 1] = NA
  gl = make_gl(masked)

  out = suppressMessages(
    impute_gl(gl, ploidity = 2, method = "random_forest",
              ntree = 50, seed = 1, nthreads = 1)
  )

  expect_equal(unname(as.matrix(out$gl)), truth)
})

test_that("random forest imputation is reproducible with a seed", {
  skip_if_not_installed("randomForest")
  set.seed(42)
  mt = matrix(sample(0:2, 40 * 8, replace = TRUE), nrow = 40)
  mt[sample(length(mt), 30)] = NA
  gl = make_gl(mt)
  run = function() {
    suppressMessages(as.matrix(
      impute_gl(gl, ploidity = 2, method = "random_forest",
                nflank = 3, ntree = 20, seed = 7, nthreads = 1)$gl
    ))
  }

  first = run()
  second = run()

  expect_equal(first, second)
  expect_true(all(first %in% 0:2))
})
