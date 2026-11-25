# Audit finding I1. Random Forest imputation wrote factor codes as dosages.

test_that("random forest imputation writes dosages, not factor codes", {
  skip_if_not_installed("randomForest")

  # Six identical loci on one chromosome, with dosages 0 and 2 only. The
  # flanking loci predict the masked calls without error, and a dosage of 1
  # can only come from a factor code.
  truth = matrix(rep(rep(c(0, 2), 20), 6), nrow = 40)
  masked = truth
  masked[1:4, 3] = NA
  gl = make_gl(masked)

  out = suppressMessages(
    impute_gl(gl, ploidity = 2, method = "random_forest",
              nflank = 2, ntree = 50, seed = 1, nthreads = 1)
  )
  imputed = as.matrix(out$gl)

  expect_false(anyNA(imputed))
  expect_true(all(imputed[1:4, 3] %in% c(0, 2)))
  expect_equal(unname(imputed[1:4, 3]), c(0, 2, 0, 2))
  expect_equal(unname(imputed[-(1:4), ]), truth[-(1:4), ])
  expect_setequal(names(out$log), c("0", "2"))
})
