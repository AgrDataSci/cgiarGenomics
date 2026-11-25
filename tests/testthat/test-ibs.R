# ibs_matrix_purrr() computes identity by state with matrix products. These
# tests compare it with the pairwise definition in ibs_dosage().

pairwise_ibs = function(gl, ploidy) {
  G = as.matrix(gl)
  n = nrow(G)
  score = matrix(NA_real_, n, n)
  overlap = matrix(0L, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      res = cgiarGenomics:::ibs_dosage(G[i, ], G[j, ], ploidy)
      score[i, j] = res$score
      overlap[i, j] = as.integer(res$overlap)
    }
  }
  diag(score) = 1
  dimnames(score) = list(rownames(G), rownames(G))
  dimnames(overlap) = dimnames(score)
  list(ibs = score, overlap = overlap)
}

test_that("diploid IBS matches the pairwise definition", {
  set.seed(3)
  mt = matrix(sample(0:2, 12 * 30, replace = TRUE), nrow = 12)
  mt[sample(length(mt), 60)] = NA
  gl = make_gl(mt)

  expect_equal(ibs_matrix_purrr(gl, ploidy = 2), pairwise_ibs(gl, 2))
})

test_that("tetraploid IBS matches the pairwise definition", {
  set.seed(4)
  mt = matrix(sample(0:4, 10 * 25, replace = TRUE), nrow = 10)
  mt[sample(length(mt), 40)] = NA
  gl = make_gl(mt, ploidy = 4)

  expect_equal(ibs_matrix_purrr(gl, ploidy = 4), pairwise_ibs(gl, 4))
})

test_that("pairs with no shared locus get NA and overlap 0", {
  mt = rbind(c(0, NA, 2),
             c(NA, 1, NA),
             c(0, 1, 2))
  gl = make_gl(mt)

  out = ibs_matrix_purrr(gl, ploidy = 2)

  expect_true(is.na(out$ibs[1, 2]))
  expect_equal(out$overlap[1, 2], 0L)
  expect_equal(unname(diag(out$ibs)), c(1, 1, 1))
  expect_equal(out$ibs[1, 3], 1)
  expect_equal(out$ibs[2, 3], 1)
})

test_that("non-integer dosages use the pairwise fallback", {
  mt = rbind(c(0, 0.5, 2), c(1, 1.5, 2), c(2, 2, 0))
  gl = make_gl(mt)

  expect_equal(ibs_matrix_purrr(gl, ploidy = 2), pairwise_ibs(gl, 2))
})
