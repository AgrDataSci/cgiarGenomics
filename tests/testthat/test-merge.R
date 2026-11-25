# Audit finding P1. merge_duplicate_inds() matched the wrong column, failed
# with single samples only and dropped locus annotations.

panel = function() {
  mt = rbind(s1 = c(0, 0, 2, 1),
             s2 = c(0, 2, 2, 1),
             s3 = c(0, 1, 2, NA),
             s4 = c(1, 1, 0, 2),
             s5 = c(2, 2, 2, 2))
  colnames(mt) = paste0("snp", 1:4)
  make_gl(mt, chromosome = c("1", "1", "2", "2"))
}

test_that("single samples only are returned under their designation", {
  gl = panel()
  dict = data.frame(sample_id = paste0("s", 1:5), designation_id = paste0("s", 1:5))

  out = merge_duplicate_inds(gl, dict)

  expect_equal(as.matrix(out), as.matrix(gl))
})

test_that("replicates merge by majority call when sample and designation IDs differ", {
  gl = panel()
  dict = data.frame(sample_id = paste0("s", 1:5),
                    designation_id = c("G2", "G2", "G2", "G9", "G1"))

  out = merge_duplicate_inds(gl, dict)
  m = as.matrix(out)

  # singles first in the order of the object, then merged designations
  expect_equal(rownames(m), c("G9", "G1", "G2"))
  expect_equal(unname(m["G9", ]), c(1, 1, 0, 2))
  expect_equal(unname(m["G1", ]), c(2, 2, 2, 2))
  # snp1 0,0,0. snp2 0,2,1 is a tie. snp3 2,2,2. snp4 1,1,NA
  expect_equal(unname(m["G2", ]), c(0, NA, 2, 1))
})

test_that("the merged object keeps locus annotations and has metrics", {
  gl = panel()
  dict = data.frame(sample_id = paste0("s", 1:5),
                    designation_id = c("G2", "G2", "G2", "G9", "G1"))

  out = merge_duplicate_inds(gl, dict)

  expect_equal(adegenet::locNames(out), adegenet::locNames(gl))
  expect_equal(adegenet::position(out), adegenet::position(gl))
  expect_equal(as.character(adegenet::chromosome(out)), c("1", "1", "2", "2"))
  expect_equal(adegenet::alleles(out), adegenet::alleles(gl))
  expect_equal(nrow(out@other$ind.metrics), 3)
  expect_equal(nrow(out@other$loc.metrics), 4)
})

test_that("all individuals in one designation merge to one row", {
  gl = panel()
  dict = data.frame(sample_id = paste0("s", 1:5), designation_id = "G1")

  out = merge_duplicate_inds(gl, dict)

  expect_equal(adegenet::nInd(out), 1)
  expect_equal(adegenet::indNames(out), "G1")
})

test_that("an individual missing from the dictionary is an error", {
  gl = panel()
  dict = data.frame(sample_id = paste0("s", 1:4), designation_id = "G1")

  expect_error(merge_duplicate_inds(gl, dict))
})
