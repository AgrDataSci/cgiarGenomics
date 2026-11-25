# Audit findings R4, R6 and R9, on small hapmap files written by the test.

calls = rbind(c("AA", "AG", "GG", NA),
              c("CC", "CT", "TT", "CT"))

test_that("a diploid hapmap is read as alternative allele dosages", {
  path = write_hapmap(calls, alleles = c("A/G", "C/T"))
  on.exit(unlink(path))

  gl = suppressMessages(read_hapmap(path, ploidity = 2))

  expect_equal(unname(as.matrix(gl)), cbind(c(0, 1, 2, NA), c(0, 1, 2, 1)))
  expect_equal(adegenet::indNames(gl), paste0("ind", 1:4))
  expect_equal(adegenet::locNames(gl), c("snp1", "snp2"))
  expect_equal(adegenet::position(gl), c(100L, 200L))
})

test_that("calls written with a separator are read", {
  sep_calls = rbind(c("A/A", "A/G", "G/G", NA),
                    c("C/C", "C/T", "T/T", "T/C"))
  path = write_hapmap(sep_calls, alleles = c("A/G", "C/T"))
  on.exit(unlink(path))

  gl = suppressMessages(read_hapmap(path, ploidity = 2, sep = "/"))

  expect_equal(unname(as.matrix(gl)), cbind(c(0, 1, 2, NA), c(0, 1, 2, 1)))
})

test_that("single-letter IUPAC calls are decoded", {
  iupac_calls = rbind(c("A", "R", "G", NA),
                      c("C", "Y", "T", "Y"))
  path = write_hapmap(iupac_calls, alleles = c("A/G", "C/T"))
  on.exit(unlink(path))

  gl = suppressMessages(read_hapmap(path, ploidity = 2))

  expect_equal(unname(as.matrix(gl)), cbind(c(0, 1, 2, NA), c(0, 1, 2, 1)))
})

test_that("a tetraploid hapmap is read", {
  tetra = rbind(c("AAAA", "AAAG", "AAGG", "GGGG"))
  path = write_hapmap(tetra, alleles = "A/G")
  on.exit(unlink(path))

  gl = suppressMessages(read_hapmap(path, ploidity = 4))

  expect_equal(as.vector(as.matrix(gl)), c(0, 1, 2, 4))
})

test_that("the diploid example file gives the hand-counted MAF", {
  path = testthat::test_path("..", "hmp_fmt", "diploid_hmp.txt")
  skip_if_not(file.exists(path))

  gl = suppressWarnings(suppressMessages(read_hapmap(path, ploidity = 2)))
  m = as.matrix(gl)
  f = colSums(m, na.rm = TRUE) / (2 * colSums(!is.na(m)))

  expect_equal(gl@other$loc.metrics$maf, unname(pmin(f, 1 - f)))
})
