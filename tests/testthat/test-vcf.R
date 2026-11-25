# Audit finding R5 and the export of missing polyploid genotypes.

write_small_vcf = function(lines) {
  path = tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3\tind4",
               lines), path)
  path
}

test_that("a VCF that mixes phased and unphased calls is read", {
  skip_if_not_installed("vcfR")
  path = write_small_vcf(c(
    "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0|1\t1|1\t./.",
    "1\t200\tsnp2\tC\tT\t.\tPASS\t.\tGT\t0|0\t1/0\t1/1\t.|."
  ))
  on.exit(unlink(path))

  gl = suppressWarnings(suppressMessages(read_vcf(path, ploidity = 2)))

  expect_equal(unname(as.matrix(gl)), cbind(c(0, 1, 2, NA), c(0, 1, 2, NA)))
  expect_equal(adegenet::locNames(gl), c("snp1", "snp2"))
})

test_that("the old phased separator argument still works", {
  skip_if_not_installed("vcfR")
  path = write_small_vcf("1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t1|0")
  on.exit(unlink(path))

  gl = suppressWarnings(suppressMessages(read_vcf(path, ploidity = 2, sep = "\\|")))

  expect_equal(as.vector(as.matrix(gl)), c(0, 1, 2, 1))
})

test_that("write_vcf writes missing genotypes at the ploidy of the data", {
  mt = cbind(c(0, 1, NA, 4), c(2, NA, 3, 0))
  gl = make_gl(mt, ploidy = 4)
  path = tempfile(fileext = ".vcf")
  on.exit(unlink(path))

  write_vcf(gl, path)
  body = grep("^#", readLines(path), invert = TRUE, value = TRUE)
  fields = strsplit(body, "\t")

  expect_equal(fields[[1]][10:13], c("0/0/0/0", "0/0/0/1", "./././.", "1/1/1/1"))
  expect_equal(fields[[2]][10:13], c("0/0/1/1", "./././.", "0/1/1/1", "0/0/0/0"))
})

test_that("a tetraploid panel survives a write and read round trip", {
  skip_if_not_installed("vcfR")
  mt = cbind(c(0, 1, NA, 4), c(2, NA, 3, 0), c(4, 4, 1, 0))
  gl = make_gl(mt, ploidy = 4)
  path = tempfile(fileext = ".vcf.gz")
  on.exit(unlink(path))

  write_vcf(gl, path)
  back = suppressWarnings(suppressMessages(read_vcf(path, ploidity = 4)))

  expect_equal(as.matrix(back), as.matrix(gl))
  expect_equal(adegenet::position(back), adegenet::position(gl))
  expect_equal(adegenet::alleles(back), adegenet::alleles(gl))
})
