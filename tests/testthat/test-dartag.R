# read_DArTag_count_dosage() on the sweetpotato example. It checks the base R
# parsing of the counts file, which pairs the Ref and Alt tag of each marker.

dose_file = testthat::test_path("..", "DartSeq_fmt", "DArT_TAG_sweet_potato",
                                "DP23-8340_Allele_Dose_Report.csv")
count_file = testthat::test_path("..", "DartSeq_fmt", "DArT_TAG_sweet_potato",
                                 "DP23-8340_Allele_match_counts_collapsed.csv")

test_that("tag pairs give the SNP alleles", {
  pair = cgiarGenomics:::compare_tags("TAAGGGGG", "TAAGGAGG")
  expect_equal(unname(pair), c("G", "A"))
  expect_equal(unname(cgiarGenomics:::compare_tags("TAAG", "TCCG")), c("N", "N"))
})

test_that("the sweetpotato DArTag files are read", {
  skip_if_not(file.exists(dose_file) && file.exists(count_file))

  gl = suppressWarnings(suppressMessages(
    read_DArTag_count_dosage(dose_file, count_file, ploidity = 6)
  ))

  expect_gt(adegenet::nLoc(gl), 0)
  expect_gt(adegenet::nInd(gl), 0)
  expect_equal(anyDuplicated(adegenet::locNames(gl)), 0)
  expect_true(all(grepl("^[ACGT]/[ACGT]$", adegenet::alleles(gl))))
  expect_equal(adegenet::alleles(gl)[adegenet::locNames(gl) == "PotVar0000025"], "G/A")
  m = as.matrix(gl)
  expect_true(all(m[!is.na(m)] %in% 0:6))
})
