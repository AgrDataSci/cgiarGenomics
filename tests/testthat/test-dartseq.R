# Audit findings R2 and R3. read_DArTSeq_SNP() dropped every locus when all
# SNPs were mapped, and could take positions from the wrong rows.

dart_file = testthat::test_path("..", "DartSeq_fmt", "DartSeq", "DArTSeq_CB_snp.csv")
chr_col = "Chrom_Common_bean_Phaseolus_acutifolius_v10"
pos_col = "ChromPosSnp_Common_bean_Phaseolus_acutifolius_v10"

test_that("a report where every SNP is mapped keeps its loci", {
  skip_if_not_installed("dartR.base")
  skip_if_not(file.exists(dart_file))

  # Build a report with mapped SNPs only, from the example file
  lines = readLines(dart_file)
  header = which(startsWith(lines, "AlleleID"))[1]
  fields = strsplit(lines[(header + 1):length(lines)], ",")
  cols = strsplit(lines[header], ",")[[1]]
  chr = vapply(fields, function(x) x[match(chr_col, cols)], character(1))
  pos = suppressWarnings(as.numeric(vapply(fields, function(x) x[match(pos_col, cols)], character(1))))
  mapped = which(!is.na(chr) & chr != "" & !is.na(pos) & pos > 0)
  mapped_file = tempfile(fileext = ".csv")
  on.exit(unlink(mapped_file), add = TRUE)
  writeLines(c(lines[seq_len(header)], lines[header + mapped]), mapped_file)

  gl = suppressWarnings(suppressMessages(
    read_DArTSeq_SNP(mapped_file, snp_id = "AlleleID",
                     chr_name = chr_col, pos_name = pos_col)
  ))

  expect_gt(adegenet::nLoc(gl), 0)
  expect_length(adegenet::position(gl), adegenet::nLoc(gl))
  expect_true(all(adegenet::position(gl) > 0))
})

test_that("unmapped SNPs are removed and positions match the loci kept", {
  skip_if_not_installed("dartR.base")
  skip_if_not(file.exists(dart_file))

  gl = suppressWarnings(suppressMessages(
    read_DArTSeq_SNP(dart_file, snp_id = "AlleleID",
                     chr_name = chr_col, pos_name = pos_col)
  ))

  expect_gt(adegenet::nLoc(gl), 0)
  expect_length(adegenet::position(gl), adegenet::nLoc(gl))
  expect_true(all(adegenet::position(gl) > 0))
  expect_false(anyNA(adegenet::chromosome(gl)))
})
