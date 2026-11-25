# Audit finding R1. Duplicated marker IDs were not renamed.

test_that("duplicated marker IDs are renamed and unique IDs are kept", {
  meta = data.frame(id = c("a", "b", "b"), chrom = c("1", "1", "1"),
                    pos = c(10, 20, 30), ref = "A", alt = "G")

  out = suppressMessages(cgiarGenomics:::process_metadata(meta))

  expect_equal(out$id, c("a", "b", "1_30"))
  expect_equal(anyDuplicated(out$id), 0)
  expect_false(any(out$filter))
})

test_that("multi-allelic markers and indels are flagged", {
  meta = data.frame(id = c("a", "b", "c", "d"), chrom = "1",
                    pos = c(10, 20, 30, 40),
                    ref = c("A", "A", "AT", "A"), alt = c("G", "G,T", "A", NA))

  out = suppressWarnings(suppressMessages(cgiarGenomics:::process_metadata(meta)))

  expect_equal(out$filter, c(FALSE, TRUE, TRUE, TRUE))
})
