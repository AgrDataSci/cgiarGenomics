# Beagle wrapper. Settings are checked before any process starts, and markers
# are sorted for Beagle and matched back by name. Audit findings I4 and I5.

beagle_panel = function() {
  set.seed(11)
  founders = matrix(sample(0:1, 6 * 40, replace = TRUE), nrow = 6)
  mt = t(vapply(seq_len(20), function(i) {
    founders[sample(6, 1), ] + founders[sample(6, 1), ]
  }, numeric(40)))
  make_gl(mt, chromosome = rep(c("1", "2"), each = 20),
          position = rep(seq_len(20) * 10000L, 2))
}

test_that("text that is not a memory setting is refused before anything runs", {
  gl = beagle_panel()
  marker = tempfile()

  expect_error(impute_beagle(gl, memory = paste0("Xmx1g; touch ", marker)), "memory")
  expect_error(impute_beagle(gl, memory = "Xmx1g -version"), "memory")
  expect_error(impute_beagle(gl, memory = c("Xmx1g", "Xms1g")), "memory")
  expect_false(file.exists(marker))
})

test_that("settings that are not whole numbers are refused", {
  gl = beagle_panel()

  expect_error(impute_beagle(gl, burnin = "3; ls"), "burnin")
  expect_error(impute_beagle(gl, iterations = 1.5), "iterations")
  expect_error(impute_beagle(gl, nthreads = "2 & ls"), "nthreads")
})

test_that("Beagle fills gaps in an unsorted panel and keeps its order", {
  skip_if_not_installed("vcfR")
  skip_if_not(check_beagle_requirements()$ready, "Java or the Beagle JAR is missing")

  gl = beagle_panel()
  truth = as.matrix(gl)
  shuffled = sample(adegenet::nLoc(gl))
  gl = gl[, shuffled]
  mt = as.matrix(gl)
  holes = c(5, 250, 611)
  mt[holes] = NA
  gl_na = make_gl(mt, chromosome = as.character(adegenet::chromosome(gl)),
                  position = adegenet::position(gl))

  out = suppressWarnings(suppressMessages(
    impute_gl(gl_na, ploidity = 2, method = "beagle",
              memory = "Xmx1g", nthreads = 1, seed = 1)
  ))
  imputed = as.matrix(out$gl)

  expect_false(anyNA(imputed))
  expect_equal(adegenet::locNames(out$gl), adegenet::locNames(gl_na))
  expect_equal(imputed[-holes], mt[-holes])
  expect_true(all(imputed %in% 0:2))
  expect_setequal(unlist(out$log, use.names = FALSE), holes)
})

test_that("Beagle returns an empty log for a complete panel", {
  skip_if_not_installed("vcfR")
  skip_if_not(check_beagle_requirements()$ready, "Java or the Beagle JAR is missing")

  gl = beagle_panel()

  expect_length(suppressWarnings(suppressMessages(
    impute_beagle(gl, memory = "Xmx1g", nthreads = 1)
  )), 0)
})

test_that("polyploid objects are refused", {
  gl = make_gl(cbind(c(0, 1, 4, 2), c(2, 2, 1, 0)), ploidy = 4)

  expect_error(impute_beagle(gl), "diploid")
})
