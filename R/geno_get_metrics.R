# Purity Metrics ----------------------------------------------------------
#' Identity-by-state matrix from allelic dosages
#'
#' Computes pairwise identity by state (IBS) between all individuals of a
#' genlight object. For each locus typed in both individuals, IBS is
#' `1 - |dosage_x - dosage_y| / ploidy`, and the pair score is the mean over
#' those loci. Used to check purity and identity of samples.
#'
#' @param gl A genlight object.
#' @param ploidy Integer. Ploidy level used to scale dosage differences.
#'
#' @return A list with two square matrices named by individual.
#'   `ibs` holds the mean IBS score (diagonal set to 1), and `overlap` holds
#'   the number of loci typed in both individuals.
#'
#' @details Runs in pure R and evaluates every pair twice, so time grows with
#'   the square of the number of individuals. Use a subset of loci for large
#'   panels, for example with [random_select_loci()].
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' ibs_matrix_purrr(gl, ploidy = 2)
ibs_matrix_purrr <- function(gl, ploidy) {
  G <- t(as.matrix(gl))
  cols <- as.data.frame(G, check.names = FALSE)
  ids  <- colnames(cols)
  score_list <- purrr::map(cols, function(x)
    purrr::map_dbl(cols, function(y) ibs_dosage(x, y, ploidy)$score)
  )
  score_mat   <- do.call(cbind, score_list)
  rownames(score_mat) <- ids
  colnames(score_mat) <- ids
  # numeric safety for diagonal (optional)
  diag(score_mat) <- 1
  
  overlaps_list <- purrr::map(cols, function(x)
    purrr::map_int(cols, function(y) ibs_dosage(x, y, ploidy)$overlap)
  )
  overlap_mat   <- do.call(cbind, overlaps_list)
  rownames(overlap_mat) <- ids
  colnames(overlap_mat) <- ids
  return(list(ibs = score_mat, overlap = overlap_mat))
}

#' IBS score between two dosage vectors (internal)
#'
#' @param x,y Numeric dosage vectors of equal length.
#' @param ploidy Ploidy level.
#' @return A list with `score` (mean IBS) and `overlap` (loci typed in both).
#' @noRd
ibs_dosage <- function(x, y, ploidy) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(ploidy) != 1 || !is.numeric(ploidy) || ploidy <= 0) {
    stop("ploidy must be a positive number")
  }
  # mask non-missing pairs
  ok <- !is.na(x) & !is.na(y)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(list(score = NA, overlap = 0L))
  
  xv <- x[ok]; yv <- y[ok]
  
  # locus-wise IBS: 1 - |Δdosage| / P
  ibs_locus <- 1 - abs(xv - yv) / ploidy
  # numeric safety clamp
  ibs_locus <- pmin(pmax(ibs_locus, 0), 1)

  list(score = mean(ibs_locus), overlap = n_ok)
}

# Missingness metrics -----------------------------------------------------

#' Get Locus Missingness
#'
#' Compute the missing rate for each locus (marker) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that locus, and 1 indicates that all samples have missing data for that locus.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of loci (markers),
#'   containing the missing rate for each locus.
#'
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_loc_missing(gl)
get_loc_missing <- function(gl) {
  # Get the number of occurrences of NAs (missing data) for each marker
  NA_counts <- adegenet::glNA(gl)
  
  # Divide the NA counts by the total number of samples to get the missing rate
  NA_counts <- NA_counts / adegenet::nInd(gl)/max(adegenet::ploidy(gl))
  
  return(NA_counts)
}

#' Get Individual Missingness
#'
#' Compute the missing rate for each individual (sample) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that individual, and 1 indicates that all loci have missing data for that individual.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of individuals (samples),
#'   containing the missing rate for each individual.
#'
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_ind_missing(gl)
get_ind_missing <- function(gl) {
  # Convert the genlight object to a matrix and identify missing genotype calls
  mt <- is.na(as.matrix(gl))
  
  # Calculate the proportion of missing data for each individual (row)
  ind_miss <- Matrix::rowSums(mt) / adegenet::nLoc(gl)
  
  return(ind_miss)
}


#' Get Overall Missingness
#'
#' Compute the overall missing rate for a genlight object.
#' The overall missing rate is the proportion of missing genotype calls
#' across all individuals and loci in the dataset.
#'
#' @param gl A genlight object.
#'
#' @return A single numeric value representing the overall missing rate.
#'
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_overall_missingness(gl)
get_overall_missingness <- function(gl) {
  # Convert the genlight object to a matrix
  mt <- as.matrix(gl)
  
  # Identify missing genotype calls
  mt <- is.na(mt)
  
  # Calculate the overall missing rate
  overall_miss <- sum(mt) / (nrow(mt) * ncol(mt))
  
  return(overall_miss)
}


# Heterozygosity metrics --------------------------------------------------

#' Observed heterozygosity by individual and by locus
#'
#' A call is heterozygous when its dosage is above 0 and below `ploidy`. The
#' function returns the share of heterozygous calls per individual and per
#' locus in one pass.
#'
#' @param gl A genlight object.
#' @param ploidy Integer. Ploidy level (default 2).
#'
#' @return A list with `het_ind` (one value per individual) and `het_loc`
#'   (one value per locus).
#' @seealso [get_ind_heterozygosity()], [get_loc_heterozygosity()]
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_heterozygosity_metrics(gl, ploidy = 2)
get_heterozygosity_metrics <- function(gl, ploidy = 2){
  # Boolean matrix of genotype calls where 0 > dosage < ploidy
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  # share of heterozygotes among the calls that are not missing
  called <- !is.na(mt)
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/Matrix::colSums(called)
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/Matrix::rowSums(called)
  return(list(het_ind = het_ind, het_loc = het_loc))
}

#' Observed heterozygosity by locus
#'
#' Share of individuals with a heterozygous call (dosage above 0 and below
#' `ploidy`) at each locus.
#'
#' @param gl A genlight object.
#' @param ploidy Integer. Ploidy level (default 2).
#'
#' @return A numeric vector with one value per locus.
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_loc_heterozygosity(gl)
get_loc_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  # share of heterozygotes among the calls that are not missing
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/Matrix::colSums(!is.na(mt))
  return(het_loc)
}

#' Observed heterozygosity by individual
#'
#' Share of loci with a heterozygous call (dosage above 0 and below
#' `ploidy`) in each individual.
#'
#' @param gl A genlight object.
#' @param ploidy Integer. Ploidy level (default 2).
#'
#' @return A numeric vector with one value per individual.
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_ind_heterozygosity(gl)
get_ind_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  # share of heterozygotes among the calls that are not missing
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/Matrix::rowSums(!is.na(mt))
  return(het_ind)
}


#' Minor allele frequency by locus
#'
#' Computes the frequency of the alternative allele with
#' [adegenet::glMean()] and folds it to the minor allele.
#'
#' @param gl A genlight object.
#' @param ploidy Integer. Ploidy level (default 2).
#'
#' @return A numeric vector with one value per locus, between 0 and 0.5.
#'
#' @details
#' [adegenet::glMean()] returns the frequency of the alternative allele for any
#' ploidy, so `ploidy` is not used in the calculation. The argument is kept so
#' that existing calls keep working. Versions up to 1.0.2 divided the frequency
#' by `ploidy` a second time and returned values that were too low.
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' get_maf(gl)
get_maf <- function(gl, ploidy = 2){
  # glMean() already returns the alternative allele frequency
  alf <- adegenet::glMean(gl)
  maf <- ifelse(alf > 0.5, 1 - alf, alf)
  return(maf)
}

#' Per-locus inbreeding coefficient (internal)
#'
#' `1 - Ho / He` with `He = 2pq`, from the metrics stored by
#' [recalc_metrics()]. Valid for diploids only.
#'
#' @param gl A genlight object with `other$loc.metrics`.
#' @return A numeric vector with one value per locus.
#' @noRd
get_inbreeding <- function(gl){
  p <- 1 - gl@other$loc.metrics$maf
  he <- 2 * gl@other$loc.metrics$maf * p
  Fis <- 1 - (gl@other$loc.metrics$loc_het/he)
  return(Fis)
}

#' Recalculate locus and individual metrics
#'
#' Rebuilds the two tables that filtering relies on and stores them in the
#' `other` slot. Every reader and every filtering or imputation step calls it,
#' so metrics always describe the current object.
#'
#' `other$loc.metrics` holds `maf`, `loc_miss`, `loc_het` and `loc_Fis`, one
#' row per locus. `other$ind.metrics` holds `ind_miss` and `ind_het`, one row
#' per individual.
#'
#' @param gl A genlight object.
#'
#' @return The same genlight object with updated `other$loc.metrics` and
#'   `other$ind.metrics`.
#'
#' @section Known issue:
#' `loc_Fis` assumes a diploid. See `dev/audit.md`, finding M3.
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' gl = recalc_metrics(gl)
#' gl@other$loc.metrics
#' gl@other$ind.metrics
recalc_metrics <- function(gl){

  gl@other$loc.metrics <- data.frame(
    maf = get_maf(gl, max(adegenet::ploidy(gl))),
    loc_miss = get_loc_missing(gl),
    loc_het = get_loc_heterozygosity(gl, max(adegenet::ploidy(gl)))
  )
  
  # Use already calculated loc stats
  gl@other$loc.metrics$loc_Fis <- get_inbreeding(gl)
  
  gl@other$ind.metrics <- data.frame(
                                     ind_miss = get_ind_missing(gl),
                                     ind_het = get_ind_heterozygosity(gl, max(adegenet::ploidy(gl)))
                                     )
  return(gl)
}


#' Overall summary of a genlight object
#'
#' Summarises the metrics stored by [recalc_metrics()]. Run that function first
#' if the object did not come from a cgiarGenomics reader.
#'
#' @param gl A genlight object with `other$loc.metrics` and `other$ind.metrics`.
#'
#' @return A list with `nind` (individuals), `nloc` (loci), `ov_miss` (mean
#'   individual missingness), `ov_het` (mean individual heterozygosity) and
#'   `ov_maf` (mean MAF).
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' gl = recalc_metrics(gl)
#' get_overall_summary(gl)
get_overall_summary <- function(gl){
  ninds <- adegenet::nInd(gl)
  nlocs <- adegenet::nLoc(gl)
  overall_missiness <- mean(gl@other$ind.metrics$ind_miss)
  overall_heterozygosity <- mean(gl@other$ind.metrics$ind_het, na.rm = T)
  overall_maf <- mean(gl@other$loc.metrics$maf, na.rm = T)
  
  out <- list(
    nind = ninds,
    nloc = nlocs,
    ov_miss = overall_missiness,
    ov_het = overall_heterozygosity,
    ov_maf = overall_maf
  )
  return(out)
}



