# Majority rule function
#' Majority dosage of a vector of calls (internal)
#'
#' @param x Dosages for one locus across replicate samples.
#' @return The most frequent dosage, or `NA` when there is a tie or no data.
#' @noRd
majority_call <- function(x) {
  freq <- table(x)
  if (length(freq) == 0){
    return(NA)
  } else if (sum(freq == max(freq)) > 1) {
    # tie between calls (e.g. 0 vs 2): no majority, set missing instead of
    # silently picking the lowest dosage
    return(NA)
  } else {
    major <- names(freq)[which.max(freq)]
    return(as.numeric(major))
  }
}

#' Consensus row from a matrix of replicate samples (internal)
#'
#' @param mt Dosage matrix, replicates in rows and loci in columns.
#' @return A vector of majority calls, one per locus.
#' @noRd
collapse_inds <- function(mt) {
  col_list_split <- split(mt, col(mt))  
  consensus <- unlist(purrr::map(col_list_split, ~majority_call(.x)))
  return(consensus)
}

#' Merge replicate samples of the same genotype
#'
#' Collapses samples that belong to the same designation (germplasm) into one
#' individual, taking the majority dosage at each locus. Ties between dosages
#' become `NA` rather than picking one.
#'
#' @param gl A genlight object.
#' @param sample_dictionary A data frame with columns `sample_id` (matching
#'   `indNames(gl)`) and `designation_id` (the genotype each sample belongs to).
#'   Every individual in `gl` must appear in it.
#'
#' @return A genlight object with one individual per designation, named by
#'   `designation_id`. Designations with one sample come first, in the order
#'   of `gl`, followed by the merged designations in alphabetical order.
#'   Chromosomes, positions and alleles are kept and metrics are recalculated.
#'
#' @export
#'
#' @examples
#' mt = matrix(c(0, 0, 2,
#'               0, 2, 2,
#'               1, 1, 0),
#'             nrow = 3, byrow = TRUE,
#'             dimnames = list(c("s1", "s2", "s3"), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' dict = data.frame(sample_id = c("s1", "s2", "s3"),
#'                   designation_id = c("G1", "G1", "G2"))
#' merged = merge_duplicate_inds(gl, dict)
#' as.matrix(merged)
merge_duplicate_inds <- function(gl, sample_dictionary) {
  ind_names <- adegenet::indNames(gl)
  
  if(!all(c("sample_id", "designation_id") %in% colnames(sample_dictionary))) {
    cli::cli_abort("`sample_dictionary` needs the columns sample_id and designation_id")
  }
  if(sum(sample_dictionary$sample_id %in% ind_names) != adegenet::nInd(gl)) {
    cli::cli_abort("sample_dictionary doesn't consider all uploaded because not all priveded sample_id match")
  }
  
  # designation of each individual, in the order of the genlight object
  designation <- as.character(sample_dictionary$designation_id[match(ind_names, sample_dictionary$sample_id)])
  if(anyNA(designation)) {
    cli::cli_abort("Every individual in `gl` needs a designation_id in `sample_dictionary`")
  }
  
  mt <- as.matrix(gl)
  n_per_designation <- table(designation)
  single_designations <- unique(designation[designation %in% names(n_per_designation)[n_per_designation == 1]])
  dup_designations <- sort(names(n_per_designation)[n_per_designation > 1])
  
  # Single samples keep their calls and take the name of their designation
  s_mt <- mt[match(single_designations, designation), , drop = FALSE]
  rownames(s_mt) <- single_designations
  
  # Replicates are collapsed to the majority call at each locus
  dup_mt <- matrix(numeric(0), nrow = 0, ncol = ncol(mt))
  if(length(dup_designations) > 0) {
    dup_mt <- do.call(rbind, lapply(dup_designations, function(d) {
      collapse_inds(mt[which(designation == d), , drop = FALSE])
    }))
    rownames(dup_mt) <- dup_designations
  }
  
  out_mt <- rbind(s_mt, dup_mt)
  colnames(out_mt) <- adegenet::locNames(gl)
  
  # keep the source ploidy (otherwise genlight infers it per individual from max dosage)
  merged_gl <- new("genlight", out_mt,
                   ploidy = max(adegenet::ploidy(gl)),
                   chromosome = gl@chromosome,
                   position = gl@position)
  if (!is.null(adegenet::alleles(gl))) {
    adegenet::alleles(merged_gl) <- adegenet::alleles(gl)
  }
  merged_gl <- recalc_metrics(merged_gl)
  return(merged_gl)  
}




