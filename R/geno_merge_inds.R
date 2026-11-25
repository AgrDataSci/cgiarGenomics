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
#' @return A genlight object with one individual per designation.
#'
#' @section Known issue:
#' Single samples are matched on `sample_id` instead of `designation_id`, the
#' function fails when no designation has replicates, and the result loses
#' chromosomes, positions, alleles and metrics. Run [recalc_metrics()] on the
#' result before filtering it. See `dev/audit.md`, finding P1.
#' @export
#'
#' @examples
#' \dontrun{
#' dict = data.frame(sample_id = c("s1", "s2", "s3"),
#'                   designation_id = c("G1", "G1", "G2"))
#' merged = merge_duplicate_inds(gl, dict)
#' }
merge_duplicate_inds <- function(gl, sample_dictionary) {
  if(sum(sample_dictionary$sample_id %in% indNames(gl))!=nInd(gl)) {
    cli::cli_abort("sample_dictionary doesn't consider all uploaded because not all priveded sample_id match")
  }
  
  single_samples <- sample_dictionary %>% 
    dplyr::group_by(designation_id) %>% 
    dplyr::filter(dplyr::n() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::pull(designation_id)


  dup_samples <- sample_dictionary[which(!sample_dictionary$sample_id %in% single_samples),]

  dup_mt <- dup_samples %>% 
    split(.$designation_id) %>% 
    purrr::map(~{
      # rows of gl that belong to this designation (not positions inside .x)
      ind_idx <- which(indNames(gl) %in% .x$sample_id)
      smt <- as.matrix(gl[ind_idx,])
      collapse_inds(smt)
    }) %>% 
    do.call(rbind, .)
  # split() orders groups by factor levels, so name rows the same way
  rownames(dup_mt) <- levels(factor(dup_samples$designation_id))
  
  
  if(length(single_samples) > 0) {
    ss_idx <- which(indNames(gl) %in% single_samples)
    s_mt <- as.matrix(gl[ss_idx,])
    out_mt <- rbind(s_mt, dup_mt)
  } else {
    out_mt <- dup_mt
  }
  colnames(out_mt) <- locNames(gl)
  # keep the source ploidy (otherwise genlight infers it per individual from max dosage)
  merged_gl <- new("genlight", out_mt, ploidy = max(adegenet::ploidy(gl)))
  return(merged_gl)  
}




