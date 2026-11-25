#' Read a text table and use common na representations for convert to NA
#'
#' @param path Path to tabular data
#' @param sep  Delimitator of the table cells
#'
#' @return Returns a dataframe
#' @export
#'
#' @examples
read_tabular_geno <- function(path, sep = c('\t',',')){
  
  # Validate params
  sep = match.arg(sep)
  
  # Missing data representation values
  missingData = c("N","NN","FAIL","FAILED","Uncallable","Unused","-","NA","",-9)

  df <- as.data.frame(data.table::fread(path,
                                        sep = sep,
                                        header = TRUE,
                                        na.strings = missingData))
  return(df)
}

#' Split genotype calls into alleles (internal)
#' @noRd
split_strings <- function(l, ploidity){
  sapply(l, function(x){
    allele_length <- nchar(x) / ploidity
    # List with each allele as element
    split_genotype <- substring(x, seq(1, nchar(x), allele_length),
                                seq(allele_length, nchar(x), allele_length))
    return(split_genotype)
  })
}



#' Replay BioFlow modifications on a genlight object
#'
#' BioFlow records quality control as a table of modifications instead of
#' saving a new object. This function replays that table. Rows with `reason`
#' starting with `imputation` are imputations. All other rows are removals,
#' where a missing `col` removes the individual in `row` and a missing `row`
#' removes the locus in `col`.
#'
#' @param gl A genlight object, the same one BioFlow used to record the
#'   modifications.
#' @param modifications A data frame with columns `reason`, `row`, `col` and,
#'   for imputations, `value`.
#'
#' @return The genlight object with the removals applied.
#'
#' @section Known issue:
#' Imputation rows are never applied, because the function returns before
#' reaching them. See `dev/audit.md`, finding B2.
#' @export
apply_bioflow_modifications <- function(gl, modifications){
  
  # Filtering modifications
  filt_mods <- modifications %>% 
    dplyr::filter(!grepl("^imputation", reason))
  
  if (dim(filt_mods)[1] == 0) {
    return(gl)
  }
  
  
  ind_out <- filt_mods %>% 
    dplyr::filter(is.na(col)) %>% 
    dplyr::pull(row)
  
  loc_out <- filt_mods %>% 
    dplyr::filter(is.na(row)) %>% 
    dplyr::pull(col)
  
  # setdiff() keeps every element when nothing is removed from a margin.
  # Negative indexing with an empty vector would drop the whole margin.
  ind_idx <- setdiff(seq_len(adegenet::nInd(gl)), ind_out)
  loc_idx <- setdiff(seq_len(adegenet::nLoc(gl)), loc_out)
  

  filt_gl <- gl[ind_idx, loc_idx]
  return(filt_gl)
  
  imp_mods <- modifications %>% 
    filter(grepl("^imputation", reason))
  
  mt <- as.matrix(filt_gl)

  mt[cbind(imp_mods$row, imp_mods$col)] <- imp_mods$value
  

  imp_gl <- new("genlight",
                mt,
                loc.names = filt_gl@loc.names,
                ind.names = filt_gl@ind.names,
                chromosome = filt_gl@chromosome,
                position = filt_gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(filt_gl)
  imp_gl <- recalc_metrics(imp_gl)
  
  return(imp_gl)
}



#' Format counts with k, M and B suffixes for plot axes (internal)
#' @noRd
addUnits <-
  function(n) {
    labels <-
      ifelse(n < 1000, n,  # less than thousands
             ifelse(n < 1e6, paste0(round(n/1e3), 'k'),  # in thousands
                    ifelse(n < 1e9, paste0(round(n/1e6), 'M'),  # in millions
                           ifelse(n < 1e12, paste0(round(n/1e9), 'B'), # in billions
                                  'too big!'))))
    return(labels)
  }

#' Midpoint of a cut() interval label (internal)
#' @noRd
get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}

#' Write genlight object to VCF file
#'
#' This function converts a genlight object to a VCF file and writes it to an output path.
#'
#' @param gl A genlight object.
#' @param output_path String. Path to the output VCF file. It could end in `.gz` for compression.
#' @param na_rep String. NA representation for genotypes in the VCF file. Default is "./.".
#'
#' @export
write_vcf <- function(gl, output_path, na_rep = "./.") {
  if (!inherits(gl, "genlight")) {
    cli::cli_abort("`gl` must be a genlight object.")
  }
  if (!is.character(output_path) || length(output_path) != 1) {
    cli::cli_abort("`output_path` must be a single character string.")
  }
  
  n_loci <- length(adegenet::locNames(gl))
  n_samples <- length(adegenet::indNames(gl))
  
  gt_matrix <- as.matrix(gl)
  
  # Extract alleles
  alleles <- adegenet::alleles(gl)
  if (is.null(alleles)) {
    ref_alleles <- rep("N", n_loci)
    alt_alleles <- rep("N", n_loci)
  } else {
    split_alleles <- strsplit(as.character(alleles), "/")
    ref_alleles <- sapply(split_alleles, function(x) {
      if (length(x) >= 1 && !is.na(x[1])) x[1] else "N"
    })
    alt_alleles <- sapply(split_alleles, function(x) {
      if (length(x) >= 2 && !is.na(x[2])) x[2] else "N"
    })
  }
  
  # Ploidy
  p_ind <- adegenet::ploidy(gl)
  if (length(p_ind) == 1) {
    p_ind <- rep(p_ind, n_samples)
  }
  
  # Translate dosages to VCF genotype format
  V <- matrix("", nrow = n_samples, ncol = n_loci)
  for (i in seq_len(n_samples)) {
    p <- p_ind[i]
    possible_dosages <- 0:p
    gt_strings <- sapply(possible_dosages, function(d) {
      paste(c(rep("0", p - d), rep("1", d)), collapse = "/")
    })
    names(gt_strings) <- as.character(possible_dosages)
    
    row_na_rep <- if (is.null(na_rep)) {
      paste(rep(".", p), collapse = "/")
    } else {
      na_rep
    }
    
    row_vals <- gt_matrix[i, ]
    row_chars <- as.character(row_vals)
    translated <- gt_strings[row_chars]
    translated[is.na(row_vals)] <- row_na_rep
    V[i, ] <- translated
  }
  
  # Calculate AC and AN
  ac_vals <- colSums(gt_matrix, na.rm = TRUE)
  if (length(unique(p_ind)) == 1) {
    an_vals <- colSums(!is.na(gt_matrix)) * p_ind[1]
  } else {
    an_vals <- colSums(!is.na(gt_matrix) * p_ind)
  }
  info_val <- paste0("AC=", ac_vals, ";AN=", an_vals)
  
  # Chromosome and position
  chrom_val <- if (!is.null(adegenet::chromosome(gl))) as.character(adegenet::chromosome(gl)) else rep(".", n_loci)
  pos_val <- if (!is.null(adegenet::position(gl))) as.character(adegenet::position(gl)) else rep(".", n_loci)
  
  chrom_val[is.na(chrom_val)] <- "."
  pos_val[is.na(pos_val)] <- "."
  
  id_val <- adegenet::locNames(gl)
  if (is.null(id_val)) id_val <- paste0("Locus_", seq_len(n_loci))
  
  # VCF data columns
  vcf_data <- data.frame(
    `#CHROM` = chrom_val,
    POS = pos_val,
    ID = id_val,
    REF = ref_alleles,
    ALT = alt_alleles,
    QUAL = rep(".", n_loci),
    FILTER = rep("PASS", n_loci),
    INFO = info_val,
    FORMAT = rep("GT", n_loci),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Genotype columns
  gt_cols <- t(V)
  colnames(gt_cols) <- adegenet::indNames(gl)
  vcf_df <- cbind(vcf_data, gt_cols)
  
  # Headers
  header <- c(
    "##fileformat=VCFv4.0",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
  )
  
  # Write
  con <- if (grepl("\\.gz$", output_path, ignore.case = TRUE)) {
    gzfile(output_path, "wt")
  } else {
    file(output_path, "wt")
  }
  writeLines(header, con = con)
  write.table(vcf_df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  close(con)
}

#' Validate method-specific arguments passed through ... (internal)
#'
#' @param dots Named list of arguments.
#' @param allowed,required Character vectors of argument names.
#' @param method Method name, used in error messages.
#' @noRd
check_method_args <- function(dots, allowed, required, method) {
  supplied <- names(dots)
  
  if (is.null(supplied)) {
    supplied <- character()
  }
  
  unnamed <- supplied == ""
  
  if (any(unnamed)) {
    stop(
      "All method-specific arguments passed through '...' must be named.",
      call. = FALSE
    )
  }
  
  unknown <- setdiff(supplied, allowed)
  
  if (length(unknown) > 0) {
    stop(
      sprintf(
        "Unknown argument(s) for method '%s': %s. Allowed arguments are: %s.",
        method,
        paste(unknown, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  missing_required <- setdiff(required, supplied)
  
  if (length(missing_required) > 0) {
    stop(
      sprintf(
        "Missing required argument(s) for method '%s': %s.",
        method,
        paste(missing_required, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

#' Count Mendelian inconsistencies in parent and offspring trios
#'
#' For each trio, counts loci where both parents carry the same homozygous
#' dosage (0 or 2) and the offspring differs from it. Only diploid 0/1/2
#' dosages are supported.
#'
#' @param gl A diploid genlight object, ideally with no missing data.
#' @param trios A data frame with columns `designation` (offspring), `mother`
#'   and `father`, holding individual names of `gl`.
#'
#' @return The `trios` data frame with added columns `pedigree_present` (all
#'   three individuals found in `gl`), `n_markers`, `n_evaluated` (loci typed
#'   in all three), `n_missing`, `n_errors` and `error_rate`.
#'
#' @section Known issue:
#' Opposite homozygous parents (0 and 2) with a homozygous offspring are not
#' counted as errors, so the count is a lower bound. See `dev/audit.md`,
#' finding P2.
#' @export
#'
#' @examples
#' mt = matrix(c(0, 0, 2,
#'               2, 2, 2,
#'               1, 0, 1),
#'             nrow = 3, byrow = TRUE,
#'             dimnames = list(c("mum", "dad", "kid"), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' trios = data.frame(designation = "kid", mother = "mum", father = "dad")
#' get_mendelian_error_counts(gl, trios)
get_mendelian_error_counts <- function(gl, trios) {
  
  if (any(adegenet::ploidy(gl) != 2)) {
    cli::cli_abort("Only diploid 0/1/2 dosages are supported.")
  }
  # # invalid trios for diploid late generations
  #   0, 2, 2,
  #   2, 0, 0,
  #   1, 0, 0,
  #   1, 2, 2
  # 
  # valid trios for F1
  # 0, 0, 0,
  # 0, 0, 1,
  # 1, 0, 1,
  # 1, 0, 2,
  # 0, 1, 0,
  # 1, 1, 0,
  # 0, 1, 1,
  # 1, 1, 1,
  # 2, 1, 1,
  # 1, 1, 2,
  # 2, 1, 2,
  # 1, 2, 0,
  # 1, 2, 1,
  # 2, 2, 1,
  # 2, 2, 2
  
  G <- as.matrix(gl)
  ids <- adegenet::indNames(gl)
  
  child_idx <- match(trios$designation, ids)
  mother_idx <- match(trios$mother, ids)
  father_idx <- match(trios$father, ids)
  
  present <- !is.na(child_idx) & !is.na(mother_idx) & !is.na(father_idx)
  
  n <- nrow(trios)
  
  res <- data.frame(n_markers = rep(NA_integer_, n), n_evaluated = NA_integer_,
                    n_missing = NA_integer_, n_errors = NA_integer_,
                    error_rate = NA_real_)
  
  if (any(present)) {
    C  <- G[child_idx[present], , drop = FALSE]
    M  <- G[mother_idx[present], , drop = FALSE]
    Fa <- G[father_idx[present], , drop = FALSE]
    ok  <- !is.na(C) & !is.na(M) & !is.na(Fa)
    # late generations: impossible only if parents share a homozygote and child differs
    err <- ok & (M == Fa) & (M != 1) & (C != M)
    n_ev  <- rowSums(ok)
    n_err <- rowSums(err)
    res[present, ] <- data.frame(ncol(G), n_ev, ncol(G) - n_ev, n_err,
                                 ifelse(n_ev > 0, n_err / n_ev, NA_real_))
  }
  
  out <- trios
  out$pedigree_present <- present
  cbind(out, res)
}




