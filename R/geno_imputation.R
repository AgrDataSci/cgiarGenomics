#' Draw one genotype call from an allele frequency
#'
#' Draws each of the `ploidity` allele copies independently, with probability
#' `q_frq` of being the alternative allele, and returns their sum.
#'
#' @param q_frq Numeric between 0 and 1. Frequency of the alternative allele.
#' @param ploidity Integer. Ploidy level (default 2).
#'
#' @return Numeric. A dosage between 0 and `ploidity`.
#' @keywords internal
#'
#' @examples
#' set.seed(1)
#' cgiarGenomics:::i_freq_impute(0.3, ploidity = 2)
i_freq_impute = function(q_frq, ploidity = 2){
  if(!is_whole_number(ploidity)){
    cli::cli_abort("`ploidity` is not an integer: {ploidity}")  
  }
  dosage = 0
  for(i_chromatid in seq(ploidity)){
    i_dosage = sample(c(1,0), size = 1, 
                       prob = c(q_frq, 1 - q_frq), replace = T)
    dosage = dosage + i_dosage
  }
  return(dosage)
}

#' Frequency imputation of all missing calls (internal)
#'
#' @return A list of linear matrix indices, named by imputed dosage.
#' @noRd
freq_impute = function(gl, mt, ploidity){
  mt = as.matrix(gl)
  # Linear index of nas
  idx_na = which(is.na(mt))
  # A complete panel has nothing to impute
  if (length(idx_na) == 0) {
    return(list())
  }
  # Get the allelic frequencies
  q_allele = adegenet::glMean(gl)
  # Locus (column) of each missing cell
  na_loc_idx = ceiling(idx_na/nrow(mt))
  
  imp = vapply(q_allele[na_loc_idx],
                function(x) as.numeric(i_freq_impute(q_frq = x, ploidity)),
                numeric(1))
  return(split(idx_na, unname(imp)))
}

#' Write imputed dosages into a matrix (internal)
#'
#' @param mt Dosage matrix.
#' @param imp_dict List of linear indices named by dosage.
#' @return The matrix with imputed values.
#' @noRd
apply_imputation = function(mt, imp_dict){
  for (dosage in names(imp_dict)) {
    # Convert the list name to a numeric value
    num_dosage = as.numeric(dosage)
    
    # Get the linear indices associated with this value
    idx = imp_dict[[dosage]]
    
    # Assign the value to these positions in the matrix
    mt[idx] = num_dosage
  }
  return(mt)
}

#' Random Forest imputation for a single marker
#'
#' Trains a random forest model on a target marker using flanking markers as predictors.
#' For missing genotypes in the target marker, predictions are made from the trained model.
#'
#' @param target Vector of genotype calls for target marker (may contain NA)
#' @param predictors Matrix of predictor markers (flanking markers, rows=samples, cols=markers)
#' @param ntree Number of trees in random forest (default=100)
#' @param seed Optional seed for reproducibility
#'
#' @return Vector of imputed genotype calls for target marker
#' @keywords internal
i_rf_impute = function(target, predictors, ntree = 100, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Identify missing and non-missing positions in target
  miss_idx = which(is.na(target))
  n_miss = length(miss_idx)
  
  # If no missing data, return as-is
  if (n_miss == 0) {
    return(NA)
  }
  
  # If no predictors available, return unchanged
  if (ncol(predictors) == 0) {
    return(NA)
  }
  
  # Find indices where target is not NA AND predictors are complete (for training)
  pred_complete = which(complete.cases(predictors))
  train_idx = pred_complete[!is.na(target[pred_complete])]
  

  if (length(train_idx) < 2) {
    # Not enough training data
    return(NA)
  }
  
  
  # Filter target and predictors to training indices
  target_train = target[train_idx]
  
  
  if (length(unique(target_train)) == 1) {
    # Target snp to be predicted only have one allele present. Use freq impute predicition
    return(NA)
  }
  
  predictors_train = predictors[train_idx, , drop = FALSE]
  
  # Train random forest on non-missing observations
  rf_model = randomForest::randomForest(
    x = predictors_train,
    y = target_train,
    ntree = ntree,
    na.action = na.omit
  )

  
  if (length(miss_idx) > 0) {
    # Predict for samples with complete predictor data
    predictors_miss_complete = predictors[miss_idx, , drop = FALSE]
    predictions = stats::predict(rf_model, newdata = predictors_miss_complete)
    
    # Round predictions to nearest integer (genotype call)
    # predictions = round(predictions)
    
    # Assign predictions only to positions with complete predictors
    target[miss_idx] = predictions
  }
  return(target)
}

#' Random Forest imputation for genotype matrix
#'
#' Applies Random Forest imputation to a genotype matrix by iterating through markers.
#' For each marker with missing data, uses flanking markers as predictors.
#' Respects chromosome boundaries to avoid using markers from different chromosomes.
#' Markers are processed in parallel using the specified number of cores.
#' Predictors are first filled by frequency imputation.
#'
#' @param gl genlight object
#' @param ploidity Ploidy level (default=2)
#' @param nflank Number of flanking markers on each side (default=100)
#' @param ntree Number of trees in random forest (default=100)
#' @param seed Optional seed for reproducibility
#' @param nthreads Number of parallel threads. Default uses get_default_beagle_threads() (half cores, max 4).
#'
#' @return List where names are dosage values and elements are linear indices of NA positions
#'
#' @keywords internal
rf_impute = function(gl, ploidity = 2, nflank = 100, ntree = 100, seed = NULL, nthreads = NULL) {
  
  require_optional("randomForest")
  
  if (is.null(nthreads)) {
    nthreads = get_default_beagle_threads()
  }
  
  # The predictors are filled by random draws, so the seed is set here too.
  # Otherwise two runs with the same seed could differ.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # impute with frequency to get a filled predictors matrix
  pred_gl = impute_gl(gl, ploidity, method = "frequency")$gl
  pred_mt = as.matrix(pred_gl)
  
  mt = as.matrix(gl)
  n_markers = ncol(mt)
  n_samples = nrow(mt)
  
  # Get chromosome information for each marker
  chr = adegenet::chromosome(gl)
  
  # Identify markers that need imputation
  markers_with_na = which(colSums(is.na(mt)) > 0)
  n_to_impute = length(markers_with_na)
  
  if (n_to_impute == 0) {
    return(list())
  }
  
  cli::cli_inform("Imputing {n_to_impute} markers with missing data using {nthreads} thread(s)...")
  
  # Worker function: imputes a single marker and returns the imputed column

  impute_one_marker = function(i_marker) {
    target = mt[, i_marker]
    
    # Get chromosome of target marker
    target_chr = chr[i_marker]
    
    # Find all markers on the same chromosome
    same_chr_idx = which(chr == target_chr)
    
    # Find position of target marker within chromosome markers
    pos_in_chr = which(same_chr_idx == i_marker)
    
    # Skip if marker not found in this chromosome
    if (length(pos_in_chr) == 0) {
      return(pred_mt[, i_marker])
    }
    
    # Determine flanking range respecting chromosome boundaries
    n_available_before = pos_in_chr - 1
    n_available_after = length(same_chr_idx) - pos_in_chr
    
    n_flank_before = min(nflank, n_available_before)
    n_flank_after = min(nflank, n_available_after)
    
    # If one side has fewer than nflank, compensate from the other side
    if (n_flank_before < nflank) {
      n_flank_after = min(nflank + (nflank - n_flank_before), n_available_after)
    } else if (n_flank_after < nflank) {
      n_flank_before = min(nflank + (nflank - n_flank_after), n_available_before)
    }
    
    # Get indices of flanking markers
    flank_idx_before = integer(0)
    flank_idx_after = integer(0)
    
    if (n_flank_before > 0 && pos_in_chr > 1) {
      start_pos = max(1, pos_in_chr - n_flank_before)
      end_pos = pos_in_chr - 1
      if (start_pos <= end_pos) {
        flank_idx_before = same_chr_idx[start_pos:end_pos]
      }
    }
    
    if (n_flank_after > 0 && pos_in_chr < length(same_chr_idx)) {
      start_pos = pos_in_chr + 1
      end_pos = min(length(same_chr_idx), pos_in_chr + n_flank_after)
      if (start_pos <= end_pos) {
        flank_idx_after = same_chr_idx[start_pos:end_pos]
      }
    }
    
    flank_idx = c(flank_idx_before, flank_idx_after)
    flank_idx = flank_idx[!is.na(flank_idx) & flank_idx > 0 & flank_idx <= n_markers]
    
    # Get predictor markers
    if (length(flank_idx) > 0) {
      predictors = pred_mt[, flank_idx, drop = FALSE]
    } else {
      predictors = matrix(numeric(0), nrow = n_samples, ncol = 0)
    }
    
    # Impute using RF
    imputed_target = i_rf_impute(
      target = as.factor(unlist(target)),
      predictors = predictors,
      ntree = ntree,
      seed = seed
    )
    
    if (length(imputed_target) == 1) {
      return(pred_mt[, i_marker])
    }
    # randomForest classifies and returns a factor. Convert its labels back to
    # dosages, otherwise the matrix receives the factor codes.
    return(as.numeric(as.character(imputed_target)))
  }
  
  # Run in parallel or sequentially depending on nthreads
  if (nthreads > 1) {
    if (.Platform$OS.type == "windows") {
      # Windows: use parLapply with a PSOCK cluster
      cl = parallel::makeCluster(nthreads)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      # Export required objects and functions to the cluster
      parallel::clusterExport(cl, varlist = c(
        "mt", "pred_mt", "chr", "n_markers", "n_samples",
        "nflank", "ntree", "seed", "i_rf_impute"
      ), envir = environment())
      parallel::clusterEvalQ(cl, library(randomForest))
      
      results = parallel::parLapply(cl, markers_with_na, impute_one_marker)
    } else {
      # Unix: use mclapply (fork-based, no export needed)
      results = parallel::mclapply(
        markers_with_na, impute_one_marker,
        mc.cores = nthreads
      )
    }
  } else {
    results = lapply(markers_with_na, impute_one_marker)
  }
  
  # Apply results back to the matrix
  for (i in seq_along(markers_with_na)) {
    mt[, markers_with_na[i]] = results[[i]]
  }
  
  # Build imputation dictionary: group by imputed dosage values
  idx_na = which(is.na(as.matrix(gl)))
  if (length(idx_na) > 0) {
    imputed_values = mt[idx_na]
    imp_dict = split(idx_na, imputed_values)
  }
  
  return(imp_dict)
}


#' Impute a gl object
#'
#' Fills missing dosages with one of three methods and returns a new genlight
#' object with recalculated metrics. Loci with no data are removed first, with
#' a warning. An object with no missing calls is returned with an empty log.
#' The function stops when a method returns a value outside the dosage range.
#'
#' * `"frequency"` draws each missing call from the allele frequency of its
#'   locus with [i_freq_impute()]. Fast, ignores linkage.
#' * `"random_forest"` predicts each missing call from flanking markers on the
#'   same chromosome. Extra arguments `nflank`, `ntree`, `seed`, `nthreads`.
#'   Needs the suggested package randomForest.
#' * `"beagle"` runs Beagle 5 on a temporary VCF, diploids only, see
#'   [impute_beagle()]. Extra arguments `jre_path`, `beagle_path`, `memory`,
#'   `burnin`, `iterations`, `seed`, `nthreads`.
#'
#' @param gl genlight object
#' @param ploidity Integer. Ploidy level (default 2). Must match the object.
#' @param method Character. `"frequency"`, `"random_forest"` or `"beagle"`.
#'   There is no default.
#' @param ... Named, method-specific arguments listed above. Unknown names are
#'   an error.
#'
#' @return List with elements:
#'   - gl: imputed genlight object
#'   - log: imputation dictionary (imputed positions grouped by dosage)
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
#' set.seed(1)
#' imp = impute_gl(gl, ploidity = 2, method = "frequency")
#' as.matrix(imp$gl)
impute_gl = function(gl, ploidity = 2, method, ...){
  
  method = match.arg(
    method,
    choices = c("frequency", "random_forest", "beagle")
  )
  dots = list(...)
  
  obj_ploidy = max(adegenet::ploidy(gl))
  if (ploidity != obj_ploidy) {
    cli::cli_warn("`ploidity` is {ploidity} but the object has ploidy {obj_ploidy}. Pass the ploidy of the data.")
  }
  
  mt = as.matrix(gl)
  loci_all_nas = colSums(!is.na(mt)) == 0
  
  if(sum(loci_all_nas) > 0){
    cli::cli_warn("There are {sum(loci_all_nas)} loci with all missing data. They are removed before imputation.")
    # Filter out the all na loci
    all_notna_idxs = which(!loci_all_nas)
    gl = gl[,all_notna_idxs]
    mt = as.matrix(gl)
  }
  
  number_imputations = sum(is.na(mt))
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  # A complete panel is a valid input, for instance a cleaned panel submitted
  # twice. Return it with an empty log.
  if (number_imputations == 0) {
    return(list(gl = recalc_metrics(gl), log = list()))
  }
  
  if(method == 'frequency'){
    imp_dict = freq_impute(gl, mt, ploidity)
  } else if(method == 'random_forest'){
    allowed = c("nflank", 'ntree', 'seed', 'nthreads')
    check_method_args(
      dots = dots,
      allowed = allowed,
      required = character(),
      method = method
    )
    cli::cli_inform("Imputing with Random Forest")
    imp_dict = do.call(rf_impute, c(list(gl = gl,
                                     ploidity = ploidity), dots))
      
  } else if(method == 'beagle'){
    allowed = c("jre_path", 'beagle_path', 'memory',
                 'burnin', 'iterations', 'seed', 'nthreads')
    check_method_args(
      dots = dots,
      allowed = allowed,
      required = character(),
      method = method
    )
    imp_dict = do.call(impute_beagle, c(list(gl = gl), dots))
  }
  
  # Every imputed value must be a dosage, and only missing cells may change
  imp_values = suppressWarnings(as.numeric(names(imp_dict)))
  if (anyNA(imp_values) || any(imp_values < 0 | imp_values > ploidity | imp_values %% 1 != 0)) {
    cli::cli_abort("Imputation returned values outside the dosage range 0 to {ploidity}: {names(imp_dict)}")
  }
  imputed_cells = unlist(imp_dict, use.names = FALSE)
  if (!all(imputed_cells %in% which(is.na(mt)))) {
    cli::cli_abort("Imputation tried to change calls that were not missing.")
  }
  
  # apply the imputation creating a new gl instance
  imp_mt = apply_imputation(mt, imp_dict)
  
  
  imp_gl = new("genlight",
            imp_mt,
            ploidy = ploidity,
            loc.names = gl@loc.names,
            ind.names = gl@ind.names,
            chromosome = gl@chromosome,
            position = gl@position)
  
  if (!is.null(adegenet::alleles(gl))) {
    adegenet::alleles(imp_gl) = adegenet::alleles(gl)
  }
  imp_gl = recalc_metrics(imp_gl)
  return(list(gl = imp_gl, log = imp_dict))
}


#' Mask a fraction of genotype calls for measure the
#' accuracy of an imputation method
#'
#' @param gl A genlight instance
#' @param fraction Share of calls to mask, between 0 and 1 (default 0.1)
#'
#' @details
#' Only calls that are not missing are masked. The result keeps individual
#' and locus names, chromosomes, positions and alleles.
#'
#' @return A genlight instance with the given fraction of calls set to `NA`.
#'   The linear indices of the masked cells are in `other$masked_cells`. Pass
#'   the object to [get_accuracy()] as `masked_gl`.
#'
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, 1,
#'               2, 2, 1, 0,
#'               0, 1, 1, 2),
#'             nrow = 4)
#' gl = new("genlight", mt, ploidy = 2)
#' set.seed(1)
#' masked = set_fraction_na(gl, fraction = 0.25)
#' as.matrix(masked)
set_fraction_na = function(gl, fraction = 0.1) {
    # Ensure fraction is between 0 and 1
    if (!is.numeric(fraction) || length(fraction) != 1 || fraction < 0 || fraction > 1) {
      stop("fraction must be a number between 0 and 1")
    }
    m = as.matrix(gl)
    
    # Only calls that are not missing can be masked
    called = which(!is.na(m))
    
    # Compute how many elements to replace
    n_replace = min(round(length(m) * fraction), length(called))
    
    # Randomly choose positions
    positions = integer(0)
    if (n_replace > 0) {
      positions = sort(called[sample.int(length(called), size = n_replace)])
    }
    
    # Replace chosen positions with NA
    m[positions] = NA
    
    masked_gl = new("genlight", m, ploidy = max(adegenet::ploidy(gl)),
                     chromosome = gl@chromosome,
                     position = gl@position)
    if (!is.null(adegenet::alleles(gl))) {
      adegenet::alleles(masked_gl) = adegenet::alleles(gl)
    }
    masked_gl = recalc_metrics(masked_gl)
    # linear indices of the masked cells, for get_accuracy()
    masked_gl@other$masked_cells = positions
    return(masked_gl)
}


#' Measure imputation accuracy on masked calls
#'
#' Compares an imputed genlight object with the reference it was masked from.
#' With `masked_gl`, only the cells that [set_fraction_na()] masked are scored,
#' which is the accuracy of the imputation. A masked cell that stays missing
#' counts as wrong. Objects are matched by individual and locus names, so loci
#' removed by [impute_gl()] are left out.
#'
#' @param ref_gl Reference genlight instance, before masking.
#' @param imp_gl Imputed genlight instance.
#' @param masked_gl The masked genlight instance that was imputed, as returned
#'   by [set_fraction_na()]. When `NULL`, every called reference cell is
#'   scored, including cells that were never masked. That inflates accuracy,
#'   and the function warns.
#'
#' @return A numeric vector with one value per locus, the share of scored
#'   cells that the imputed object reproduces. Loci with no scored cell are
#'   `NA`. The attribute `n_scored` holds the number of scored cells per locus.
#'
#' @export
#'
#' @examples
#' mt = matrix(rep(c(0, 1, 2, 1), 6), nrow = 8,
#'             dimnames = list(paste0("ind", 1:8), paste0("snp", 1:3)))
#' gl = new("genlight", mt, ploidy = 2)
#' set.seed(1)
#' masked = set_fraction_na(gl, fraction = 0.25)
#' imputed = impute_gl(masked, ploidity = 2, method = "frequency")
#' get_accuracy(gl, imputed$gl, masked_gl = masked)
get_accuracy = function(ref_gl, imp_gl, masked_gl = NULL){
  # Get the gt matix of the two gls
  ref_m = as.matrix(ref_gl)
  imp_m = as.matrix(imp_gl)
  
  # Imputation removes loci with no data, so match the two by name
  has_names = function(m) !is.null(rownames(m)) && !is.null(colnames(m))
  if (has_names(ref_m) && has_names(imp_m)) {
    inds = intersect(rownames(ref_m), rownames(imp_m))
    locs = intersect(colnames(ref_m), colnames(imp_m))
    ref_m = ref_m[inds, locs, drop = FALSE]
    imp_m = imp_m[inds, locs, drop = FALSE]
  }
  if (!identical(dim(ref_m), dim(imp_m))) {
    cli::cli_abort("`ref_gl` and `imp_gl` do not have the same individuals and loci.")
  }
  
  if (is.null(masked_gl)) {
    cli::cli_warn(c(
      "Every called reference cell is scored, including cells that were never masked, so accuracy is inflated.",
      "i" = "Pass the output of `set_fraction_na()` as `masked_gl` to score the masked cells only."
    ))
    target = !is.na(ref_m)
  } else {
    mask_m = as.matrix(masked_gl)
    if (has_names(ref_m) && has_names(mask_m)) {
      if (!all(rownames(ref_m) %in% rownames(mask_m)) || !all(colnames(ref_m) %in% colnames(mask_m))) {
        cli::cli_abort("`masked_gl` lacks individuals or loci of `ref_gl`.")
      }
      mask_m = mask_m[rownames(ref_m), colnames(ref_m), drop = FALSE]
    }
    if (!identical(dim(ref_m), dim(mask_m))) {
      cli::cli_abort("`ref_gl` and `masked_gl` do not have the same individuals and loci.")
    }
    # Cells with a known answer that the imputation had to fill
    target = is.na(mask_m) & !is.na(ref_m)
  }
  
  # A cell left missing by the imputation counts as wrong
  hit = target & !is.na(imp_m) & (ref_m == imp_m)
  n_scored = colSums(target)
  accuracy = colSums(hit, na.rm = TRUE)/n_scored
  accuracy[n_scored == 0] = NA
  attr(accuracy, "n_scored") = n_scored
  
  return(accuracy)
}

#' Impute missing genotypes using Beagle
#'
#' This function converts a genlight object to a VCF file, runs Beagle for imputation,
#' reads the imputed VCF back, and returns the imputation dictionary.
#'
#' Java and the Beagle JAR are auto-detected by default:
#' - Java is found via JAVA_HOME or system PATH
#' - Beagle JAR is found in the package's inst/extdata/ directory
#'
#' Both can be overridden manually if needed.
#'
#' @param gl A genlight object of exclusively diploid ploidity level.
#' @param jre_path String. Path to Java executable or JRE directory. 
#'   If NULL (default), auto-detected via find_java().
#' @param beagle_path String. Path to the Beagle JAR file. 
#'   If NULL (default), auto-detected via find_beagle_jar().
#' @param memory String. Java memory setting, such as `"Xmx4g"` or
#'   `"Xmx512m"`. Any other text is refused.
#' @param burnin Integer. Number of burn-in iterations. Default is 3.
#' @param iterations Integer. Number of phasing iterations. Default is 12.
#' @param seed Integer. Seed for the random number generator. Default is -99999.
#' @param nthreads Integer. Number of threads to use. Default uses get_default_beagle_threads().
#'
#' @return A list where names are dosage values and elements are linear indices of imputed positions.
#' @export
impute_beagle = function(gl, jre_path = NULL, beagle_path = NULL, memory = "Xmx4g",
                          burnin = 3, iterations = 12, seed = -99999, nthreads = NULL) {
  
  # These settings reach the command line. Check them before anything runs,
  # so that only a Java memory option and whole numbers can be passed.
  if (!is.character(memory) || length(memory) != 1 ||
      !grepl("^-?Xm[sx][0-9]+[kKmMgG]?$", memory)) {
    cli::cli_abort(c(
      "`memory` is not a valid Java memory setting.",
      "i" = "Use a value such as \"Xmx4g\" or \"Xmx512m\"."
    ))
  }
  whole_numbers = list(burnin = burnin, iterations = iterations, seed = seed)
  if (!is.null(nthreads)) {
    whole_numbers$nthreads = nthreads
  }
  for (arg_name in names(whole_numbers)) {
    value = whole_numbers[[arg_name]]
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value != round(value)) {
      cli::cli_abort("`{arg_name}` must be a single whole number.")
    }
  }
  if (burnin < 0 || iterations < 0 || (!is.null(nthreads) && nthreads < 1)) {
    cli::cli_abort("`burnin` and `iterations` cannot be negative, and `nthreads` must be at least 1.")
  }
  
  # Beagle round-trips through VCF, so vcfR is needed. Check up front rather
  # than after a potentially long Beagle run.
  require_optional("vcfR")
  
  if (!inherits(gl, "genlight")) {
    cli::cli_abort("`gl` must be a genlight object.")
  }
  
  p_ind = adegenet::ploidy(gl)
  if (any(p_ind != 2)) {
    cli::cli_abort("The genlight object must be exclusively diploid.")
  }

  # Auto-detect Beagle JAR if not provided
  if (is.null(beagle_path)) {
    beagle_path = find_beagle_jar()
    if (is.null(beagle_path)) {
      cli::cli_abort(c(
        "Beagle JAR not found.",
        "i" = "Place the Beagle .jar file in: {system.file('extdata', package = 'cgiarGenomics')}",
        "i" = "Or provide the path explicitly via the `beagle_path` argument."
      ))
    }
    cli::cli_inform("Using bundled Beagle JAR: {basename(beagle_path)}")
  } else {
    if (!file.exists(beagle_path)) {
      cli::cli_abort("`beagle_path` does not exist: {beagle_path}")
    }
  }
  
  # Auto-detect Java if not provided
  if (is.null(jre_path)) {
    java_exe = find_java()
    if (is.null(java_exe)) {
      cli::cli_abort(c(
        "Java not found on the system.",
        "i" = "Install Java 8+ from https://adoptium.net/temurin/releases/",
        "i" = "Or provide the path explicitly via the `jre_path` argument."
      ))
    }
    cli::cli_inform("Using Java: {java_exe}")
  } else {
    # Legacy behavior: resolve from jre_path
    java_exe = file.path(jre_path, "bin", "java")
    if (.Platform$OS.type == "windows") {
      java_exe = paste0(java_exe, ".exe")
    }
    if (!file.exists(java_exe)) {
      if (file.exists(jre_path) && file.info(jre_path)$isdir == FALSE) {
        java_exe = jre_path
      } else {
        cli::cli_abort("Could not find java executable at {java_exe} or {jre_path}")
      }
    }
  }
  
  # Default threads

  if (is.null(nthreads)) {
    nthreads = get_default_beagle_threads()
  }
  
  if (!startsWith(memory, "-")) {
    memory = paste0("-", memory)
  }
  
  input_vcf = tempfile(pattern = "beagle_in_", fileext = ".vcf.gz")
  out_prefix = tempfile(pattern = "beagle_out_")
  output_vcf = paste0(out_prefix, ".vcf.gz")
  output_log = paste0(out_prefix, ".log")
  
  on.exit({
    if (file.exists(input_vcf)) unlink(input_vcf)
    if (file.exists(output_vcf)) unlink(output_vcf)
    if (file.exists(output_log)) unlink(output_log)
  }, add = TRUE)
  
  # Results are matched back by name, so names must exist and be unique
  loc_names = adegenet::locNames(gl)
  ind_names = adegenet::indNames(gl)
  if (is.null(loc_names) || is.null(ind_names) ||
      anyDuplicated(loc_names) > 0 || anyDuplicated(ind_names) > 0) {
    cli::cli_abort("Beagle imputation needs unique locus and individual names.")
  }
  chrom = adegenet::chromosome(gl)
  pos = adegenet::position(gl)
  if (is.null(chrom) || is.null(pos) || anyNA(chrom) || anyNA(pos)) {
    cli::cli_abort("Beagle imputation needs a chromosome and a position for every locus.")
  }
  
  # Beagle needs markers sorted by position within each chromosome
  locus_order = order(as.character(chrom), pos)
  
  # Convert the gl to a vcf using write_vcf (located in utils.R)
  write_vcf(gl[, locus_order], input_vcf, na_rep = "./.")
  
  # Pass the arguments one by one. Paths are quoted and the other values were
  # checked above, so nothing else can reach the command line.
  beagle_args = c(
    memory,
    "-jar", shQuote(beagle_path),
    paste0("gt=", shQuote(input_vcf)),
    paste0("out=", shQuote(out_prefix)),
    paste0("burnin=", format(burnin, scientific = FALSE)),
    paste0("iterations=", format(iterations, scientific = FALSE)),
    paste0("seed=", format(seed, scientific = FALSE)),
    paste0("nthreads=", format(nthreads, scientific = FALSE))
  )
  
  cli::cli_inform("Running Beagle on {length(loc_names)} loci and {length(ind_names)} individuals.")
  
  # Execute call
  status = system2(java_exe, beagle_args)
  if (status != 0) {
    cli::cli_abort("Beagle execution failed with exit status {status}")
  }
  
  if (!file.exists(output_vcf)) {
    cli::cli_abort("Expected Beagle output VCF file not found at {output_vcf}")
  }
  
  # Read back with cgiarGenomics read_vcf
  imputed_gl = read_vcf(output_vcf, ploidity = 2, na_reps = ".")
  imp_m = as.matrix(imputed_gl)
  
  # Put the output in the order of the input before using cell positions
  lost_loc = setdiff(loc_names, colnames(imp_m))
  lost_ind = setdiff(ind_names, rownames(imp_m))
  if (length(lost_loc) > 0 || length(lost_ind) > 0) {
    cli::cli_abort(c(
      "The Beagle output lacks {length(lost_loc)} loci and {length(lost_ind)} individuals of the input, so results cannot be matched.",
      "i" = "Loci that share a position, or whose alleles are not A, C, G or T, cause this. Remove them first."
    ))
  }
  imp_m = imp_m[ind_names, loc_names, drop = FALSE]
  
  orig_m = as.matrix(gl)
  called = !is.na(orig_m)
  if (any(imp_m[called] != orig_m[called], na.rm = TRUE)) {
    cli::cli_abort("The Beagle output changed calls that were not missing. Check the reference and alternative alleles.")
  }
  
  idx_na = which(!called)
  imp_dict = list()
  if (length(idx_na) > 0) {
    imputed_values = imp_m[idx_na]
    imp_dict = split(idx_na, imputed_values)
  }
  return(imp_dict)
}
