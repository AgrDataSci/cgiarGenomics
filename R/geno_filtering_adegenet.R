#' Filter, then sample loci in proportion to chromosome size
#'
#' Filters individuals and loci by missingness and MAF with
#' [apply_sequence_filtering()], then draws `size` loci at random, allocated to
#' each chromosome in proportion to its number of loci. Loci left over by
#' rounding go to the chromosome with the most loci. Useful to thin a panel
#' before [ibs_matrix_purrr()] or other costly steps.
#'
#' @param gl A genlight object with metrics from [recalc_metrics()].
#' @param ind_miss Numeric. Keep individuals with missingness `<=` this value.
#' @param loc_miss Numeric. Keep loci with missingness `<=` this value.
#' @param maf Numeric. Keep loci with MAF `>` this value.
#' @param size Integer. Number of loci to return.
#' @param seed Integer or `NULL` (default). Seed for the random draw.
#'
#' @return A genlight object with at most `size` loci. Metrics in `other` are
#'   not recalculated, so run [recalc_metrics()] before filtering it again.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gl_small = random_select_loci(gl, ind_miss = 0.2, loc_miss = 0.2,
#'                               maf = 0.05, size = 1000, seed = 1)
#' }
random_select_loci <- function(gl, ind_miss, loc_miss, maf, size, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  filter_1 <- list("ind_miss", "<=", ind_miss)
  filter_2 <- list("loc_miss", "<=", loc_miss)
  filter_3 <- list("maf", ">", maf)
  
  filtering_steps <- list(
    filter_1,
    filter_2,
    filter_3)
  
  filt_seq <- lapply(filtering_steps, function(x){
    setNames(as.list(x), c("param", "operator", "threshold"))
  })
  
  filt_gl <- apply_sequence_filtering(gl, filt_seq)
  
  sites <- data.frame(loc_id = adegenet::locNames(filt_gl$gl),
                      CHROM = filt_gl$gl@chromosome)
  
  group_sizes <- sites %>% 
    dplyr::group_by(CHROM) %>% 
    dplyr::summarize(total = dplyr::n()) %>% 
    dplyr::mutate(q =floor((size*total)/dim(sites)[1]))
  
  remaining <- size - sum(group_sizes$q)
  
  max_total <- group_sizes %>% 
    dplyr::filter(total == max(group_sizes$total)) %>% 
    dplyr::pull(CHROM)
  
  out <- sites %>%
    dplyr::group_by(CHROM) %>%
    dplyr::group_modify(function(sdf, i_chr) {
      k <- group_sizes %>% dplyr::filter(CHROM == i_chr$CHROM) %>%  dplyr::pull(q)
      # with a tie, the first of the largest chromosomes takes the remainder
      if(i_chr$CHROM == max_total[1]){
        k = k + remaining
      }
      dplyr::slice_sample(sdf, n = k, replace = FALSE)
    }) %>%
    dplyr::ungroup()
  
    return(filt_gl$gl[,as.vector(out$loc_id)])
  }
  
#' Evaluate one filter on a genlight object
#'
#' Finds the individuals or loci whose metric passes a comparison. The metric
#' is looked up by name in `other$loc.metrics` or `other$ind.metrics`, which
#' decides whether the filter acts on loci or on individuals. The object is not
#' changed. [apply_sequence_filtering()] applies the result.
#'
#' @param gl A genlight object with metrics from [recalc_metrics()].
#' @param parameter Character. Metric name, one of `maf`, `loc_miss`,
#'   `loc_het`, `loc_Fis`, `ind_miss` or `ind_het`.
#' @param threshold Numeric, at most 1.
#' @param comparison_operator Character. One of `">"`, `">="`, `"<"` or `"<="`.
#'   Elements for which `metric <operator> threshold` is `TRUE` are kept.
#'
#' @return A list with `param`, `operator`, `threshold`, `filter_margin`
#'   (`"loc"` or `"ind"`), `index` (positions of the elements kept) and
#'   `filter_out` (names of the elements removed).
#' @seealso [apply_sequence_filtering()]
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = recalc_metrics(new("genlight", mt, ploidy = 2))
#' filter_gl(gl, parameter = "loc_miss", threshold = 0.2,
#'           comparison_operator = "<=")
filter_gl <- function(gl, parameter, threshold, comparison_operator){
  
  # Verify if parameter exist on the gl and get the margin (ind, loc)
  filter_margin <- get_parameter_margin(gl, parameter)
  comparison_operator <- match.arg(comparison_operator, choices = c(">", ">=", "<", "<="))
  comparison_func <- match.fun(comparison_operator)
  # Verify if threshold is a float value
  if(threshold > 1 ){
    cli::cli_abort("`threshold`: {threshold} is greater than 1. Metrics are rates between 0 and 1.")
  }
  
  if(filter_margin == 'loc'){
    index <- which(comparison_func(gl@other$loc.metrics[parameter], threshold))
    filter_out <- gl@loc.names[setdiff(seq_along(gl@loc.names), index)]
    
  } else {
    index <- which(comparison_func(gl@other$ind.metrics[parameter], threshold))
    filter_out <- gl@ind.names[setdiff(seq_along(gl@ind.names), index)]
    
  }
  
  out <- list(param = parameter,
              operator = comparison_operator,
              threshold = threshold,
              filter_margin = filter_margin,
              index = index,
              filter_out = filter_out)
  
  return(out)
}

#' Tell whether a metric applies to loci or individuals
#'
#' @param gl A genlight object with metrics from [recalc_metrics()].
#' @param param_name Character. Metric name, as in [filter_gl()]. Partial
#'   names are matched with [match.arg()].
#'
#' @return `"loc"` for a locus metric, `"ind"` for an individual metric.
#'   Errors when the name matches no metric.
#' @export
get_parameter_margin <- function(gl, param_name){
  # Get the expected parameters
  loc_metric_names <- colnames(gl@other$loc.metrics)
  ind_metric_names <- colnames(gl@other$ind.metrics)
  
  param_name = match.arg(param_name, choices = c(loc_metric_names, ind_metric_names))
  
  if(param_name %in% loc_metric_names){
    by = 'loc'
  } else {
    by = 'ind'
  }
  return(by)
}

#' Apply a sequence of filters to a genlight object
#'
#' Runs filters in the given order, recalculating metrics after each one, so
#' every step sees the object left by the previous step. Two steps are always
#' added at the end to remove loci and individuals with only missing data.
#' When no locus or no individual passes a step, the function stops with an
#' error of class `cgiarGenomics_empty_filter`, because the step would remove
#' everything.
#'
#' @param gl A genlight object with metrics from [recalc_metrics()].
#' @param filt_sequence A list of filter steps. Each step is a list with
#'   elements `param`, `operator` and `threshold`, as in [filter_gl()].
#'
#' @return A list with `gl` (the filtered genlight object) and `filt_log` (one
#'   [filter_gl()] result per step, named `<param>_<step>`).
#'
#' @seealso [filter_gl()], [get_parameter_margin()]
#' @export
#'
#' @examples
#' mt = matrix(c(0, 1, 2, NA,
#'               2, 2, 1, 0,
#'               0, 1, NA, 2),
#'             nrow = 4,
#'             dimnames = list(paste0("ind", 1:4), paste0("snp", 1:3)))
#' gl = recalc_metrics(new("genlight", mt, ploidy = 2))
#' steps = list(
#'   list(param = "loc_miss", operator = "<=", threshold = 0.3),
#'   list(param = "ind_miss", operator = "<=", threshold = 0.5)
#' )
#' out = apply_sequence_filtering(gl, steps)
#' out$gl
#' names(out$filt_log)
apply_sequence_filtering <- function(gl, filt_sequence){
  if(!rlang::is_bare_list(filt_sequence)){
    cli::cli_abort("Provide a list of filter operations in `filt_sequence`")
  }
  if(!inherits(gl, "genlight")){
    cli::cli_abort("`gl` is not a genlight class")
  }
  
  # Allways add at the end a filter step to remove all NA loci and ind
  
  locNA_filt <- list("loc_miss", "<", 1)
  indNA_filt <- list("ind_miss", "<", 1)
  
  
  allNA_steps <- list(
    locNA_filt,
    indNA_filt
  )
  
  
  filt_NA_seq <- lapply(allNA_steps, function(x){
    setNames(as.list(x), c("param", "operator", "threshold"))
  })
  filt_sequence <- append(filt_sequence, filt_NA_seq)
  # Duplicate the gl object to perform the filtering
  working_gl <- gl
  filtering_log <- list()
  previous_margin <- ""
  
  for (i_step in 1:length(filt_sequence)){
    filt_step <- filt_sequence[[i_step]]
    param <- filt_step[['param']]
    threshold <- filt_step[['threshold']]
    operator <- filt_step[['operator']]
  
    i_filt_out <- filter_gl(working_gl,
                            parameter = param,
                            threshold = threshold,
                            comparison_operator = operator)
    
    i_margin <- get_parameter_margin(gl, param)
    
    # A step that nothing passes would empty the object. Stop and say so,
    # instead of skipping the step and returning unfiltered data.
    if(length(i_filt_out$index) == 0){
      what <- if(i_margin == "loc") "locus" else "individual"
      cli::cli_abort(
        c("No {what} passes filter step {i_step} ({param} {operator} {threshold}).",
          "i" = "The step would remove everything, so filtering stopped. Relax the threshold or drop the step."),
        class = "cgiarGenomics_empty_filter"
      )
    }
    if(i_margin == "loc"){
      working_gl <- working_gl[,i_filt_out$index]
    } else {
      working_gl <- working_gl[i_filt_out$index,]
    }
    working_gl <- recalc_metrics(working_gl)
    filtering_log[[glue::glue("{param}_{i_step}")]] <- i_filt_out
  }
  
  return(list(gl = working_gl, filt_log = filtering_log))
}

