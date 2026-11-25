#' Flag markers that cannot enter a genlight object (internal)
#'
#' Markers without a position go to chromosome `unk` with consecutive
#' positions. Markers are flagged when they share a position with an earlier
#' marker, lack reference or alternative alleles, are multi-allelic or are
#' indels.
#'
#' @param geno_metadata Dataframe with the expected columns id, chrom, pos, ref, alt
#'
#' @return Same input dataframe but added a `filter` column if true the marker don't meet
#' the minimal requirements: no duplicated markers, physical position, bi-allelic,
#'  allelic information, SNP.
#'
#' @section Known issue:
#' Duplicated marker IDs are not renamed. See `dev/audit.md`, finding R1.
#' @noRd
process_metadata <- function(geno_metadata){
  # Expected fields id Chrom pos ref alt
  meta_columns <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  # Check if the columns in the input file match the expected column names
  if (length(intersect(meta_columns, colnames(geno_metadata))) != 5) {
    cli::cli_abort("metadata file doesn't have the expected column names,
                   have: {colnames(geno_metadata)}")
  }
  
  geno_metadata['filter'] <- FALSE
  
  # Flag markers with physical location
  geno_usable_idx <- which(!complete.cases(geno_metadata[,c('chrom','pos')]),)
  #geno_metadata[geno_usable_idx,'filter'] <- TRUE
  no_pos <- length(geno_usable_idx)
  
  if (no_pos > 0){
    cli::cli_inform("Were found {no_pos} markers
                   without physicall location, put into unk chrom and consecutive position")
    
    geno_metadata[geno_usable_idx, 'chrom'] <- "unk"
    geno_metadata[geno_usable_idx, 'pos'] <- seq(1, no_pos)
  }
  
  # Flag colocalized markers
  
  dup_pos_idx <- which(duplicated(geno_metadata[,c('chrom','pos')]))
  no_dup_pos <- length(dup_pos_idx)
  if (no_dup_pos > 0 ){
    cli::cli_inform("Were found {no_dup_pos} markers with duplicated position")
    geno_metadata[dup_pos_idx,'filter'] <- TRUE
  }
  
  # Check duplicated ids
  id_dup = duplicated(geno_metadata['id'])
  
  if(sum(id_dup) > 0){
    cli::cli_inform("Were found {sum(id_dup)} markers with duplicated id. \n
                  Renamed with chrom_pos nomenclature")  
  }
  
  geno_metadata <- geno_metadata %>%
    dplyr::mutate(
      id = dplyr::if_else(dplyr::row_number() %in% id_dup,  paste0(chrom, "_", pos), id),
      allele_count = stringr::str_count(alt, ",") + 1,
      ref_len = dplyr::if_else(stringr::str_detect(ref, "^[ACGT]$"),nchar(ref),NA),
      alt_len = dplyr::if_else(stringr::str_detect(alt, "^[ACGT]$"),nchar(stringr::str_replace(alt, ",", "")),NA) 
    ) 
  
  # Flag no reference alleles
  ref_alt_na_idx <- which(is.na(geno_metadata['ref_len']) | is.na(geno_metadata['alt_len']))
  no_ref_alt <- length(ref_alt_na_idx)
  
  if(no_ref_alt > 0){
    cli::cli_inform("Were found {no_ref_alt} without ref and alt allele data")
    geno_metadata[ref_alt_na_idx,'filter'] <- TRUE
  }
  
  
    # flag multiallelic 
  multi_allelic_idx <- which(geno_metadata['allele_count'] > 1)
  multi_allelic <- length(multi_allelic_idx)
  
  if(multi_allelic > 0){
    cli::cli_warn("Were found {multi_allelic} multi-allelic markers")
    geno_metadata[multi_allelic_idx,'filter'] <- TRUE
  }

  # flag indels
  
  indel_idx <- which(geno_metadata['ref_len'] > 1)
  indel_idx <- unique(c(indel_idx,which(geno_metadata['alt_len'] > 1)))
  
  
  indel <- length(indel_idx)
  
  if(indel > 0){
    cli::cli_inform("Were found {indel} indel markers")
    geno_metadata[indel_idx,'filter'] <- TRUE
  }
  
  return(geno_metadata)
}

#' From locus genotype call data, get the allelic dosage given the alleles and ploidity
#'
#' Generates every possible call for the two alleles with [get_all_gt_calls()],
#' converts each to the dosage of the alternative (second) allele with
#' [convert_gt_to_dosage()], and looks up the observed calls with
#' [replace_strings_with_integers()]. Calls that are not possible for these
#' alleles, including `NA`, become `NA`.
#'
#' @param l A vector of genotype calls, e.g., `c("AG", "GG", "AA")`.
#' @param alleles String. Reference and alternative allele separated by `/`,
#'   e.g., `"A/G"`.
#' @param ploidity Integer. The ploidity level of the organism.
#' @param sep String. Separator between alleles within a call (default `""`).
#'
#' @return A named integer vector with the dosage of the alternative allele
#'   for each call.
#' @export
#'
#' @examples
#' get_allelic_dosage(c("AG", "GG", "AA"), "A/G", 2)  # dosages 1, 2, 0
#' get_allelic_dosage(c("AG", "XX", "AA"), "A/G", 2)  # dosages 1, NA, 0
#' get_allelic_dosage(c("A/G", "G/G"), "A/G", 2, sep = "/")  # dosages 1, 2
get_allelic_dosage <- function(l, alleles, ploidity, sep = "") {
  alleles_c <- unlist(stringr::str_split(alleles, "/"))
  # All possible genotype calls
  possible_gt_calls <- get_all_gt_calls(alleles_c, ploidity, sep)
  # Get dosage given the alternative allele
  possible_dosage <- convert_gt_to_dosage(possible_gt_calls, alleles_c[2], ploidity,sep)
  dosages <- replace_strings_with_integers(possible_dosage, l)
  return(dosages)
}


#' Get all possible genotype calls given a unique set of alleles
#'
#' Generates all ordered combinations of `alleles` for the given ploidity, so
#' heterozygotes appear in every order (`AB` and `BA`).
#'
#' @param alleles Character vector of unique alleles, e.g., `c("A", "B")`.
#' @param ploidity Integer. The ploidity level of the organism.
#' @param sep String. Separator placed between alleles (default `""`).
#'
#' @return Character vector of all possible genotype calls.
#' @export
#'
#' @examples
#' get_all_gt_calls(c("A", "B"), 2)  # "AA" "AB" "BA" "BB"
#' get_all_gt_calls(c("A", "B"), 2, sep = "/")  # "A/A" "A/B" "B/A" "B/B"
#' length(get_all_gt_calls(c("A", "B", "C"), 3))  # 27
get_all_gt_calls <- function(alleles, ploidity, sep = "") {
  sep <- gsub("\\\\", "", sep)
  generate_calls <- function(prefix, ploidity, del = sep) {
    if (ploidity == 0) {
      if(nchar(sep) > 0){
        out <- substr(prefix, 1, nchar(prefix)-1)  
      } else {
        out <- substr(prefix, 1, nchar(prefix))  
      }
      
      return(out)
    }
    calls <- c()
    for (allele in alleles) {
      call <- paste0(prefix, allele,del)
      calls <- c(calls, generate_calls(call, ploidity - 1))
    }
    return(calls)
  }
  generate_calls("", ploidity)
}



#' Given a list of genotype calls, get the dosage of each one
#'
#' @param locus Character vector of genotype calls, e.g., `c("AG", "GG", "AA")`.
#' @param alt_allele String. The alternative allele, e.g., `"A"`.
#' @param ploidity Integer. The ploidity level of the organism, default is 2.
#' @param sep String. Separator between alleles within a call (default `""`).
#'
#' @return A named vector with the dosage of `alt_allele` in each call.
#' @export
#'
#' @examples
#' convert_gt_to_dosage(c("AG", "GG", "AA"), "A")  # 1, 0, 2
#' convert_gt_to_dosage(c("AAA", "GGG"), "A", 3)  # 3, 0
convert_gt_to_dosage <- function(locus, alt_allele, ploidity = 2,sep="") {
  sep <- gsub("\\\\", "", sep)
  l <- sapply(locus,
              genocall_to_allelic_dosage,
              alt_allele = alt_allele,
              ploidity = ploidity,
              sep=sep)
  return(l)
}

#' Genotype call to allelic dosage of alternative allele
#'
#' Counts the copies of `alt_allele` in one genotype call. Phased separators
#' (`|`) and `sep` are removed first. A missing call returns `NA`.
#'
#' @param genotype_call String. Genotype call, e.g., `"AG"`, `"AAA"`.
#' @param alt_allele String. Alternative allele, e.g., `"A"`.
#' @param ploidity Integer. Ploidity level of the organism, default is 2.
#' @param sep String. Separator between alleles within the call (default `""`).
#'
#' @return Integer. The dosage of `alt_allele`, or `NA`.
#' @export
#'
#' @examples
#' genocall_to_allelic_dosage("AG", "A")  # 1
#' genocall_to_allelic_dosage("GG", "A")  # 0
#' genocall_to_allelic_dosage("AAA", "A", 3)  # 3
#' genocall_to_allelic_dosage("A/G", "G", sep = "/")  # 1
#' genocall_to_allelic_dosage(NA, "A")  # NA
genocall_to_allelic_dosage <- function(genotype_call, alt_allele, ploidity = 2,sep="") {
  if (!is.na(nchar(genotype_call))) {
    # remove separators (and normalize phasing if present)
    genotype_call <- gsub("\\|", sep, genotype_call)
    if(sep != ""){
      genotype_call <- gsub(sep, "", genotype_call, fixed = TRUE)
    }
    # Genotype call successfully genotyped
    allele_length <- nchar(genotype_call) / ploidity
    
    # List with each allele as element
    split_genotype <- substring(genotype_call, seq(1, nchar(genotype_call), allele_length),
                                seq(allele_length, nchar(genotype_call), allele_length))
    
    # Matches of alt allele are the dosage
    dosage <- length(which(split_genotype == alt_allele))
    return(dosage)
  } else {
    # Genotype call missed
    return(NA)
  }
}


#' Replace a list of strings with their corresponding integer values
#'
#' This function takes two inputs: a named list or vector with string keys and integer values,
#' and a list of strings to be replaced. It replaces each string in the second list with the
#' corresponding integer value from the first list, based on the string keys.
#'
#' If a string in the second list does not have a corresponding key in the first list,
#' it will be replaced with NA.
#'
#' @param lookup_table A named list or vector with string keys and integer values.
#' @param strings_to_replace A list of strings to be replaced with their corresponding integer values.
#'
#' @return A list of integers, where each string in the input list has been replaced with its
#'         corresponding integer value from the lookup table, or NA if no match was found.
#' @export
#'
#' @examples
#' lookup <- c(A = 1, B = 2, C = 3)
#' strings <- c("B", "A", "D", "C")
#' replace_strings_with_integers(lookup, strings)  # Returns list(2, 1, NA, 3)
replace_strings_with_integers <- function(lookup_table, strings_to_replace) {
  # Use match to find the indices of the strings in the lookup table
  indices <- match(strings_to_replace, names(lookup_table))
  
  # Replace the strings with the corresponding integer values
  # or NA if no match was found
  integer_values <- lookup_table[indices]
  
  return(integer_values)
}


