
#' read_hapmap
#'
#' Reads a tab-separated hapmap file (uncompressed) and converts it into a
#' genlight object of alternative allele dosages. The first allele in the
#' `alleles` column is taken as reference and the second as alternative.
#' Markers that are not bi-allelic SNPs, or share a position with another
#' marker, are removed. Samples start at column 12.
#'
#' @param path Path to the hapmap file (uncompressed).
#' @param ploidity Integer, ploidity level of the organism (default: 2).
#' @param sep String, separator between alleles in a call, one of `""`
#'   (default, as in `AG`), `"/"` or `"|"`.
#'
#' @details
#' Diploid files that write each call as one IUPAC code (`A`, `R`, `Y`), as
#' TASSEL does, are decoded to two-letter genotypes.
#'
#' @return A genlight object with metrics from [recalc_metrics()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gl = read_hapmap("tests/hmp_fmt/diploid_hmp.txt", ploidity = 2)
#' gl4 = read_hapmap("tests/hmp_fmt/tetraploid_hmp.txt", ploidity = 4)
#' }
read_hapmap <- function(path, ploidity = 2, sep = "") {
  # Validate arguments
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  if (!rlang::is_integerish(ploidity)) {
    cli::cli_abort("`ploidity` must be a round number not {ploidity}")
  }
  
  if (!sep %in% c("", "/", "|")) {
    stop("sep must be '', '/', or '|'")
  }
  
  
  # Read the genotype data from the tabular file
  table <- read_tabular_geno(path)
  
  # Define the expected column names for the hapmap file
  hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos')
  
  # Check if the first 4 columns in the input file match the expected column names
  if (length(intersect(hapmap_snp_attr, colnames(table)[1:4])) != 4) {
    cli::cli_abort("Hapmap file doesn't have  first four standard column names: {hapmap_snp_attr}, have: {colnames(table)[1:4]}")
  }
  
  nrows <- dim(table)[1]
  ncols <- dim(table)[2]
  individuals <- colnames(table)[12:ncols]
  loci <- paste0(c(table[1:nrows, 3]), "_", c(table[1:nrows, 4]))
  
  # Print a data integrity message
  int_message <- paste0(
    "Input data should have loci as rows and individuals as columns \n",
    nrows, " Loci, confirming first 5: \n",
    paste(loci[1:5], collapse = " "), "\n",
    length(individuals), " Individuals, confirming first 5:\n",
    paste(individuals[1:5], collapse = " ")
  )
  
  cli::cli_inform(int_message)
  
  # Check if individual labels and loci IDs are unique
  if (length(unique(individuals)) != length(individuals)) {
    cli::cli_abort("Fatal Error: Individual labels are not unique, check and edit your input file\n")
  }
  
  # Metadata processing, get the first allele as reference and the remaining
  # as alternative
  meta <- table[1:11] %>% 
    dplyr::mutate(alt = purrr::map_chr(stringr::str_split(alleles, '/'), \(.x) paste(.x[-1], collapse = ",")),
           ref = stringr::str_split(alleles, '/', simplify = TRUE)[,1],
           ) %>% 
    dplyr::mutate(alt = dplyr::na_if(alt, "")) %>% 
    dplyr::select('rs#', chrom, pos, ref, alt) %>% 
    dplyr::rename(id = 'rs#')
  
  meta <- process_metadata(meta)
  
  
  # Read marker columns and transpose to obtain samples x snps
  mt <- t(table[!meta$filter, c(12:dim(table)[2]), drop = FALSE])
  
  # Count the alleles in each call, without the separator
  call_len <- nchar(if (sep == "") mt else gsub(sep, "", mt, fixed = TRUE))
  max_dosage <- max(call_len, na.rm = TRUE)
  
  # Diploid hapmap files from TASSEL write each call as one IUPAC code
  if (ploidity == 2 && max_dosage == 1) {
    iupac <- c(A = "AA", C = "CC", G = "GG", T = "TT",
               R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
    decoded <- unname(iupac[toupper(mt)])
    not_decoded <- sum(!is.na(mt) & is.na(decoded))
    cli::cli_inform("Calls are single IUPAC codes, read as diploid genotypes.")
    if (not_decoded > 0) {
      cli::cli_warn("{not_decoded} call{?s} with an unknown IUPAC code set to missing.")
    }
    mt[] <- decoded
    sep <- ""
    max_dosage <- 2
  }
  
  if(max_dosage < ploidity){
    cli::cli_warn("Max dosage ({max_dosage}) lower than ploidy lvl ({ploidity})")
  }
  
  if (max_dosage > ploidity){
    cli::cli_abort("Max dosage ({max_dosage}) higher than ploidy lvl ({ploidity})")
  }
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep='/')
  
  # One column per locus. vapply keeps the matrix shape with one sample or
  # one locus.
  gt <- vapply(seq_len(ncol(mt)), function(col) {
    as.numeric(get_allelic_dosage(mt[, col], allele_set[col], ploidity, sep))
  }, numeric(nrow(mt)))
  gt <- matrix(gt, nrow = nrow(mt))
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}

#' read_vcf
#'
#' Reads a VCF file (compressed or uncompressed) and converts it into a
#' genlight object of alternative allele dosages. Markers that are not
#' bi-allelic SNPs, or share a position with another marker, are removed.
#' Requires the suggested package vcfR.
#'
#' @param path String. Path to the VCF file. It could be compressed.
#' @param ploidity Integer. Ploidity level of the organism. (Default = 2)
#' @param na_reps Vector. Genotype strings to treat as missing
#'   (default `c("-", "./.")`). vcfR already converts fully missing calls.
#' @param sep Kept for compatibility and not used. Phased (`|`) and unphased
#'   (`/`) calls are both read, also when one file mixes them.
#'
#' @return A genlight object with metrics from [recalc_metrics()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gl = read_vcf("tests/vcf_fmt/diploid.vcf.gz", ploidity = 2)
#' gl
#' }
read_vcf <- function(path, ploidity = 2, na_reps = c("-", "./."), sep="/") {
  
  require_optional("vcfR")
  
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  if (!rlang::is_integerish(ploidity)) {
    cli::cli_abort("`ploidity` must be a round number not {ploidity}")
  }
  # Read the VCF file
  vcf <- vcfR::read.vcfR(path)
  
  # Get the metadata from the VCF file
  meta_vcf <- as.data.frame(vcfR::getFIX(vcf))
  
  meta <- meta_vcf %>% 
    dplyr::rename(id = ID,
           chrom = CHROM,
           pos = POS,
           ref = REF,
           alt = ALT) %>% 
  dplyr::select(id, chrom, pos, ref, alt)
  
  meta <- process_metadata(meta)
  mt <- t(vcfR::extract.gt(vcf, return.alleles = TRUE)[!meta$filter, , drop = FALSE])
  
  # A file can mix phased (|) and unphased (/) calls. Write all of them with
  # one separator, so both are read whatever `sep` says.
  mt_dim <- dim(mt)
  mt_names <- dimnames(mt)
  mt <- gsub("|", "/", mt, fixed = TRUE)
  dim(mt) <- mt_dim
  dimnames(mt) <- mt_names
  sep <- "/"
  
  na_reps <- unique(c(na_reps, gsub("|", "/", na_reps, fixed = TRUE)))
  mt[mt %in% na_reps] <- NA
  # calls made of dots only are missing at any ploidy (./. or ./././.)
  mt[grepl("^[./]+$", mt)] <- NA
  
  # check ploidity
  call_len <- nchar(gsub("/", "", mt, fixed = TRUE))
  max_dosage <- max(call_len, na.rm = TRUE)
  if(max_dosage < ploidity){
    cli::cli_warn("Max dosage ({max_dosage}) lower than ploidy lvl ({ploidity})")
  }
  
  if (max_dosage > ploidity){
    cli::cli_abort("Max dosage ({max_dosage}) higher than ploidy lvl ({ploidity})")
  }
  
  individuals <- rownames(mt)
  
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep="/")

  gt <- vapply(seq_len(ncol(mt)), function(col) {
    as.numeric(get_allelic_dosage(mt[, col], allele_set[col], ploidity, sep))
  }, numeric(nrow(mt)))
  gt <- matrix(gt, nrow = nrow(mt))
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}


#' read_DArTSeq_SNP
#'
#' Reads a standard DArTSeq SNP report (CSV) with
#' `dartR.base::gl.read.dart()` and moves chromosome and position from the
#' locus metadata into the genlight slots. Unmapped SNPs are removed. Diploid
#' only. Requires the suggested package dartR.base.
#'
#' @param path Path to the DArTSeq CSV file.
#' @param snp_id Column name of the SNP ID data in the CSV file.
#' @param chr_name Column name of the chromosome where the SNP is located in the CSV file.
#' @param pos_name Column name of the physical position where the SNP is located in the CSV file.
#'
#' @return A genlight object with mapped SNPs and metrics from [recalc_metrics()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' read_DArTSeq_SNP("tests/DartSeq_fmt/DartSeq/DArTSeq_CB_snp.csv",
#'                  snp_id = "AlleleID", chr_name = "Chrom", pos_name = "ChromPos")
#' }
read_DArTSeq_SNP <- function(path, snp_id, chr_name, pos_name) {
  
  require_optional("dartR.base")
  
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  
  # Read the DArTSeq CSV file
  gl <- dartR.base::gl.read.dart(path)
  
  if(!hasName(gl@other$loc.metrics, snp_id)){
    cli::cli_abort("The input snp_id: {snp_id} column doesn't exist in the dartSeq file")
  }
  
  if(!hasName(gl@other$loc.metrics, chr_name)){
    cli::cli_abort("The input chr_name: {chr_name} column doesn't exist in the dartSeq file")
  }
  
  if(!hasName(gl@other$loc.metrics, pos_name)){
    cli::cli_abort("The input pos_name: {pos_name} column doesn't exist in the dartSeq file")
  }
  
  # DArTSeq stores position and chromosome information in the 'other' slot
  # It is necessary to modify it to be compatible with the native genlight object
  pos <- c(gl@other$loc.metrics[, pos_name])
  chrom <- c(gl@other$loc.metrics[, chr_name])
  
  # Remove SNPs where the tag is unmapped (no position data)
  no_pos <- which(is.na(pos) | pos <= 0)
  no_chrom <- which(is.na(chrom))
  no_loc <- unique(c(no_pos, no_chrom))
  
  if(length(no_loc) > 0){
    cli::cli_inform("Were detected {length(no_loc)} \n
                    SNPs withoud position data, removed.")
  }
  
  # Remove SNPs without position data. Skip when every SNP is mapped, since
  # indexing with an empty negative vector would drop every locus.
  if(length(no_loc) > 0){
    gl <- gl[, -c(no_loc)]
    pos <- pos[-c(no_loc)]
    chrom <- chrom[-c(no_loc)]
  }
  
  # Assign the position and chromosome of the loci that were kept
  adegenet::position(gl) <- pos
  adegenet::chromosome(gl) <- chrom
  
  gl <- recalc_metrics(gl)
  return(gl)
  
}

#' extract_dart_gt
#'
#' Extracts the genotype block from a DArT report read without a header, and
#' names its columns by sample.
#'
#' @param data Dataframe. DArT report read with `header = FALSE`.
#' @param sample.name.row Integer. Row number which the sample names are (default = 2)
#' @param first.data.row Integer. Row number of the first genotype call (default = 9)
#' @param first.data.col Integer. Col number of the first genotype call (default = 6)
#'
#' @return Dataframe. The genotype block, with the sample names as column
#'   names. Values are kept as read.
#' @keywords internal
extract_dart_gt <- function(data,
                            sample.name.row= 2,
                            first.data.row = 9,
                            first.data.col=6){
  
  if (!is.data.frame(data)) {
    cli::cli_abort("`data` must be a data frame.")
  }
  if (!rlang::is_integerish(sample.name.row)) {
    cli::cli_abort("`sample.name.row` must be a round number not {sample.name.row}")
  }
  if (!rlang::is_integerish(first.data.row)) {
    cli::cli_abort("`first.data.row` must be a round number not {first.data.row}")
  }
  if (!rlang::is_integerish(first.data.col)) {
    cli::cli_abort("`first.data.col` must be a round number not {first.data.col}")
  }
  
  # rows with data
  rows <- first.data.row:nrow(data)
  cols <- first.data.col:ncol(data)
  # sample names
  id <- as.character(data[sample.name.row,cols])
  # Samples with duplicated id?
  dupes <- unique(id[which(duplicated(id))])
  if(length(dupes)>0){
    cli::cli_abort("The samples: {dupes} are duplicated, fix it.")
  }
  out_data <- data[rows, cols]
  colnames(out_data) <- id
  return(out_data)
}

#' read_one_row_dart_gt_data
#'
#' Reads a DArT report where each variant takes a single row, such as a DArTag
#' dosage report, and splits it into genotype calls and marker metadata.
#'
#' @param path Path to the DArT CSV file.
#' @param sample.name.row Integer. Row holding the sample names (default 2).
#' @param first.data.row Integer. First row with genotype data (default 9).
#'   The row above it must hold the metadata column names.
#' @param first.data.col Integer. First column with genotype data (default 6).
#' @param snp_id_col Character. Name of the metadata column with marker IDs.
#'
#' @return A list with `gt` (a data frame of integer dosages, samples in
#'   columns, plus a `marker_id` column) and `dosage_meta` (the metadata
#'   columns).
#' @keywords internal
read_one_row_dart_gt_data <- function(path,
                                      sample.name.row= 2,
                                      first.data.row = 9,
                                      first.data.col=6,
                                      snp_id_col){
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  
  data_raw <- read.csv(path, header=F, check.names=F)
  data <- extract_dart_gt(data_raw,
                          sample.name.row,
                          first.data.row,
                          first.data.col)
  
  
  # Markers with duplicated id
  meta_data <- data_raw[first.data.row:nrow(data_raw), 1:first.data.col-1]
  colnames(meta_data) <- data_raw[first.data.row-1, 1:first.data.col-1]
  
  if(!hasName(meta_data, snp_id_col)){
    cli::cli_abort("The input `snp_id_col`: {snp_id_col} column doesn't exist in the file")
  }
  
  marker_ids <- as.vector(unlist(meta_data[snp_id_col]))
  l_dupes <- unique(marker_ids[which(duplicated(marker_ids))])
  
  if(length(l_dupes) > 0){
    cli::cli_warn("{length(l_dupes)} marker ID{?s} appear{?s/} more than once in the dosage report. The first row of each is kept.")
    keep <- !duplicated(marker_ids)
    data <- data[keep, , drop = FALSE]
    meta_data <- meta_data[keep, , drop = FALSE]
    marker_ids <- marker_ids[keep]
  }
  
  geno <- apply(data,2,as.integer)
  
  geno <- as.data.frame(geno)
  geno$marker_id <- marker_ids
  
  return(list(gt = geno, dosage_meta = meta_data))
}

#' Parse reference and alternative alleles from a DArTag counts file (internal)
#'
#' Pivots the allele sequences of each marker to wide format and finds the
#' single base where reference and alternative tags differ with
#' `compare_tags()`.
#'
#' @return A data frame with one row per marker, the tag sequences and the
#'   columns `ref_allele` and `alt_allele`.
#' @noRd
parse_counts_metadata <- function(path,
                                  sample.name.row= 7,
                                  first.data.row = 9,
                                  first.data.col=6,
                                  snp_id_col = 'CloneID',
                                  allele_id_col = 'AlleleID',
                                  allele_sep = "\\|",
                                  allele_sequence_col = 'AlleleSequence',
                                  ref_name = 'Ref',
                                  alt_name = 'Alt'){
  
  data_raw <- read.csv(path, header=F, check.names=F)
  
  meta_data <- data_raw[first.data.row:nrow(data_raw), 1:first.data.col-1]
  colnames(meta_data) <- data_raw[first.data.row-1, 1:first.data.col-1]
  
  if(!hasName(meta_data, snp_id_col)){
    cli::cli_abort("The input `snp_id_col`: {snp_id_col} column doesn't exist in the file")
  }
  if(!hasName(meta_data, allele_id_col)){
    cli::cli_abort("The input `allele_id_col`: {allele_id_col} column doesn't exist in the file")
  }
  if(!hasName(meta_data, allele_sequence_col)){
    cli::cli_abort("The input `allele_sequence_col`: {allele_sequence_col} column doesn't exist in the file")
  }
  
  meta_data <- meta_data %>% 
    dplyr::mutate(allele_type = stringr::str_split_fixed(!!rlang::sym(allele_id_col), allele_sep, 2)[, 2])
  
  # Wider the metadata table
  meta_data_wide <- meta_data %>% 
    dplyr::select(!!rlang::sym(snp_id_col),
           !!rlang::sym(allele_sequence_col),
           allele_type) %>% 
    tidyr::pivot_wider(
      id_cols = !!rlang::sym(snp_id_col),
      names_from = allele_type,
      values_from = !!rlang::sym(allele_sequence_col)
    )
  
  alleles <- purrr::map2_dfr(as.vector(unlist(meta_data_wide[,ref_name])),
                             as.vector(unlist(meta_data_wide[,alt_name])),
                             \(.x,.y) compare_tags(.x, .y))
  
  wide_alleles <- cbind(meta_data_wide, alleles)
  
  
  return(wide_alleles)
}

#' Read a DArTag dosage report with its allele counts file
#'
#' Builds a genlight object from the two files DArT delivers for DArTag
#' panels, which also work for polyploids. Dosages come from the dosage
#' report. Reference and alternative alleles come from the tag sequences in
#' the counts file, where the single differing base defines the SNP. Markers
#' without allele sequences are removed with a warning.
#'
#' @param dosage_path Path to the dosage report (CSV).
#' @param counts_path Path to the allele match counts file (CSV), collapsed.
#' @param ploidity Integer. Ploidy of the organism.
#' @param sample.name.row Integer. Row with sample names (default 7).
#' @param first.data.row Integer. First data row (default 9).
#' @param first.data.col Integer. First data column (default 6).
#' @param snp_id_col_dosage Character. Marker ID column in the dosage report
#'   (default `"MarkerID"`).
#' @param snp_id_col_count Character. Marker ID column in the counts file
#'   (default `"CloneID"`).
#' @param allele_id_col Character. Allele ID column in the counts file
#'   (default `"AlleleID"`).
#' @param allele_sep Regular expression separating marker and allele type in
#'   `allele_id_col` (default `"\\|"`).
#' @param allele_sequence_col Character. Tag sequence column (default
#'   `"AlleleSequence"`).
#' @param ref_name,alt_name Character. Allele type labels for reference and
#'   alternative tags (default `"Ref"` and `"Alt"`).
#' @param dosage_col_chr,dosage_col_pos Character. Chromosome and position
#'   columns in the dosage report (default `"Chrom"` and `"ChromPos"`).
#'
#' @return A genlight object with metrics from [recalc_metrics()].
#' @export
#'
#' @examples
#' \dontrun{
#' gl = read_DArTag_count_dosage(
#'   dosage_path = "tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_Dose_Report.csv",
#'   counts_path = "tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_match_counts_collapsed.csv",
#'   ploidity = 6)
#' }
read_DArTag_count_dosage <- function(dosage_path,
                                     counts_path,
                                     ploidity,
                                     sample.name.row = 7,
                                     first.data.row = 9,
                                     first.data.col = 6,
                                     snp_id_col_dosage = 'MarkerID',
                                     snp_id_col_count = 'CloneID',
                                     allele_id_col = 'AlleleID',
                                     allele_sep = "\\|",
                                     allele_sequence_col = 'AlleleSequence',
                                     ref_name = 'Ref',
                                     alt_name = 'Alt',
                                     dosage_col_chr = 'Chrom',
                                     dosage_col_pos = 'ChromPos'){
  
  if (!file.exists(dosage_path)){
    cli::cli_abort("`dosage_path` don't exist. Verify if is writed properly {dosage_path}")
  }
  if (!file.exists(counts_path)){
    cli::cli_abort("`counts_path` don't exist. Verify if is writed properly {counts_path}")
  }
  # Parse metadata (Allele sequences)
  meta_data <- parse_counts_metadata(
    counts_path,
    sample.name.row,
    first.data.row,
    first.data.col,
    snp_id_col_count,
    allele_id_col,
    allele_sep,
    allele_sequence_col,
    ref_name,
    alt_name)
  
  # Parse genotype call information
  dosage_data <- read_one_row_dart_gt_data(
    dosage_path, 
    sample.name.row,
    first.data.row,
    first.data.col,
    snp_id_col_dosage)
  
  gt <- dosage_data$gt
  
  dosage_meta <- dosage_data$dosage_meta
  
  if(!hasName(dosage_meta, dosage_col_chr)){
    cli::cli_abort("The input `dosage_col_chr`: {dosage_col_chr} column doesn't exist in the file")
  }
  if(!hasName(dosage_meta, dosage_col_pos)){
    cli::cli_abort("The input `dosage_col_pos`: {dosage_col_pos} column doesn't exist in the file")
  }
  
  
  dosage_meta[,dosage_col_pos] <- as.integer(dosage_meta[,dosage_col_pos])
  
  # Check if markers reported in dosage are in counts 
  
  orphan_snps <- which(!gt$marker_id %in% meta_data[, snp_id_col_count])
  if(length(orphan_snps) > 0){
    cli::cli_warn("The markers: {gt$marker_id[orphan_snps]} without Allele sequence data, removed.")
  }
  
  
  merged_data <- merge(gt, meta_data,
                       by.x = "marker_id",
                       by.y = snp_id_col_count)
  
  
  
  merged_data <- merge(merged_data, dosage_meta,
                       by.x = "marker_id",
                       by.y = snp_id_col_dosage)
  
  
  nrows <- dim(gt)[1]
  ncols <- dim(gt)[2]
  individuals <- colnames(gt)[which(colnames(gt) != "marker_id")]
  
  
  loci <- paste0(c(merged_data[1:5,dosage_col_chr]), "_", c(merged_data[1:5,dosage_col_pos]))
  
  report_bullets <- cli::cli_bullets(c("Input data should have loci as rows and individuals as columns.",
                                       "i" = "{nrows} Loci, first 5: {loci}",
                                       "i" = "{length(individuals)} Individuals, first 5: {individuals[1:5]}"))
  
  cli::cli_inform(report_bullets)
  
  # Metadata processing, get the first allele as reference and the remaining
  # as alternative
  meta <- merged_data[,c('marker_id', dosage_col_chr, dosage_col_pos, 'ref_allele', 'alt_allele')]
  colnames(meta) <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  
  meta <- process_metadata(meta)
  
  # Read marker columns and transpose to obtain samples x snps
  mt <- t(merged_data[!meta$filter, individuals])
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep='/')
  
  max_dosage <- max(mt, na.rm = TRUE)
  if(max_dosage < ploidity){
    cli::cli_warn("Max dosage ({max_dosage}) lower than ploidy lvl ({ploidity})")
  } 
  if (max_dosage > ploidity){
    cli::cli_abort("Max dosage ({max_dosage}) higher than ploidy lvl ({ploidity})")
  }
  
  gl <- new("genlight",
            mt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}

#' Find the single SNP between reference and alternative tags (internal)
#'
#' @param ref,alt Character. Tag sequences of equal length.
#' @return A named character vector with `ref_allele` and `alt_allele`, or `N`
#'   for both when the tags differ at zero or several bases.
#' @noRd
compare_tags <- function(ref, alt){
  bases <- c("A", "C", "G", "T")
  hap.ref <- unlist(strsplit(ref, split = ""))
  hap.alt <- unlist(strsplit(alt, split = ""))
  k <- which(hap.ref != hap.alt & hap.ref %in% bases & 
               hap.alt %in% bases)
  if (length(k) == 1) {
    return(c(ref_allele = hap.ref[k], alt_allele = hap.alt[k]))
  }
  else {
    return(c(ref_allele = "N", alt_allele = "N"))
  }
}


