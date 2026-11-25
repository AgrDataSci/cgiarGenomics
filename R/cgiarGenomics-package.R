#' cgiarGenomics: reading, quality control and imputation of SNP marker data
#'
#' Reads bi-allelic SNP marker data from hapmap, VCF, DArTSeq and DArTag files
#' for diploid and polyploid crops, and converts them to allele dosages in an
#' adegenet genlight object. Provides quality metrics, sequential filtering,
#' imputation, merging of replicate samples, identity and Mendelian checks,
#' and export to VCF. It is the genomic data engine of Bioflow, the CGIAR
#' breeding analytics pipeline.
#'
#' Start with [read_hapmap()], [read_vcf()], [read_DArTSeq_SNP()] or
#' [read_DArTag_count_dosage()], then [apply_sequence_filtering()],
#' [impute_gl()] and [write_vcf()].
#'
#' @keywords internal
"_PACKAGE"

#' @import adegenet
#' @importFrom methods new
#' @importFrom stats complete.cases na.omit predict setNames
#' @importFrom utils hasName read.csv write.table
NULL
