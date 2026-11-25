# Function reference

Every exported function of cgiarGenomics 1.0.2, grouped by the step of the workflow it serves, plus the documented helpers marked *internal*. Call internal helpers with `cgiarGenomics:::`, since they may change without notice. The help pages in R (`?read_vcf`) hold the same content in more detail. Functions marked ⚠ have a known issue described in [`dev/audit.md`](dev/audit.md), with the finding number in brackets.

All readers return a [genlight](https://github.com/thibautjombart/adegenet) object of alternative allele dosages, with individuals in rows and loci in columns. Each one stores quality metrics in the `other` slot, which the filtering functions read.

| Slot | Columns | One row per |
|---|---|---|
| `gl@other$loc.metrics` | `maf`, `loc_miss`, `loc_het`, `loc_Fis` | locus |
| `gl@other$ind.metrics` | `ind_miss`, `ind_het` | individual |

## Read genotype data

| Function | What it does | Ploidy |
|---|---|---|
| `read_hapmap(path, ploidity = 2, sep = "")` | Reads a tab-separated hapmap file. The first allele in `alleles` is the reference. Samples start at column 12. ⚠ (R4, R6) | any |
| `read_vcf(path, ploidity = 2, na_reps = c("-", "./."), sep = "/")` | Reads a VCF, compressed or not, through vcfR. Use `sep = "\\|"` for phased files. ⚠ (R5) | any |
| `read_DArTSeq_SNP(path, snp_id, chr_name, pos_name)` | Reads a standard DArTSeq SNP report through dartR.base and keeps mapped SNPs. ⚠ (R2, R3) | diploid |
| `read_DArTag_count_dosage(dosage_path, counts_path, ploidity, ...)` | Reads a DArTag dosage report with its collapsed allele counts file. Alleles come from the single base that differs between tags. | any |
| `read_one_row_dart_gt_data(path, sample.name.row, first.data.row, first.data.col, snp_id_col)` *internal* | Splits a one-row-per-variant DArT report into dosages and metadata. Used by `read_DArTag_count_dosage()`. | any |
| `extract_dart_gt(data, sample.name.row = 2, first.data.row = 9, first.data.col = 6)` *internal* | Extracts the genotype block of a DArT report read without header. | any |

Every reader drops markers that are not bi-allelic SNPs, that lack alleles, or that share a position with an earlier marker. Markers without a position are kept on chromosome `unk`. A message reports each case.

## Inspect quality

| Function | Returns |
|---|---|
| `recalc_metrics(gl)` | The object with fresh `loc.metrics` and `ind.metrics`. Run it after any change made outside cgiarGenomics. |
| `get_overall_summary(gl)` | List with `nind`, `nloc`, `ov_miss`, `ov_het`, `ov_maf`, from stored metrics. |
| `get_loc_missing(gl)` | Missing rate per locus. |
| `get_ind_missing(gl)` | Missing rate per individual. |
| `get_overall_missingness(gl)` | Share of missing calls in the whole matrix. |
| `get_maf(gl, ploidy = 2)` | Minor allele frequency per locus. ⚠ (M1, values are divided by ploidy twice) |
| `get_loc_heterozygosity(gl, ploidy = 2)` | Share of heterozygous calls per locus. ⚠ (M2) |
| `get_ind_heterozygosity(gl, ploidy = 2)` | Share of heterozygous calls per individual. ⚠ (M2) |
| `get_heterozygosity_metrics(gl, ploidy = 2)` | Both of the above in one list (`het_ind`, `het_loc`). ⚠ (M2) |

## Filter

| Function | What it does |
|---|---|
| `apply_sequence_filtering(gl, filt_sequence)` | Applies a list of filter steps in order, recalculating metrics after each one, then removes loci and individuals with no data. Returns `gl` and `filt_log`. ⚠ (F1, a step that nothing passes is skipped) |
| `filter_gl(gl, parameter, threshold, comparison_operator)` | Evaluates one step without changing the object. Returns the kept positions and removed names. |
| `get_parameter_margin(gl, param_name)` | Says whether a metric acts on loci (`"loc"`) or individuals (`"ind"`). |
| `random_select_loci(gl, ind_miss, loc_miss, maf, size, seed)` | Filters, then samples `size` loci in proportion to chromosome size. ⚠ (F3) |

A filter step is a list with `param`, `operator` and `threshold`.

```r
steps = list(
  list(param = "loc_miss", operator = "<=", threshold = 0.2),
  list(param = "ind_miss", operator = "<=", threshold = 0.3),
  list(param = "maf",      operator = ">",  threshold = 0.05)
)
out = apply_sequence_filtering(gl, steps)
```

Valid `param` values are `maf`, `loc_miss`, `loc_het`, `loc_Fis`, `ind_miss` and `ind_het`. Operators are `>`, `>=`, `<` and `<=`, and thresholds cannot exceed 1.

## Impute

| Function | What it does |
|---|---|
| `impute_gl(gl, ploidity = 2, method, ...)` | Imputes missing dosages with `"frequency"`, `"random_forest"` or `"beagle"`. Returns `gl` and `log`, the imputed cells grouped by dosage. ⚠ (I1 for random forest, I6) |
| `impute_beagle(gl, jre_path = NULL, beagle_path = NULL, memory = "Xmx4g", burnin = 3, iterations = 12, seed = -99999, nthreads = NULL)` | Runs Beagle 5 on a temporary VCF. Diploids only. Usually called through `impute_gl(method = "beagle")`. ⚠ (I4, I5) |
| `i_freq_impute(q_frq, ploidity = 2)` *internal* | Draws one dosage from an allele frequency. |
| `set_fraction_na(gl, fraction = 0.1)` | Masks a random share of calls to test imputation. ⚠ (I3) |
| `get_accuracy(ref_gl, imp_gl)` | Share of reference calls reproduced, per locus. ⚠ (I2, inflated) |

Method-specific arguments go through `...` by name. Random forest takes `nflank`, `ntree`, `seed` and `nthreads`. Beagle takes `jre_path`, `beagle_path`, `memory`, `burnin`, `iterations`, `seed` and `nthreads`. Unknown names are an error.

## Identity, replicates and pedigree

| Function | What it does |
|---|---|
| `ibs_matrix_purrr(gl, ploidy)` | Pairwise identity by state and the number of loci shared, for purity and identity checks. Slow on many individuals. |
| `merge_duplicate_inds(gl, sample_dictionary)` | Collapses replicate samples of one designation by majority call. ⚠ (P1) |
| `get_mendelian_error_counts(gl, trios)` | Counts Mendelian inconsistencies per trio. Diploid only. ⚠ (P2, a lower bound) |

## Write

| Function | What it does |
|---|---|
| `write_vcf(gl, output_path, na_rep = "./.")` | Writes a genlight object to VCF 4.0, gzipped when the path ends in `.gz`. Genotypes are unphased. |

## BioFlow interface

| Function | What it does |
|---|---|
| `apply_bioflow_modifications(gl, modifications)` | Replays the removals that BioFlow recorded in its modifications table. ⚠ (B1, B2) |

## Set-up and dependencies

| Function | What it does |
|---|---|
| `check_dependencies()` | Lists the optional packages (vcfR, randomForest, dartR.base) and what each unlocks. |
| `install_optional_deps(which, ask = interactive())` | Installs the optional packages, including the Bioconductor ones dartR.base needs. |
| `check_beagle_requirements()` | Checks Java, the Beagle JAR and vcfR. Returns a list with `ready` and a `message`. |
| `find_java()` | Path to a working Java, from `JAVA_HOME` or the system path, or `NULL`. |
| `get_java_version(java_exe = "java")` | Java version line, or `NULL`. |
| `find_beagle_jar()` | Path to the Beagle JAR shipped in `inst/extdata`, or `NULL`. |
| `get_default_beagle_threads()` | Half the cores, between 1 and 4. |

## Low-level helpers

These turn genotype strings into dosages, and the readers call them. Only `get_allelic_dosage()` is exported. The other four are internal.

| Function | Example | Result |
|---|---|---|
| `get_allelic_dosage(l, alleles, ploidity, sep = "")` | `get_allelic_dosage(c("AG", "GG"), "A/G", 2)` | 1, 2 |
| `get_all_gt_calls(alleles, ploidity, sep = "")` | `get_all_gt_calls(c("A", "B"), 2)` | `AA AB BA BB` |
| `convert_gt_to_dosage(locus, alt_allele, ploidity = 2, sep = "")` | `convert_gt_to_dosage(c("AG", "AA"), "A")` | 1, 2 |
| `genocall_to_allelic_dosage(genotype_call, alt_allele, ploidity = 2, sep = "")` | `genocall_to_allelic_dosage("AAA", "A", 3)` | 3 |
| `replace_strings_with_integers(lookup_table, strings_to_replace)` | `replace_strings_with_integers(c(A = 1, B = 2), c("B", "C"))` | 2, `NA` |

## Internal functions

Not exported, listed so that contributors can find them. `process_metadata()` flags markers in every reader. `read_tabular_geno()` reads hapmap and other tables with common missing codes. `recalc_metrics()` relies on `get_inbreeding()`. Imputation uses `freq_impute()`, `rf_impute()`, `i_rf_impute()` and `apply_imputation()`. `parse_counts_metadata()` and `compare_tags()` serve the DArTag reader. `majority_call()` and `collapse_inds()` serve the merge. `require_optional()` checks optional packages. `read_genolcsv()` and the three plot functions in `R/QC_plots_adegenet.R` are unfinished or legacy and cannot be used as they stand.
