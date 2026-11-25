# Code audit of cgiarGenomics

Audit of `main` at commit `fff3880` (4 November 2026), 25 November 2026. It covers the ten files in `R/`, the package metadata, the documentation and the test material.

## How this audit was done

Every function in `R/` was read line by line. Claims about adegenet and vcfR behaviour were checked against their current source code (adegenet `R/glFunctions.R` and `R/glHandle.R`, vcfR `src/extract_gt.cpp`). R was not available in the audit environment, so no code was run. Each finding below says how sure it is. **Confirmed** means the code path and the upstream source leave no doubt. **Likely** means a short test in R should confirm it before the fix. Every fix proposed here needs a `testthat` test that fails before the fix and passes after it.

## Summary

The package does what BioFlow needs on the common path, which is reading a diploid hapmap or unphased VCF, filtering by missingness and imputing by frequency. Off that path it has defects that change results without an error. The five that matter most are listed first.

| # | Severity | Area | Finding |
|---|---|---|---|
| 1 | Critical | Metrics | `get_maf()` divides by ploidy twice, so every MAF filter uses the wrong value |
| 2 | Critical | Imputation | Random Forest imputation writes factor codes instead of dosages |
| 3 | Critical | BioFlow | `apply_bioflow_modifications()` drops every individual or every locus when only one margin has removals |
| 4 | Critical | Reading | `read_DArTSeq_SNP()` drops every locus when all SNPs are mapped |
| 5 | High | Filtering | A filter step that no element passes is skipped instead of applied |
| 6 | High | Reading | Duplicated marker IDs are not renamed |
| 7 | High | Reading | `read_hapmap()` with `sep = "/"` or `"|"` always aborts |
| 8 | High | Imputation | `get_accuracy()` scores unmasked cells, so accuracy is inflated |
| 9 | High | Licence | DESCRIPTION says GPL (>= 2), LICENSE says MIT |
| 10 | High | Merge | `merge_duplicate_inds()` matches the wrong column and drops metadata |

The full list follows, grouped by area, with file and line.

## Reading data

**R1. Duplicated marker IDs are not renamed (High, confirmed).** `R/geno_read_base.R:54` tests `dplyr::row_number() %in% id_dup`, but `id_dup` is a logical vector. `%in%` coerces `TRUE` to 1, so only row 1 is renamed, and only when any duplicate exists. The genlight object then carries duplicated locus names. Fix with `id_dup_idx = which(id_dup)` and `dplyr::row_number() %in% id_dup_idx`.

**R2. `read_DArTSeq_SNP()` removes every locus when none are unmapped (Critical, confirmed).** `R/geno_input_adegenet.R:360` runs `gl[, -c(no_loc)]`. When `no_loc` is empty, `-integer(0)` selects nothing. Guard with `if (length(no_loc) > 0)`.

**R3. `read_DArTSeq_SNP()` misaligns positions after removing loci (High, likely).** adegenet subsets `@other` only when individuals are subset (`glHandle.R`, the `treatOther` block), not loci. After line 360, `gl@other$loc.metrics` keeps its original rows, and lines 363 and 364 assign positions from rows that no longer match. Take `pos` and `chrom` from the vectors already built at lines 346 and 347, subset by the same index.

**R4. `read_hapmap()` aborts with separated genotype calls (High, confirmed).** Lines 192 to 201 measure `nchar()` of the raw call. With `sep = "/"`, a diploid call `A/G` has three characters, which exceeds the ploidy of 2 and triggers the abort at line 200. Strip the separator before counting, as `read_vcf()` does at line 272.

**R5. Phased VCF calls abort under the default separator (Medium, confirmed).** `read_vcf()` defaults to `sep = "/"`. A phased call `A|G` keeps its `|`, measures three characters and aborts at line 280. Files that mix phased and unphased calls cannot be read with any single `sep`. Normalise both separators before the ploidy check.

**R6. Single-letter IUPAC hapmap calls become missing (Medium, confirmed).** TASSEL exports diploid hapmap calls as single IUPAC codes (`R`, `Y`). These never match the expanded calls from `get_all_gt_calls()`, so they turn into `NA` with only a ploidy warning. Either decode IUPAC codes or stop with a clear message.

**R7. Multi-allelic markers are reported under the wrong reason (Low, confirmed).** `R/geno_read_base.R:57` returns `NA` for any `alt` with a comma, so multi-allelic markers are counted as "without ref and alt allele data" as well as multi-allelic. The filter outcome is correct, the message is not.

**R8. `read_genolcsv()` cannot run (Low, confirmed).** Line 52 references an undefined `x`, and the function assumes hapmap columns. It is not exported, but its roxygen block carries `@export`. Remove it or finish it.

**R9. Wrong counts in progress messages (Low, confirmed).** `read_hapmap()` reports `nrows - 1` loci and `ncols - 4` individuals at lines 163 and 165. The true values are `nrows` and `ncols - 11`. `extract_dart_gt()` reports the wrong argument in its data frame check at line 389.

**R10. Duplicated DArTag markers are found but not handled (Low, confirmed).** `read_one_row_dart_gt_data()` builds `l_dupes` at line 456 and never uses it. Duplicates then multiply rows in the two `merge()` calls of `read_DArTag_count_dosage()`.

## Metrics

**M1. `get_maf()` divides by ploidy twice (Critical, confirmed).** `adegenet::glMean()` already returns allele frequencies, as `glSum() / (sum(ploidy) - glNA())` in `glFunctions.R`. `R/geno_get_metrics.R:158` multiplies by `1/ploidy` again. In a diploid the stored MAF is half the true value, and a locus fixed for the alternative allele gets MAF 0.5. Every MAF filter in `apply_sequence_filtering()` and `random_select_loci()`, and the Fis estimate, inherit the error. Fix with `alf = adegenet::glMean(gl)`. This changes results for every user, so record it in `NEWS.md`.

**M2. Heterozygosity is diluted by missing data (Medium, confirmed).** Lines 137, 138, 145 and 152 divide by all individuals or all loci, including missing calls. Divide by the number of non-missing calls instead.

**M3. Fis assumes diploid Hardy-Weinberg (Medium, confirmed).** `get_inbreeding()` uses `2pq` at line 165 for every ploidy, and inherits the MAF error from M1. For polyploids, the expected heterozygosity is `1 - p^k - q^k`, where `k` is ploidy. Report Fis as diploid only until that is in place.

**M4. Metrics are slow on large panels (Low, confirmed).** `recalc_metrics()` converts the genlight to a dense matrix three times per call, and filtering calls it after every step. `ibs_matrix_purrr()` computes each pair twice in pure R. Convert once per call, and compute only the upper triangle.

## Filtering

**F1. A filter that nothing passes is silently skipped (High, confirmed).** In `apply_sequence_filtering()`, the subset runs only `if (length(i_filt_out$index) > 0)`. When no locus passes a threshold, the step keeps every locus. In `filter_gl()` at lines 77 and 81, `filter_out` is also empty in that case. Abort with a clear message, or return an empty object, but never skip.

**F2. Thresholds above 1 are refused (Low, confirmed).** Line 71 blocks thresholds above 1, which is right for rates but undocumented, and the message says "greater or equal". It also blocks no useful case today. State the rule in the documentation.

**F3. `random_select_loci()` fails when two chromosomes tie for the most loci (Medium, confirmed).** `max_total` then has length 2, and `if (i_chr$CHROM == max_total)` is an error in R 4.2 and later. Use `max_total[1]`. The function also requires `seed` with no default.

## Imputation

**I1. Random Forest imputation writes factor codes (Critical, confirmed).** `rf_impute()` passes the target as a factor, so `randomForest()` classifies and returns a factor. Line 282 assigns that factor into a numeric matrix, which stores its integer codes. With levels `0`, `1` and `2`, every imputed dosage is shifted by one, and dosage 2 becomes 3. Convert with `as.numeric(as.character(imputed_target))` before returning.

**I2. Accuracy counts cells that were never masked (High, confirmed).** `get_accuracy()` compares every non-missing reference cell, most of which were never masked and match by construction. Reported accuracy is close to 1 whatever the method. `set_fraction_na()` should return the masked positions, and `get_accuracy()` should score only those.

**I3. `set_fraction_na()` drops names, positions and alleles (Medium, confirmed).** The genlight built at line 419 keeps only chromosomes. Beagle then cannot run on the masked object, and results cannot be matched by name. It can also pick cells that were already missing. Keep all slots and sample from non-missing cells.

**I4. `impute_beagle()` fails when nothing is missing (Medium, confirmed).** `imp_dict` is created only inside `if (length(idx_na) > 0)` at line 587, so the return at the end errors. Initialise it as `list()`.

**I5. Beagle input is not sorted by position (Medium, likely).** `write_vcf()` writes markers in genlight order with no `##contig` lines. Beagle requires markers sorted within each chromosome and rejects unsorted input. Sort before writing and map the results back to the original order.

**I6. `impute_gl()` prints debug output and miscounts (Low, confirmed).** Line 327 prints `dots` on every call. Line 329 subtracts all-missing loci after they were already removed, so the reported count is too low. Loci removed for being all missing are also gone from the output without being listed in `log`.

**I7. `ploidity` can disagree with the object (Low, confirmed).** `impute_gl()` takes `ploidity` as an argument instead of reading `adegenet::ploidy(gl)`. A wrong value gives wrong all-missing detection and wrong frequency draws.

## Individuals and pedigree

**P1. `merge_duplicate_inds()` matches the wrong column (High, confirmed).** Line 34 compares `sample_id` with `single_samples`, which holds `designation_id` values, and line 50 does the same. The function works only when single samples have `sample_id == designation_id`. When no designation is duplicated, `do.call(rbind, list())` returns `NULL` and line 46 fails. The merged object loses chromosomes, positions, alleles and metrics, so filtering it afterwards fails. Compare on `designation_id`, handle the no-duplicate case and rebuild the full genlight with `recalc_metrics()`.

**P2. Mendelian errors miss opposite homozygotes (Medium, confirmed).** `get_mendelian_error_counts()` flags an error only when both parents share a homozygote (line 328). A child of `0 x 2` parents must be heterozygous in an F1, and that case is never counted. The documentation also names the argument `trio` and promises a vector, while the function takes `trios` and returns a data frame.

## BioFlow interface

**B1. `apply_bioflow_modifications()` can drop a whole margin (Critical, confirmed).** Lines 56 and 57 index with `[-ind_out]` and `[-loc_out]`. When one of them is empty, the result is empty, so a modification table with only locus removals removes every individual, and the reverse. Use `setdiff(seq_len(n), ind_out)`.

**B2. Imputation modifications are never applied (High, confirmed).** The function returns at line 61, so the imputation branch below it is dead code. It also calls `filter()` without the `dplyr::` prefix and builds a genlight without ploidy. If BioFlow relies on this branch, imputed values recorded in BioFlow are not replayed.

## Packaging and documentation

**D1. Conflicting licences (High, confirmed).** `DESCRIPTION` declares `GPL (>= 2)` and `LICENSE` holds the MIT text. The owners decide which one applies. For MIT, use `License: MIT + file LICENSE` and the two-line LICENSE that R expects, with the full text in `LICENSE.md`.

**D2. NAMESPACE is written by hand while documentation uses roxygen (High, confirmed).** Twenty exported functions have no help page, among them `get_maf()`, `merge_duplicate_inds()` and `read_DArTag_count_dosage()`. Ten help pages document unexported functions. `roxygen2` will not overwrite a NAMESPACE it did not write, so the two drift apart. Eight functions also carry `@export` in roxygen while NAMESPACE keeps them internal, among them `get_all_gt_calls()` and `i_freq_impute()`. This branch adds roxygen blocks to every exported function and leaves those tags as they were. Decide which helpers to export, then move NAMESPACE to roxygen.

**D3. Examples cannot run (High, confirmed).** Several examples load `data(example_genlight)`, which does not exist, and `read_vcf()` downloads from the old upstream URL. `R CMD check` fails on them. Ship a small dataset in `inst/extdata` and use it in all examples.

**D4. No automated tests or CI (High, confirmed).** `tests/` holds data only, `.gitignore` excludes `tests/*.R`, and the two scripts in `test_scripts/` print results instead of asserting them. Add `tests/testthat/` with one test per finding above, and a GitHub Actions `R-CMD-check` workflow.

**D5. Undeclared and unneeded dependencies (Medium, confirmed).** `plot_overall_missingness()` calls `ggpubr`, which is not declared. The three plot functions in `R/QC_plots_adegenet.R` use an older object layout (`geno_obj$geno`) and are not exported, yet they pull ggplot2, ggrepel and scales into Imports. `parallel` is used under Suggests without a `requireNamespace()` check. Remove or update the plots and move their packages to Suggests.

**D6. Metadata is out of date (Medium, fixed in PR #2).** `DESCRIPTION` now has `Authors@R` (Sergio Andres Cruz, author and maintainer), a current date, a full description, and `URL` and `BugReports`. There is still no `CITATION.cff` and no DOI.

**D7. Beagle JAR in the package (Low, confirmed).** `inst/extdata/beagle.27Feb25.75f.jar` (310 KB) is fine on GitHub and Apache 2.0 allows redistribution, provided the licence and notice travel with it. Add the Apache licence text next to the JAR. CRAN will not accept the JAR inside the package.

**D8. Spelling of the ploidy argument (Low).** Reading and imputation functions use `ploidity`, while metric and IBS functions use `ploidy`. Accept `ploidy` as an alias in the next minor version and deprecate `ploidity` later, so BioFlow does not break.

## Proposed order of work

1. Fix M1, I1, B1 and R2, each with a test, and release as a patch with a `NEWS.md` entry that names the changed results.
2. Fix F1, R1, R4, I2 and P1 with tests.
3. Settle the licence (D1), move NAMESPACE to roxygen (D2), add example data, tests and CI (D3, D4).
4. Clean dependencies and metadata (D5, D6, D7), then the remaining medium and low findings.
