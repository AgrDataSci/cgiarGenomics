# Code audit of cgiarGenomics

Audit of `main` at commit `fff3880` (4 November 2026), 25 November 2026. It covers the ten files in `R/`, the package metadata, the documentation and the test material.

## How this audit was done

Every function in `R/` was read line by line. Claims about adegenet and vcfR behaviour were checked against their current source code (adegenet `R/glFunctions.R` and `R/glHandle.R`, vcfR `src/extract_gt.cpp`). R was not available in the audit environment, so no code was run. Each finding below says how sure it is. **Confirmed** means the code path and the upstream source leave no doubt. **Likely** means a short test in R should confirm it before the fix. Every fix proposed here needs a `testthat` test that fails before the fix and passes after it.

## Status of fixes

PR #3 on the AgrDataSci fork (25 November 2026) brings the package to version 1.0.3. It fixes the findings marked "fixed in PR #3" below and the additional findings of the external review, each with tests in `tests/testthat/`. The `tests` workflow runs the unit tests and `R CMD check` on GitHub Actions for every pull request, with vcfR, randomForest, dartR.base and Java installed, so the Beagle, Random Forest and DArTSeq paths are exercised.

The tests were written against the fixed code. They were not run against version 1.0.2, so their failure before the fix rests on reading the code and on the external review, which ran version 1.0.2 and observed the defects.

| Still open | Finding |
|---|---|
| `loc_Fis` assumes a diploid at every ploidy | M3 |
| Metrics build dense matrices on large panels | M4 |
| Mendelian errors use the later-generation rule only | P2 |
| `apply_bioflow_modifications()` does not replay imputation | B2 |
| Loci removed for having no data are not listed in the imputation log | I6 |
| `ploidity` and `ploidy` spellings | D8 |
| GPL text for the bundled Beagle JAR, and the decision to keep bundling it | D7 |
| Provenance and reuse terms of the example files in `tests/` | E8 |
| A DOI and a tagged release | D6 |
| Windows job in CI, large-panel benchmarks, end-to-end run through Bioflow | E9 |

The summary table below describes version 1.0.2.

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

**R1. Duplicated marker IDs are not renamed (High, confirmed, fixed in PR #3).** `R/geno_read_base.R:54` tests `dplyr::row_number() %in% id_dup`, but `id_dup` is a logical vector. `%in%` coerces `TRUE` to 1, so only row 1 is renamed, and only when any duplicate exists. The genlight object then carries duplicated locus names. Fix with `id_dup_idx = which(id_dup)` and `dplyr::row_number() %in% id_dup_idx`.

**R2. `read_DArTSeq_SNP()` removes every locus when none are unmapped (Critical, confirmed, fixed in PR #3).** `R/geno_input_adegenet.R:360` runs `gl[, -c(no_loc)]`. When `no_loc` is empty, `-integer(0)` selects nothing. Guard with `if (length(no_loc) > 0)`.

**R3. `read_DArTSeq_SNP()` misaligns positions after removing loci (High, likely, fixed in PR #3).** adegenet subsets `@other` only when individuals are subset (`glHandle.R`, the `treatOther` block), not loci. After line 360, `gl@other$loc.metrics` keeps its original rows, and lines 363 and 364 assign positions from rows that no longer match. Take `pos` and `chrom` from the vectors already built at lines 346 and 347, subset by the same index.

**R4. `read_hapmap()` aborts with separated genotype calls (High, confirmed, fixed in PR #3).** Lines 192 to 201 measure `nchar()` of the raw call. With `sep = "/"`, a diploid call `A/G` has three characters, which exceeds the ploidy of 2 and triggers the abort at line 200. Strip the separator before counting, as `read_vcf()` does at line 272.

**R5. Phased VCF calls abort under the default separator (Medium, confirmed, fixed in PR #3).** `read_vcf()` defaults to `sep = "/"`. A phased call `A|G` keeps its `|`, measures three characters and aborts at line 280. Files that mix phased and unphased calls cannot be read with any single `sep`. Normalise both separators before the ploidy check.

**R6. Single-letter IUPAC hapmap calls become missing (Medium, confirmed, fixed in PR #3).** TASSEL exports diploid hapmap calls as single IUPAC codes (`R`, `Y`). These never match the expanded calls from `get_all_gt_calls()`, so they turn into `NA` with only a ploidy warning. Either decode IUPAC codes or stop with a clear message.

**R7. Multi-allelic markers are reported under the wrong reason (Low, confirmed, fixed in PR #3).** `R/geno_read_base.R:57` returns `NA` for any `alt` with a comma, so multi-allelic markers are counted as "without ref and alt allele data" as well as multi-allelic. The filter outcome is correct, the message is not.

**R8. `read_genolcsv()` cannot run (Low, confirmed, fixed in PR #3).** Line 52 references an undefined `x`, and the function assumes hapmap columns. It is not exported, but its roxygen block carries `@export`. Remove it or finish it.

**R9. Wrong counts in progress messages (Low, confirmed, fixed in PR #3).** `read_hapmap()` reports `nrows - 1` loci and `ncols - 4` individuals at lines 163 and 165. The true values are `nrows` and `ncols - 11`. `extract_dart_gt()` reports the wrong argument in its data frame check at line 389.

**R10. Duplicated DArTag markers are found but not handled (Low, confirmed, fixed in PR #3).** `read_one_row_dart_gt_data()` builds `l_dupes` at line 456 and never uses it. Duplicates then multiply rows in the two `merge()` calls of `read_DArTag_count_dosage()`.

## Metrics

**M1. `get_maf()` divides by ploidy twice (Critical, confirmed, fixed in PR #3).** `adegenet::glMean()` already returns allele frequencies, as `glSum() / (sum(ploidy) - glNA())` in `glFunctions.R`. `R/geno_get_metrics.R:158` multiplies by `1/ploidy` again. In a diploid the stored MAF is half the true value, and a locus fixed for the alternative allele gets MAF 0.5. Every MAF filter in `apply_sequence_filtering()` and `random_select_loci()`, and the Fis estimate, inherit the error. Fix with `alf = adegenet::glMean(gl)`. This changes results for every user, so record it in `NEWS.md`.

**M2. Heterozygosity is diluted by missing data (Medium, confirmed, fixed in PR #3).** Lines 137, 138, 145 and 152 divide by all individuals or all loci, including missing calls. Divide by the number of non-missing calls instead.

**M3. Fis assumes diploid Hardy-Weinberg (Medium, confirmed, open).** `get_inbreeding()` uses `2pq` at line 165 for every ploidy, and inherits the MAF error from M1. For polyploids, the expected heterozygosity is `1 - p^k - q^k`, where `k` is ploidy. Report Fis as diploid only until that is in place.

**M4. Metrics are slow on large panels (Low, confirmed, partly fixed in PR #5).** `ibs_matrix_purrr()` now uses matrix products. `recalc_metrics()` still builds dense matrices. `recalc_metrics()` converts the genlight to a dense matrix three times per call, and filtering calls it after every step. `ibs_matrix_purrr()` computes each pair twice in pure R. Convert once per call, and compute only the upper triangle.

## Filtering

**F1. A filter that nothing passes is silently skipped (High, confirmed, fixed in PR #3).** In `apply_sequence_filtering()`, the subset runs only `if (length(i_filt_out$index) > 0)`. When no locus passes a threshold, the step keeps every locus. In `filter_gl()` at lines 77 and 81, `filter_out` is also empty in that case. Abort with a clear message, or return an empty object, but never skip.

**F2. Thresholds above 1 are refused (Low, confirmed, documented in PR #3).** Line 71 blocks thresholds above 1, which is right for rates but undocumented, and the message says "greater or equal". It also blocks no useful case today. State the rule in the documentation.

**F3. `random_select_loci()` fails when two chromosomes tie for the most loci (Medium, confirmed, fixed in PR #3).** `max_total` then has length 2, and `if (i_chr$CHROM == max_total)` is an error in R 4.2 and later. Use `max_total[1]`. The function also requires `seed` with no default.

## Imputation

**I1. Random Forest imputation writes factor codes (Critical, confirmed, fixed in PR #3).** `rf_impute()` passes the target as a factor, so `randomForest()` classifies and returns a factor. Line 282 assigns that factor into a numeric matrix, which stores its integer codes. With levels `0`, `1` and `2`, every imputed dosage is shifted by one, and dosage 2 becomes 3. Convert with `as.numeric(as.character(imputed_target))` before returning.

**I2. Accuracy counts cells that were never masked (High, confirmed, fixed in PR #3).** `get_accuracy()` compares every non-missing reference cell, most of which were never masked and match by construction. Reported accuracy is close to 1 whatever the method. `set_fraction_na()` should return the masked positions, and `get_accuracy()` should score only those.

**I3. `set_fraction_na()` drops names, positions and alleles (Medium, confirmed, fixed in PR #3).** The genlight built at line 419 keeps only chromosomes. Beagle then cannot run on the masked object, and results cannot be matched by name. It can also pick cells that were already missing. Keep all slots and sample from non-missing cells.

**I4. `impute_beagle()` fails when nothing is missing (Medium, confirmed, fixed in PR #3).** `imp_dict` is created only inside `if (length(idx_na) > 0)` at line 587, so the return at the end errors. Initialise it as `list()`.

**I5. Beagle input is not sorted by position (Medium, likely, fixed in PR #3).** `write_vcf()` writes markers in genlight order with no `##contig` lines. Beagle requires markers sorted within each chromosome and rejects unsorted input. Sort before writing and map the results back to the original order.

**I6. `impute_gl()` prints debug output and miscounts (Low, confirmed, fixed in PR #3, except that removed loci are not listed in `log`).** Line 327 prints `dots` on every call. Line 329 subtracts all-missing loci after they were already removed, so the reported count is too low. Loci removed for being all missing are also gone from the output without being listed in `log`.

**I7. `ploidity` can disagree with the object (Low, confirmed, PR #3 adds a warning).** `impute_gl()` takes `ploidity` as an argument instead of reading `adegenet::ploidy(gl)`. A wrong value gives wrong all-missing detection and wrong frequency draws.

## Individuals and pedigree

**P1. `merge_duplicate_inds()` matches the wrong column (High, confirmed, fixed in PR #3).** Line 34 compares `sample_id` with `single_samples`, which holds `designation_id` values, and line 50 does the same. The function works only when single samples have `sample_id == designation_id`. When no designation is duplicated, `do.call(rbind, list())` returns `NULL` and line 46 fails. The merged object loses chromosomes, positions, alleles and metrics, so filtering it afterwards fails. Compare on `designation_id`, handle the no-duplicate case and rebuild the full genlight with `recalc_metrics()`.

**P2. Mendelian errors miss opposite homozygotes (Medium, confirmed, open, documented as the intended conservative rule).** `get_mendelian_error_counts()` flags an error only when both parents share a homozygote (line 328). A child of `0 x 2` parents must be heterozygous in an F1, and that case is never counted. The documentation also names the argument `trio` and promises a vector, while the function takes `trios` and returns a data frame.

## BioFlow interface

**B1. `apply_bioflow_modifications()` can drop a whole margin (Critical, confirmed, fixed in PR #3).** Lines 56 and 57 index with `[-ind_out]` and `[-loc_out]`. When one of them is empty, the result is empty, so a modification table with only locus removals removes every individual, and the reverse. Use `setdiff(seq_len(n), ind_out)`.

**B2. Imputation modifications are never applied (High, confirmed, open, PR #3 removes the dead code and adds a warning).** The function returns at line 61, so the imputation branch below it is dead code. It also calls `filter()` without the `dplyr::` prefix and builds a genlight without ploidy. If BioFlow relies on this branch, imputed values recorded in BioFlow are not replayed.

## Packaging and documentation

**D1. Conflicting licences (High, confirmed, fixed in PR #3).** `DESCRIPTION` declares `GPL (>= 2)` and `LICENSE` holds the MIT text. The owners decide which one applies. For MIT, use `License: MIT + file LICENSE` and the two-line LICENSE that R expects, with the full text in `LICENSE.md`.

**D2. NAMESPACE is written by hand while documentation uses roxygen (High, confirmed, fixed in PR #3).** Twenty exported functions have no help page, among them `get_maf()`, `merge_duplicate_inds()` and `read_DArTag_count_dosage()`. Ten help pages document unexported functions. `roxygen2` will not overwrite a NAMESPACE it did not write, so the two drift apart. Eight functions also carry `@export` in roxygen while NAMESPACE keeps them internal, among them `get_all_gt_calls()` and `i_freq_impute()`. This branch adds roxygen blocks to every exported function and leaves those tags as they were. Decide which helpers to export, then move NAMESPACE to roxygen.

**D3. Examples cannot run (High, confirmed, fixed in PR #3).** Several examples load `data(example_genlight)`, which does not exist, and `read_vcf()` downloads from the old upstream URL. `R CMD check` fails on them. Ship a small dataset in `inst/extdata` and use it in all examples.

**D4. No automated tests or CI (High, confirmed, fixed in PR #3).** `tests/` holds data only, `.gitignore` excludes `tests/*.R`, and the two scripts in `test_scripts/` print results instead of asserting them. Add `tests/testthat/` with one test per finding above, and a GitHub Actions `R-CMD-check` workflow.

**D5. Undeclared and unneeded dependencies (Medium, confirmed, fixed in PR #3).** `plot_overall_missingness()` calls `ggpubr`, which is not declared. The three plot functions in `R/QC_plots_adegenet.R` use an older object layout (`geno_obj$geno`) and are not exported, yet they pull ggplot2, ggrepel and scales into Imports. `parallel` is used under Suggests without a `requireNamespace()` check. Remove or update the plots and move their packages to Suggests.

**D6. Metadata is out of date (Medium, fixed in PR #2).** `DESCRIPTION` now has `Authors@R` (Sergio Andres Cruz, author and maintainer, and Lorena Guimarães Batista, author), a current date, a full description, and `URL` and `BugReports`. PR #3 adds `CITATION.cff`. There is still no DOI or tagged release.

**D7. Beagle JAR in the package (High, corrected in PR #3).** The first version of this audit said Beagle is under Apache 2.0. That was wrong. Its authors release Beagle 5.5 under the GNU General Public License, version 3 or later. `inst/extdata/README.md` and the README now say so, with the source, version and SHA-256 of `beagle.27Feb25.75f.jar`. The GPL text is still not shipped next to the JAR, and the rights holder has not yet decided whether to keep bundling it or to fetch it at set-up. CRAN will not accept the JAR inside the package.

**D8. Spelling of the ploidy argument (Low, open).** Reading and imputation functions use `ploidity`, while metric and IBS functions use `ploidy`. Accept `ploidy` as an alias in the next minor version and deprecate `ploidity` later, so BioFlow does not break.

## External review of 25 November 2026

An external review installed version 1.0.2 (R 4.5.0 on Debian 13), ran `R CMD check` and exercised the example files and small panels with known answers. It confirmed M1, I1, F1, R1, R2, R4, I2, B1 and P1 by running them, and added the findings below. On the diploid example file, at a MAF cutoff of 5%, version 1.0.2 excluded 63 eligible markers of 3,553 and kept one that should have failed.

**E1. Beagle settings reached a shell command (High, fixed in PR #3).** `impute_beagle()` pasted `memory` into a command line run by `system()`, and Bioflow passes that value from a text field. `memory` must now match a Java memory setting, the other settings must be whole numbers, and the arguments go to `system2()` one by one with quoted paths. Invalid settings stop before any process starts. Bioflow should still bound memory, threads and run time on the server.

**E2. Frequency imputation failed on a complete panel (High, fixed in PR #3).** `freq_impute()` indexed with an empty list when no call was missing. `impute_gl()` now returns the object with an empty log.

**E3. Random Forest failed with one predictor and could not be reproduced from its seed (Medium, fixed in PR #3).** A single flanking marker dropped the matrix to a vector. The predictors are filled by random draws, which the method seed did not cover. `impute_gl()` now also stops when any method returns a value outside the dosage range, or tries to change a call that was not missing.

**E4. Beagle results were matched by position in the matrix (High, fixed in PR #3).** The reader drops loci that share a position or lack A, C, G or T alleles, so the Beagle output could have fewer loci than the input, and values then landed in the wrong cells. Results are now matched by individual and locus name, and the function stops when a locus is missing from the output or a called genotype changed.

**E5. Missing polyploid genotypes were exported as diploid (Medium, fixed in PR #3).** `write_vcf()` wrote `./.` at every ploidy. The default now follows the ploidy of each individual. `read_vcf()` reads calls made of dots only as missing at any ploidy.

**E6. A VCF with one variant could not be read (Low, fixed in PR #3).** `vcfR::getFIX()` returns a vector in that case. Hapmap and VCF readers now keep the matrix shape with one sample or one locus.

**E7. `R CMD check` failed (High, fixed in PR #3).** The check ended with 1 error, 4 warnings and 1 note. Examples called unexported functions and a dataset that does not exist, 20 of 34 exported functions had no help page, ggpubr and BiocManager were used without being declared, and the code needs R 4.1 while `DESCRIPTION` said 3.5. roxygen now writes `NAMESPACE` and `man/`, the unusable plot helpers and `read_genolcsv()` are removed, and `DESCRIPTION` declares R 4.1 and BiocManager.

**E8. No provenance for the example files (Medium, open).** The files in `tests/` hold genotype calls with provider-style sample, plate and project identifiers, and nothing records their source or reuse terms. The data stewards need to confirm that they can be redistributed, or replace them with synthetic panels. The new unit tests use invented names, except the DArTSeq test, which reads the example report. Readers also print the first sample names in their messages.

**E9. Checks not yet done (Medium, open).** No Windows job, no benchmark of memory and time at realistic panel sizes, and no end-to-end run through the Bioflow preview, report and export entry points.

## Next steps

1. Bioflow maintainers review `NEWS.md`, since MAF, heterozygosity, Fis, Random Forest, replicate merging and no-pass filters now give different results or stop with an error, and decide how the application reports an empty filter.
2. Settle E8 and D7 before any external release.
3. Rerun a sample of past analyses that used a positive MAF cutoff, thinning, Fis, Random Forest or replicate merging, and report what changed.
4. Decide on M3, B2 and D8, then add `CITATION.cff`, a DOI and a tagged release.
