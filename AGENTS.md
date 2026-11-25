# Instructions for AI agents

cgiarGenomics is an R package that reads SNP marker data (hapmap, VCF, DArTSeq, DArTag), converts it to allele dosages in an adegenet genlight object, and runs quality control, filtering, imputation and VCF export. It is the genomic data engine behind Bioflow, the CGIAR breeding analytics pipeline. People share this repository with an assistant either to **process their marker data** or to **change the package**. Find out which before you start.

## Before you write any code

1. Check the installed version with `packageVersion("cgiarGenomics")` and read the help of every function you use. Do not rely on memory for argument names.
2. Read [`REFERENCE.md`](REFERENCE.md) for what each function takes and returns, and the "Known issues" below before recommending a function.
3. Look for an existing function before writing your own. The readers already handle allele parsing, ploidy, duplicated positions and missing codes, and a hand-made parser is the most common source of errors.
4. Run `check_dependencies()` first when the user reads VCF or DArTSeq files or imputes with Random Forest or Beagle. Those need optional packages, and Beagle needs Java (`check_beagle_requirements()`).

## Which function for which task

| Task | Function | Notes |
|---|---|---|
| Read a hapmap file | `read_hapmap(path, ploidity, sep = "")` | Tab-separated, samples from column 12, calls written in full (`AG`, `AAGG`). Keep `sep = ""`. |
| Read a VCF | `read_vcf(path, ploidity, sep = "/")` | Needs vcfR. Use `sep = "\\|"` for phased files. |
| Read DArTag reports (any ploidy) | `read_DArTag_count_dosage(dosage_path, counts_path, ploidity)` | Needs both the dosage report and the collapsed counts file. Check the row and column defaults against the file. |
| Read a DArTSeq SNP report | `read_DArTSeq_SNP(path, snp_id, chr_name, pos_name)` | Needs dartR.base, and SNPRelate and snpStats from Bioconductor. |
| Quality summary | `get_overall_summary(gl)` | Reads metrics stored in `gl@other`. Run `recalc_metrics(gl)` first on objects not built by this package. |
| Per-locus or per-individual metrics | `gl@other$loc.metrics`, `gl@other$ind.metrics` | Or `get_loc_missing()`, `get_ind_missing()`, `get_maf()`, `get_heterozygosity_metrics()`. |
| Filter loci and individuals | `apply_sequence_filtering(gl, steps)` | Steps are lists of `param`, `operator`, `threshold`. Order matters, metrics are recalculated after each step. |
| Impute missing calls | `impute_gl(gl, ploidity, method = "frequency")` | `method` has no default. Prefer `"beagle"` for diploids with dense markers sorted by position (known issue I5), `"frequency"` otherwise. |
| Benchmark an imputation method | `set_fraction_na()` then `get_accuracy()` | See known issue I2 before reporting accuracy. |
| Collapse replicate samples | `merge_duplicate_inds(gl, sample_dictionary)` | See known issue P1. |
| Identity or purity between samples | `ibs_matrix_purrr(gl, ploidy)` | Slow, quadratic in individuals. Thin loci first with `random_select_loci()`. |
| Pedigree checks | `get_mendelian_error_counts(gl, trios)` | Diploid only, and the count is a lower bound. |
| Export | `write_vcf(gl, path)` | Unphased GT only. Sort by chromosome and position before using the file with Beagle or other tools. |

## Known issues

The audit in [`dev/audit.md`](dev/audit.md) lists every defect found in version 1.0.2, with location and fix. Version 1.0.3 fixes M1 (MAF), I1 (Random Forest dosages), B1 (`apply_bioflow_modifications()`) and R2 (DArTSeq reading). Check `packageVersion("cgiarGenomics")` first. On 1.0.2 or earlier, do not filter on `maf`, do not use `method = "random_forest"`, and tell the user to update. Until the remaining ones are fixed, act as follows, and tell the user when one applies to their analysis.

- **A filter step that nothing passes is skipped (F1).** After `apply_sequence_filtering()`, check `filt_log` and report any step whose `index` is empty.
- **Duplicated marker IDs survive (R1).** Check `anyDuplicated(adegenet::locNames(gl))` after reading.
- **Imputation accuracy is inflated (I2).** Score only the cells that `set_fraction_na()` masked, by comparing the matrices yourself.
- **`apply_bioflow_modifications()` ignores imputation rows (B2).** It replays removals only. Imputed values recorded in Bioflow are not written back.
- **Separators in hapmap calls abort the reader (R4).** Remove `/` or `|` from hapmap calls before reading, or keep `sep = ""`.

## Pitfalls

- **Dosage means copies of the alternative allele.** The alternative is the second allele in hapmap `alleles` and the `ALT` column in VCF. A dosage of 0 is homozygous for the reference.
- **The argument is spelled `ploidity`** in readers and imputation, and `ploidy` in metric and IBS functions. Pass the same value to both and make sure it matches the data.
- **Markers are removed while reading.** Readers drop multi-allelic markers, indels, markers without alleles and markers sharing a position. Report the counts from the reader messages, so the user knows what was lost.
- **Metrics go stale.** Any subset made outside the package (`gl[, keep]`) leaves `gl@other` describing the old object. Run `recalc_metrics()` before filtering or summarising.
- **Missingness is a rate.** All thresholds are between 0 and 1, and `filter_gl()` refuses values above 1.
- **Beagle is diploid only.** For polyploids, use frequency imputation and say that it ignores linkage.
- **Sample names may identify germplasm or farmers.** Do not put real sample names or private data in examples, issues or test files.

## Changing the package

- Follow the style of the file you edit.
- `NAMESPACE` is written by hand. When you export a function, add it to `NAMESPACE` as well as tagging it `@export`.
- Document with roxygen and run `devtools::document()`. Never edit files in `man/` by hand. Markdown is enabled in roxygen.
- Every fix needs a test in `tests/testthat/` that fails before the fix and passes after it. Run `devtools::test()` and `devtools::check()`, and report what was actually run. Never claim a change works without running it.
- Changes that alter results (dosages, metrics, filters, imputed values) go in `NEWS.md` with the reason.
- Bioflow calls these functions directly. Keep function names, argument names and return structures unchanged, or agree the change with the Bioflow team and record it in `NEWS.md`.
- When you fix a finding from `dev/audit.md`, mark it as fixed there with the pull request number, and remove the "Known issue" section from the function's roxygen.
- Work on a branch and open a pull request. Do not merge without review and passing checks.

## Where to find more

- Function reference: [`REFERENCE.md`](REFERENCE.md)
- Audit and fix plan: [`dev/audit.md`](dev/audit.md)
- Bioflow: https://github.com/Breeding-Analytics/bioflow
- adegenet genlight class: https://github.com/thibautjombart/adegenet
- Issues: https://github.com/Breeding-Analytics/cgiarGenomics/issues
