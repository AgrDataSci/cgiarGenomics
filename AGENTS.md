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
| Read a hapmap file | `read_hapmap(path, ploidity, sep = "")` | Tab-separated, samples from column 12, calls written in full (`AG`, `AAGG`). Set `sep` when calls carry `/` or `|`. Single IUPAC codes are decoded for diploids. |
| Read a VCF | `read_vcf(path, ploidity)` | Needs vcfR. Reads phased and unphased calls. |
| Read DArTag reports (any ploidy) | `read_DArTag_count_dosage(dosage_path, counts_path, ploidity)` | Needs both the dosage report and the collapsed counts file. Check the row and column defaults against the file. |
| Read a DArTSeq SNP report | `read_DArTSeq_SNP(path, snp_id, chr_name, pos_name)` | Needs dartR.base, and SNPRelate and snpStats from Bioconductor. |
| Quality summary | `get_overall_summary(gl)` | Reads metrics stored in `gl@other`. Run `recalc_metrics(gl)` first on objects not built by this package. |
| Per-locus or per-individual metrics | `gl@other$loc.metrics`, `gl@other$ind.metrics` | Or `get_loc_missing()`, `get_ind_missing()`, `get_maf()`, `get_heterozygosity_metrics()`. |
| Filter loci and individuals | `apply_sequence_filtering(gl, steps)` | Steps are lists of `param`, `operator`, `threshold`. Order matters, metrics are recalculated after each step. A step that nothing passes stops with an error. |
| Impute missing calls | `impute_gl(gl, ploidity, method = "frequency")` | `method` has no default. Prefer `"beagle"` for diploids with dense mapped markers, `"frequency"` otherwise. Set a seed with `set.seed()` for `"frequency"` and report it. |
| Benchmark an imputation method | `set_fraction_na()` then `get_accuracy(ref, imputed, masked_gl = masked)` | Always pass `masked_gl`, and report the number of scored cells from the `n_scored` attribute. |
| Collapse replicate samples | `merge_duplicate_inds(gl, sample_dictionary)` | Individuals are renamed to their `designation_id`. Ties between replicates become missing. |
| Identity or purity between samples | `ibs_matrix_purrr(gl, ploidy)` | Matrix products in base R, despite the name. Memory grows with the square of the number of individuals. |
| Pedigree checks | `get_mendelian_error_counts(gl, trios)` | Diploid only, and the count is a lower bound. |
| Export | `write_vcf(gl, path)` | Unphased GT only, in the marker order of the object. Sort by chromosome and position first for tools that need it. |

## Known issues

Version 2.0.0 fixes the defects that the two reviews of November 2026 found in version 1.0.2. [`dev/audit.md`](dev/audit.md) records each finding with its status and [`NEWS.md`](NEWS.md) lists the results that change. Check `packageVersion("cgiarGenomics")` first.

On version 1.x, MAF, heterozygosity, Fis, Random Forest imputation, replicate merging and filters that nothing passes give wrong results. Tell the user to update before running the analysis, and to rerun earlier analyses that used them.

These stay open in version 2.0.0. Tell the user when one applies to their analysis.

- **`loc_Fis` assumes a diploid (M3).** Do not interpret or filter on `loc_Fis` in polyploids.
- **`apply_bioflow_modifications()` replays removals only (B2).** It warns when the table holds imputation rows and does not write the imputed values back.
- **Mendelian error counts are a lower bound for F1 trios (P2).** The rule is the conservative one for later generations.
- **Loci with no data disappear during imputation.** `impute_gl()` removes them with a warning, and `log` does not list them. Compare `locNames()` before and after.
- **Metrics use dense matrices (M4).** `recalc_metrics()` and `ibs_matrix_purrr()` hold the full dosage matrix in memory. Thin or split very large panels.

## Pitfalls

- **Dosage means copies of the alternative allele.** The alternative is the second allele in hapmap `alleles` and the `ALT` column in VCF. A dosage of 0 is homozygous for the reference.
- **The argument is spelled `ploidity`** in readers and imputation, and `ploidy` in metric and IBS functions. Pass the same value to both and make sure it matches the data.
- **Markers are removed while reading.** Readers drop multi-allelic markers, indels, markers without alleles and markers sharing a position. Report the counts from the reader messages, so the user knows what was lost.
- **Metrics go stale.** Any subset made outside the package (`gl[, keep]`) leaves `gl@other` describing the old object. Run `recalc_metrics()` before filtering or summarising.
- **Missingness is a rate.** All thresholds are between 0 and 1, and `filter_gl()` refuses values above 1.
- **Beagle is diploid only.** For polyploids, use frequency imputation and say that it ignores linkage. Beagle needs unique locus names, a chromosome and position for every locus, and alleles written as A, C, G or T.
- **Beagle settings are checked.** `memory` must look like `"Xmx4g"`. Do not build it from free text typed by a user.
- **Beagle has its own licence.** The bundled JAR is under the GNU GPL version 3 or later, see `inst/extdata/README.md`. Tell users who redistribute the package.
- **Sample names may identify germplasm or farmers.** Do not put real sample names or private data in examples, issues or test files.

## Changing the package

- Use `=` for assignment. Keep `<-` only where an assignment sits inside a function call, as in `expect_warning(out <- f(x))`. Tests use invented sample and marker names, built with the helpers in `tests/testthat/helper-genlight.R`.
- Prefer base R. The package imports only cli, data.table, adegenet, methods, stats and utils. Do not add a dependency for what base R does in a few lines, and discuss any new import in the pull request.
- Document with roxygen and run `devtools::document()`. It writes `NAMESPACE` and `man/`, so never edit either by hand. Export a function by tagging it `@export`. Imports are declared in `R/cgiarGenomics-package.R`. Markdown is enabled in roxygen.
- Every fix needs a test in `tests/testthat/` that fails before the fix and passes after it. Run `devtools::test()` and `devtools::check()`, and report what was actually run. The `tests` workflow runs both on every pull request. Never claim a change works without running it.
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
