# cgiarGenomics 1.0.4

## Fewer dependencies

* The package now imports cli, data.table, adegenet, methods, stats and utils only. dplyr, tidyr, purrr, stringr, magrittr, Matrix, rlang and glue are no longer needed. Each call was rewritten in base R, with the same results.
* Assignments in the package code use `=` instead of `<-`.

## Faster

* `ibs_matrix_purrr()` computes identity by state with matrix products instead of a loop over every pair, which ran each pair twice. Results are unchanged for whole dosages. Other dosages fall back to a loop that visits each pair once. The name is kept for compatibility.

## Results that change

* `random_select_loci()` returns the selected loci in the order of the object, instead of grouped by chromosome in random order. For a given `seed`, the loci drawn may differ from earlier versions, because the random draws are made in base R.

# cgiarGenomics 1.0.3

This version follows two reviews of version 1.0.2 in November 2026, the code audit in `dev/audit.md` and an external review that ran the package. Many fixes change results, so numbers produced with earlier versions will differ. Finding numbers refer to `dev/audit.md`.

## Results that change

* `get_maf()` and the `maf` column stored by `recalc_metrics()` now hold the true minor allele frequency. Earlier versions divided the allele frequency by ploidy a second time, so a diploid MAF was half the true value and a locus fixed for the alternative allele got 0.5. Every filter on `maf` in `apply_sequence_filtering()` and `random_select_loci()` now keeps a different set of loci at the same threshold. `ov_maf` in `get_overall_summary()` and `loc_Fis` change as well (M1).
* Heterozygosity is now the share of heterozygotes among called genotypes. Earlier versions divided by all individuals or all loci, missing calls included, which lowered `loc_het`, `ind_het`, `ov_het` and raised `loc_Fis` wherever data were missing. A locus or individual with no call now gets `NaN` (M2).
* `impute_gl(method = "random_forest")` now writes the predicted dosages. Earlier versions wrote the integer codes of the predicted factor, so a predicted 0 was stored as 1 and a predicted 2 as 3 (I1).
* `apply_sequence_filtering()` now stops with an error of class `cgiarGenomics_empty_filter` when no locus or no individual passes a step. Earlier versions skipped the step and returned unfiltered data (F1). `filter_gl()` now lists every name in `filter_out` in that case.
* `merge_duplicate_inds()` now groups individuals by `designation_id`. Earlier versions compared sample IDs with designation IDs and worked only when the two were equal. Individuals are named by designation, and the result keeps chromosomes, positions and alleles and has metrics. The function no longer fails when no designation has replicates (P1).
* `apply_bioflow_modifications()` keeps every individual when the table removes only loci, and every locus when it removes only individuals (B1). It warns when the table holds imputation rows, which it does not apply (B2).
* Duplicated marker IDs are renamed to `chrom_pos` when read. Earlier versions renamed the first marker instead (R1).
* `get_accuracy()` takes a new argument `masked_gl` and then scores only the cells that were masked. Without it the function behaves as before and warns that accuracy is inflated. The result carries the number of scored cells in the attribute `n_scored` (I2).
* `set_fraction_na()` masks called genotypes only, keeps names, chromosomes, positions and alleles, and stores the masked cells in `other$masked_cells` (I3). The cells it picks for a given seed differ from earlier versions.
* `write_vcf()` writes a missing genotype at the ploidy of the individual, `./././.` for a tetraploid. The default of `na_rep` changes from `"./."` to `NULL`. Diploid output is unchanged (E5).
* `impute_beagle()` sorts markers by chromosome and position before running Beagle and matches the result back by individual and locus name. Earlier versions wrote markers in object order and matched results by position in the matrix, which misplaced values when Beagle or the reader dropped a locus (I5, E4).
* Random Forest imputation sets its seed before filling the predictors, so two runs with the same `seed` agree (E3).

## Bug fixes

* `read_DArTSeq_SNP()` keeps all loci when every SNP is mapped. Earlier versions returned an object with no loci. Positions and chromosomes come from the loci that were kept (R2, R3).
* `read_hapmap()` reads calls written with `/` or `|` when `sep` is set, and decodes diploid files with one IUPAC code per call. Its progress message reports the right numbers of loci and individuals (R4, R6, R9).
* `read_vcf()` reads files that mix phased and unphased calls, files with one variant, and missing calls at any ploidy. `sep` is kept for compatibility and no longer used (R5, E6).
* `impute_gl()` returns an object with no missing calls unchanged, with an empty log. Earlier versions failed with the frequency method and with Beagle (E2, I4). It no longer prints its extra arguments, reports the right number of missing calls, warns when `ploidity` disagrees with the object, and stops when a method returns a value outside the dosage range (I6, I7).
* Random Forest imputation works with a single flanking marker (E3).
* `random_select_loci()` works when two chromosomes tie for the most loci, and `seed` defaults to `NULL` (F3).
* `read_DArTag_count_dosage()` keeps the first row of a marker ID that appears more than once in the dosage report, with a warning (R10).
* Reader messages no longer count multi-allelic markers as markers without alleles (R7).

## Security

* `impute_beagle()` checks its settings before starting any process. `memory` must be a Java memory setting such as `"Xmx4g"`, and `burnin`, `iterations`, `seed` and `nthreads` must be whole numbers. The Beagle command is run with `system2()` and quoted paths. Earlier versions pasted `memory` into a shell command (E1).

## Package

* The licence is MIT. `DESCRIPTION` declares `MIT + file LICENSE`, and the full text is in `LICENSE.md` (D1).
* The bundled Beagle 5.5 JAR is under the GNU General Public License, version 3 or later, and not under Apache 2.0 as earlier notes said. `inst/extdata/README.md` gives its source, version and checksum (D7).
* roxygen now writes `NAMESPACE` and `man/`. Every exported function has a help page with a runnable example (D2, D3).
* Removed the unexported `read_genolcsv()` and the plot helpers `plot_missingnes_by_marker()`, `plot_MAF()` and `plot_overall_missingness()`, none of which could run. ggplot2, ggrepel and scales are no longer imported (R8, D5).
* The package requires R 4.1 or later, which the code already needed. BiocManager is declared in `Suggests`.
* Unit tests in `tests/testthat/` cover every fix, and a GitHub Actions workflow runs them and `R CMD check` on each pull request (D4).
