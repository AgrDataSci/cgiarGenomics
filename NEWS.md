# cgiarGenomics 1.0.3

This version fixes the four critical findings of the code audit in `dev/audit.md`. Three of them change results, so numbers produced with earlier versions will differ.

## Results that change

* `get_maf()` and the `maf` column stored by `recalc_metrics()` now hold the true minor allele frequency. Earlier versions divided the allele frequency by ploidy a second time, so a diploid MAF was half the true value and a locus fixed for the alternative allele got 0.5. Every filter on `maf` in `apply_sequence_filtering()` and `random_select_loci()` now keeps a different set of loci at the same threshold. `ov_maf` in `get_overall_summary()` and `loc_Fis` change as well (audit finding M1).
* `impute_gl(method = "random_forest")` now writes the predicted dosages. Earlier versions wrote the integer codes of the predicted factor, which shifted imputed dosages, so a predicted 0 was stored as 1 (I1).
* `apply_bioflow_modifications()` keeps every individual when the table removes only loci, and every locus when it removes only individuals. Earlier versions returned an empty margin in both cases (B1).

## Bug fixes

* `read_DArTSeq_SNP()` keeps all loci when every SNP is mapped. Earlier versions returned an object with no loci. Positions and chromosomes now come from the loci that were kept (R2, R3).

## Package

* The licence is MIT. `DESCRIPTION` declares `MIT + file LICENSE`, and the full text is in `LICENSE.md` (D1).
* First unit tests in `tests/testthat/`, one file per fixed finding, and a GitHub Actions workflow that runs them.
* Every function has roxygen documentation. Run `devtools::document()` to rebuild `man/`.
