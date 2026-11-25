# Version 1.0.3, fixes from a code audit and an external review, tests and CI

This PR brings cgiarGenomics from 1.0.2 to 1.0.3. It fixes the defects found by two reviews in November 2026, adds unit tests and a CI workflow, and documents the package for people and AI assistants. Function names, exported functions and existing arguments are unchanged.

## Why

Version 1.0.2 gives wrong results on several paths that Bioflow exposes, without raising an error. An external review ran 1.0.2 on the diploid hapmap example. At a 5% MAF cutoff it excluded 63 eligible markers of 3,553 and kept one that should have failed. Random Forest imputation wrote a dosage of 3 in a diploid. A filter that no locus passed was skipped, and all loci were kept.

## What changes results

- **MAF.** `get_maf()` divided the allele frequency by ploidy twice. Stored `maf`, `ov_maf`, `loc_Fis` and every MAF filter change.
- **Heterozygosity.** Now the share of heterozygotes among called genotypes. Missing calls no longer lower `loc_het`, `ind_het` and `ov_het`.
- **Random Forest imputation.** Writes the predicted dosages instead of factor codes, and is reproducible from its seed.
- **Empty filters.** `apply_sequence_filtering()` stops with an error of class `cgiarGenomics_empty_filter` when nothing passes a step.
- **Replicate merging.** `merge_duplicate_inds()` groups by `designation_id`, names individuals by designation, keeps positions and alleles, and works when no designation has replicates.
- **Duplicated marker IDs.** Renamed to `chrom_pos` when read.
- **Imputation accuracy.** `get_accuracy()` takes an optional `masked_gl` and then scores only the masked cells. `set_fraction_na()` masks called genotypes only and keeps names, positions and alleles.
- **Beagle.** Markers are sorted by position before Beagle runs, and results are matched back by name.
- **VCF export.** Missing polyploid genotypes are written at their ploidy, `./././.` for a tetraploid.

## Bug fixes

- `read_DArTSeq_SNP()` no longer drops every locus when all SNPs are mapped.
- `read_hapmap()` reads calls with `/` or `|` and decodes TASSEL IUPAC calls.
- `read_vcf()` reads mixed phased and unphased calls, and files with one variant.
- `impute_gl()` accepts a panel with no missing calls and stops when a method returns an impossible dosage.
- `apply_bioflow_modifications()` keeps the other margin when a table removes only loci or only individuals.
- `random_select_loci()` works when two chromosomes tie for the most loci.

## Security

`impute_beagle()` pasted the `memory` setting into a shell command, and Bioflow passes that value from a text field in the QA form. Settings are now checked before any process starts, and Beagle runs through `system2()` with quoted arguments. We still recommend that shared Bioflow instances bound memory, threads and run time on the server.

## What Bioflow needs to handle

- Catch `cgiarGenomics_empty_filter` in the QA preview and show the user which step removed everything.
- `write_vcf()` now defaults to `na_rep = NULL`. Diploid output is unchanged.
- `merge_duplicate_inds()` returns individuals named by `designation_id`, singles first, then merged designations in alphabetical order.
- The package now requires R 4.1, which the code already needed.
- Past analyses that used a positive MAF cutoff, thinning, Fis, heterozygosity, Random Forest or replicate merging will differ. A rerun of a sample of them would show how much.

## Package and documentation

- Licence set to MIT in `DESCRIPTION`, with the full text in `LICENSE.md`.
- `Authors@R` lists Sergio Andres Cruz and Lorena Guimarães Batista. New `CITATION.cff`, so GitHub shows "Cite this repository".
- The bundled Beagle 5.5 JAR is under GPL 3 or later, not Apache 2.0 as the old note said. `inst/extdata/README.md` gives its source, version and checksum.
- roxygen now writes `NAMESPACE` and `man/`. Every exported function has a help page with a runnable example.
- Removed `read_genolcsv()` and three plot helpers. None was exported and none could run.
- New `NEWS.md`, `REFERENCE.md`, `dev/audit.md`, a rewritten `README.md`, and `AGENTS.md`, `SKILL.md` and `llms.txt` so that AI assistants pick the right function and avoid known pitfalls.

## Checks

- 116 test expectations in `tests/testthat/` pass, none skipped, with vcfR, randomForest, dartR.base and Java installed.
- `R CMD check --no-manual` ends with Status OK on R 4.6.2, Ubuntu.
- The `tests` workflow in `.github/workflows/` runs both on every pull request and fails on warnings.

## Still open

`dev/audit.md` lists every finding with its status. The main open items are these.

1. Provenance and reuse terms of the example files in `tests/`, which hold provider-style sample and project identifiers.
2. Whether to keep bundling the Beagle JAR, and shipping the GPL text next to it.
3. `loc_Fis` assumes a diploid at every ploidy.
4. `apply_bioflow_modifications()` does not replay imputation rows. It now warns when they are present.
5. No Windows CI job, no large-panel benchmark and no end-to-end run through Bioflow yet.
6. ORCIDs and affiliations for `CITATION.cff`, a DOI and a tagged release.

Happy to split this into smaller PRs if that is easier to review.
