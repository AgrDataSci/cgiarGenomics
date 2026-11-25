---
name: cgiargenomics
description: Read, check, filter, impute and export SNP marker data (hapmap, VCF, DArTSeq, DArTag) in R with cgiarGenomics, diploid or polyploid. Use for genomic data steps in Bioflow or for any question about this package.
---

# cgiarGenomics

cgiarGenomics turns marker files from genotyping providers into a genlight object of alternative allele dosages, then measures quality, filters, imputes and writes VCF. It is the genomic data engine of Bioflow, the CGIAR breeding analytics pipeline. Load the files below only when the question needs them.

| File | Open it when |
|---|---|
| [`AGENTS.md`](AGENTS.md) | Always, before writing code. Task to function table, known issues and pitfalls. |
| [`REFERENCE.md`](REFERENCE.md) | You need the arguments or return value of a function. |
| [`dev/audit.md`](dev/audit.md) | A result looks wrong, or the user asks about reliability or wants to fix the package. |
| [`README.md`](README.md) | You need the supported formats with examples, or installation. |

## Standard workflow

```r
library("cgiarGenomics")

check_dependencies()                       # optional packages for VCF, DArTSeq, RF

gl = read_vcf("file.vcf.gz", ploidity = 2) # or read_hapmap(), read_DArTag_count_dosage()
get_overall_summary(gl)

steps = list(
  list(param = "loc_miss", operator = "<=", threshold = 0.2),
  list(param = "ind_miss", operator = "<=", threshold = 0.3)
)
filtered = apply_sequence_filtering(gl, steps)
filtered$filt_log                          # check that no step has an empty index

imputed = impute_gl(filtered$gl, ploidity = 2, method = "frequency")
write_vcf(imputed$gl, "clean.vcf.gz")
```

## Rules

1. Run code before reporting results, and report the numbers of loci and individuals kept at each step.
2. Check `packageVersion("cgiarGenomics")`. Version 1.0.2 and earlier give wrong `maf` values and wrong dosages from `method = "random_forest"`. Ask the user to update before using either.
3. Tell the user when a known issue from `AGENTS.md` applies to their analysis.
4. Keep `ploidity` equal to the ploidy of the data in every call.
5. Never invent function names or arguments. Check `REFERENCE.md` or the R help.

## Cite

Cite the package with `citation("cgiarGenomics")`, and Bioflow when the work ran through it.
