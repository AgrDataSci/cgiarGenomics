# Small genlight objects shared by the tests. Sample and marker names are
# invented.

make_gl = function(mt, ploidy = 2, chromosome = NULL, position = NULL) {
  if (is.null(rownames(mt))) {
    rownames(mt) = paste0("ind", seq_len(nrow(mt)))
  }
  if (is.null(colnames(mt))) {
    colnames(mt) = paste0("snp", seq_len(ncol(mt)))
  }
  if (is.null(chromosome)) {
    chromosome = rep("1", ncol(mt))
  }
  if (is.null(position)) {
    position = seq_len(ncol(mt)) * 100L
  }
  gl = methods::new("genlight", mt, ploidy = ploidy,
                    chromosome = chromosome, position = position)
  adegenet::alleles(gl) = rep("A/G", ncol(mt))
  gl
}

# Write a small hapmap file and return its path. `calls` is a matrix of
# genotype calls with loci in rows and samples in columns.
write_hapmap = function(calls, alleles, chrom = "1", pos = seq_len(nrow(calls)) * 100L) {
  colnames(calls) = paste0("ind", seq_len(ncol(calls)))
  meta = data.frame(
    `rs#` = paste0("snp", seq_len(nrow(calls))), alleles = alleles,
    chrom = chrom, pos = pos, strand = "+", `assembly#` = "v1", center = "x",
    protLSID = "x", assayLSID = "x", panelLSID = "x", QCcode = "x",
    check.names = FALSE
  )
  path = tempfile(fileext = ".hmp.txt")
  utils::write.table(cbind(meta, calls), path, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  path
}
