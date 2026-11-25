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
  adegenet::alleles(gl) = rep("a/g", ncol(mt))
  gl
}
