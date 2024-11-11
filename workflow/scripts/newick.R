suppressPackageStartupMessages(library(ape))

build_newick <- function(distance_matrix_file, newick_file) {
  distance_matrix <- as.matrix(read.table(distance_matrix_file, sep="\t", header=TRUE, row.names=1))

  hc <- hclust(as.dist(distance_matrix), method="complete")

  my_tree <- as.phylo(hc)
  write.tree(phy=my_tree, file=newick_file)
}


build_newick(snakemake@input[[1]], snakemake@output[[1]])
