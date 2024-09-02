# Script to extract rooting information from a table info column.
# The roots were inferred with a maximum likelihood approach
# from MSA and a given phylogenetic tree.
# The roots were added to a vcf file.
# I extracted the column containing the root information
# and here I am cleaning this up to make it available as
# a .tsv file that can be merged with other tables from other
# populations. In a way, I am transfering the roots from
# from one population to another.

extract_root_info <- function(
    table_name,
    output_name) {

  table_ <- read.table(table_name, sep = "\t", header = TRUE)
  table_["roots"] <- sapply(
                            strsplit(table_$id, "_SNP_"),
                            function(x) x[2])
  table_$roots[is.na(table_$roots)] <- "unknown"

  # Write the table
  write.table(table_, output_name, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Input paths
roots_path <- "data/roots/"
output_path <- "results/extracted_roots/"

# Create output path
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

chrms = c("Chr2L", "Chr2R", "Chr3L", "Chr3R")
for (chr in chrms) {
  roots_name = paste(roots_path, "dgrp2dm6_", chr, "_rooted_snps.tsv", sep = "")
  output_name = paste(output_path, "dgrp2dm6_", chr, "_extracted_roots.tsv", sep = "")
  extract_root_info(roots_name, output_name)
}