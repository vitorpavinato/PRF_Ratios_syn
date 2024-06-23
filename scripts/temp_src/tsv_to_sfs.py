"""
Main function to convert .TSV files to .SFS files
"""

from implementations import *


def main() -> None:

    """
    Main function
    """

    # Call the function
    tsv_lines_chr2L = process_functional_snps_tsv(
        inputfile="data/dgrp2/NC_Chr2L_tables.tsv",
        outputfile="data/dgrp2/NC_Chr2L_rooted.tsv",
        functional_effect="SYNONYMOUS_CODING"
    )

    for line in tsv_lines_chr2L:
        print(line)

    # Make a matrix for each codon_pair
    # matrix should be total_count x derived_count
    # for each codon_pair
    # Create the matrix for each codon_pair
    # for codon_pair, total_counts in snp_codon_pairs_dict.items():
    #     matrix = []
    #     for total_count, derived_counts in total_counts.items():
    #         row = [derived_count for derived_count in derived_counts]
    #         matrix.append(row)
    #     snp_codon_pairs_dict[codon_pair] = matrix

    # # Print the dictionary
    # for codon_pair, matrix in snp_codon_pairs_dict.items():
    #     print(f"Codon Pair: {codon_pair}")
    #     for row in matrix:
    #         print(row)
    #     print()




if __name__ == "__main__":
    main()
