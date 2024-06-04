"""
Implementations of the functions to obtain unfolded SFS from tsv files.
For the moment, only make sense for synonymous and neutral mutations.
"""

from typing import List
from processing_tsv_utils import *


def make_header() -> List[str]:
    """
    Function to make the header
    """

    header = ["chrom", "pos", "aainfo",
              "refallele", "altallele",
              "refcount", "altcount", "totalcount",
              "refcodon", "altcodon", "effect"]

    header = "\t".join(str(item) for item in header)

    return header


def write_tsv_file(lines: List[List[str]], header: List[str], outputfile: str) -> None:
    """
    Write the processed .TSV to another .TSV file
    """

    # Prompt message
    print("Exporting the processed TSV to a .TSV file")
    with open(outputfile, "w", encoding="utf-8") as output_file:
        output_file.write(header + "\n")
        for line in lines:
            output_file.write(("\t".join(str(item) for item in line)) + "\n")
    print("tsv file exported to: " + outputfile)


def synonymous_tsv(inputfile: str,
                   outputfile: str = None,
                   filter_by_effect: str = "SYNONYMOUS_CODING"
                   ) -> List[List[str]]:

    """
    This function takes a tsv file and returns a tsv file
    with synonymous mutations
    """

    tsv_lines = []

    # Open the input file
    with open(inputfile, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue

            # Strip the line
            line_stripped = line.strip()

            # Split the line by tabs
            line_split = line_stripped.split("\t")

            # Process snp coordinates, total count and effect
            snp_fields = processing_snp_fields(line_split)

            if snp_fields[2] == filter_by_effect:

                # Process information used for rooting SNPs
                root, snp_alleles, allele_counts, allele_codons = processing_fields_for_rooting(line_split)

                # Swap reference and alternative codons if root is root_alt:
                if root == "root_alt":
                    snp_alleles, allele_counts, allele_codons = root_snps(snp_alleles, allele_counts, allele_codons)

                # Make codon change keys
                codon_change = make_codon_change_keys(allele_codons)

                # Assemble the new line
                new_line = snp_fields + allele_counts + [root] + snp_alleles + allele_codons + [codon_change]

                # Write the new line to the output file
                tsv_lines.append(new_line)

    # Output file if specified
    if outputfile is not None:
        header = make_header()
        write_tsv_file(tsv_lines, header, outputfile)

    return tsv_lines


def main() -> None:

    """
    Main function
    """

    # Call the function
    tsv_lines = synonymous_tsv(
                                inputfile="data/dgrp2/NC_Chr2L_tables.tsv",
                                outputfile="data/dgrp2/NC_Chr2L_rooted_2.tsv",
                                filter_by_effect="SYNONYMOUS_CODING"
    )

    # Make dictionary:
    snp_codon_pairs_dict = {}

    for i, line in enumerate(tsv_lines):
        print(i, line)
        total_count = line[3]
        derived_count = line[5]
        codon_pair = line[11]

        if codon_pair not in snp_codon_pairs_dict:
            snp_codon_pairs_dict[codon_pair] = {}

        if total_count not in snp_codon_pairs_dict[codon_pair]:
            snp_codon_pairs_dict[codon_pair][total_count] = []

        snp_codon_pairs_dict[codon_pair][total_count].append(derived_count)

    # Sort the total_count for each codon_pair
    for codon_pair, total_counts in snp_codon_pairs_dict.items():
        for total_count, derived_counts in total_counts.items():
            snp_codon_pairs_dict[codon_pair][total_count] = sorted(derived_counts)
            

    # Print the dictionary
    print(snp_codon_pairs_dict['ATT->ATA'])

    # Print the maximum and minimum values for each codon_pair
    for codon_pair, total_counts in snp_codon_pairs_dict.items():
        # print(codon_pair, [int(total_count) for total_count in total_counts.keys()]) 
        if total_counts:
            max_total_count = max(int(total_count) for total_count in total_counts.keys())
            min_total_count = min(int(total_count) for total_count in total_counts.keys())
            print(f"Codon Pair: {codon_pair}, Max Total Count: {max_total_count}, Min Total Count: {min_total_count}")
        else:
            print(f"Codon Pair: {codon_pair}, Max Total Count: None, Min Total Count: None")

    
    # Make a matrix for each codon_pair
    # matrix should be total_count x derived_count
    # for each codon_pair
    # Create the matrix for each codon_pair
    for codon_pair, total_counts in snp_codon_pairs_dict.items():
        matrix = []
        for total_count, derived_counts in total_counts.items():
            row = [derived_count for derived_count in derived_counts]
            matrix.append(row)
        snp_codon_pairs_dict[codon_pair] = matrix

    # Print the dictionary
    for codon_pair, matrix in snp_codon_pairs_dict.items():
        print(f"Codon Pair: {codon_pair}")
        for row in matrix:
            print(row)
        print()




if __name__ == "__main__":
    main()
