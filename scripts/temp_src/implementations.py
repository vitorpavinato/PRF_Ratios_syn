"""
Rail road implementation of process_tsv functions.
Path one is for functional SNPs: SYNONYMOUS_CODING or NON_SYNONYMOUS_CODING
Path two is for non-functional SNPs: INTERGENIC, UPSTREAM, DOWNSTREAM etc.
"""

from typing import List
from process_tsv_utils import *


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
    Write the processed .TSV file.
    """

    # Prompt message
    print("Exporting the processed .TSV file lines to: " + outputfile)
    with open(outputfile, "w", encoding="utf-8") as output_file:
        output_file.write(header + "\n")
        for line in lines:
            output_file.write(("\t".join(str(item) for item in line)) + "\n")


def process_functional_snps_tsv(inputfile: str,
                                outputfile: str = None,
                                functional_effect: str = "SYNONYMOUS_CODING"
                                ) -> List[List[str]]:

    """
    This function takes a .TSV file and returns another .TSV file.
    It returns a subset of the input .TSV for functional SNPs
    with only specified mutations in functional_effect:
    SYNONYMOUS_CODING or NON_SYNONYMOUS_CODING.
    """

    tsv_lines = []

    # Open the input file
    with open(inputfile, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue

            # Strip andn split the line
            line_split = line.strip().split("\t")

            # # Split the line by tabs
            # line_split = line_stripped.split("\t")

            # Process SNP signature: chromosome, pos, total count and effect
            snp_fields = process_snp_signature(line_split)

            # Retain only functional SNPs specified by functional_effect
            if snp_fields[2] == functional_effect:

                # Process information used for rooting SNPs
                root, snp_alleles, allele_counts, allele_codons = process_alleles(line_split)

                # Swap reference and alternative codons if root is root_alt:
                if root == "root_alt":
                    snp_alleles, allele_counts, allele_codons = root_snp(snp_alleles, allele_counts, allele_codons)

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


def process_nonfunctional_snps_tsv(inputfile: str,
                                   outputfile: str = None,
                                   functional_effect: str = "INTERGENIC"
                                   ) -> List[List[str]]:

    """
    This function takes a .TSV file and returns another .TSV file.
    It returns a subset of the input .TSV for functional SNPs
    with only specified mutations in functional_effect:
    INTERGENIC, UPSTREAM, DOWNSTREAM etc.
    """

    tsv_lines = []

    # Open the input file
    with open(inputfile, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue

            # Strip andn split the line
            line_split = line.strip().split("\t")

            # # Split the line by tabs
            # line_split = line_stripped.split("\t")

            # Process SNP signature: chromosome, pos, total count and effect
            snp_fields = process_snp_signature(line_split)

            # Retain only functional SNPs specified by functional_effect
            if snp_fields[2] == functional_effect:

                # Process information used for rooting SNPs
                root, snp_alleles, allele_counts, allele_codons = process_alleles(line_split)

                # Swap reference and alternative codons if root is root_alt:
                if root == "root_alt":
                    snp_alleles, allele_counts, allele_codons = root_snp(snp_alleles, allele_counts, allele_codons)

                # Make mutational context keys (NOT COMPLETED YET)
                mutational_change = make_mutational_change_keys(allele_codons)

                # Assemble the new line
                new_line = snp_fields + allele_counts + [root] + snp_alleles + allele_codons + [mutational_change]

                # Write the new line to the output file
                tsv_lines.append(new_line)

    # Output file if specified
    if outputfile is not None:
        header = make_header()
        write_tsv_file(tsv_lines, header, outputfile)

    return tsv_lines
