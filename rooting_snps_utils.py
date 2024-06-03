"""
Set of functions for rooting SNPs
"""

from pathlib import Path


def root_snps(inputfile: str,
              outputfile: str,
              filter_by_effect: str = None
              ) -> None:

    """
    Function to return a tsv file with the rooted SNPs.
    It can also return only SNPs that have the specified effect.
    """

    header = ["chrom", "pos", "aainfo",
              "ref", "alt",
              "refcount", "altcount", "totalcount",
              "refcodon", "altcodon", "effect"]

    with open(outputfile, "w", encoding="utf-8") as output_file:
        output_file.write("\t".join(header) + "\n")

    # Open the input and output files
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "a", encoding="utf-8") as output_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            line_stripped = line.strip()
            line_split = line_stripped.split("\t")

            # Get the chromosome and position
            chrm = line_split[0]
            pos = line_split[1]
            
            # Get the root
            root = line_split[2]

            # Get the reference and alternative
            refallele = line_split[3]
            altallele = line_split[4]

            # Get the allele counts
            refcount = line_split[10]
            altcount = line_split[11]
            totalcount = line_split[12]

            # Get the allele codons
            refcodon = line_split[17].upper()
            altcodon = line_split[18].upper()

            # Get the effect
            effect = line_split[19]

            # Swap reference and alternative codons if root is root_alt:
            if root == "root_alt":

                print("Swapping reference and alternative alleles fields")

                # Swap allele fields because of the root
                refallele, altallele = altallele, refallele
                refcount, altcount = altcount, refcount
                refcodon, altcodon = altcodon, refcodon

            new_line = [chrm, pos, root,
                        refallele, altallele,
                        refcount, altcount, totalcount,
                        refcodon, altcodon,
                        effect]

            # Filter lines by effect
            if filter_by_effect is not None:
                if line_split[19] == filter_by_effect:

                    # Write the new line to the output file
                    output_file.write("\t".join(new_line) + "\n")

            else:

                # Write the new line to the output file
                output_file.write("\t".join(new_line) + "\n")


def main() -> None:

    """
    Main function
    """

    # Call the function
    root_snps(
        inputfile=Path("data/dgrp2/NC_Chr2L_tables.tsv"),
        outputfile=Path("data/dgrp2/NC_Chr2L_rooted.tsv"),
        filter_by_effect="SYNONYMOUS_CODING"
    )


if __name__ == "__main__":
    main()
