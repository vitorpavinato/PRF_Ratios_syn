"""
Set of functions for processing .TSV file lines.
"""

from typing import Tuple, List


# Functions:
def process_snp_signature(line: List[str]) -> List[str]:
    """
    Function to process a line of the .TSV file.
    It extracts SNP coordinates, total counts and effect.
    """

    # Get the chromosome and position
    chrm = line[0]
    pos = line[1]

    # Get the total count
    totalcount = line[12]

    # Get the effect
    effect = line[19]

    snp_fields = [chrm, pos, effect, totalcount]

    return snp_fields


def process_alleles(line: List[str]) -> Tuple[str, List[str], List[str], List[str]]:
    """
    Function to process a line of the .TSV file.
    It extracts the allele information used for rooting SNPs:
    root, snp_alleles, refcount, altcount, refcodon, altcodon.
    """

    # Get the root
    root = line[2]

    # Get the reference and alternative
    refallele = line[3]
    altallele = line[4]
    snp_alleles = [refallele, altallele]

    # Get the allele counts
    refcount = line[10]
    altcount = line[11]
    allele_counts = [refcount, altcount]

    # Get the allele codons
    refcodon = line[17].upper()
    altcodon = line[18].upper()
    allele_codons = [refcodon, altcodon]

    return root, snp_alleles, allele_counts, allele_codons


def root_snp(snp_alleles: List[str],
             allele_counts: List[str],
             allele_codons: List[str]
             ) -> Tuple[List[str], List[str], List[str]]:

    """
    Function to root SNPs in a tsv table
    """

    # Swap reference and alternative codons if root is root_alt:
    print("Swapping reference and alternative alleles fields")

    # Swap allele fields because of the root
    snp_alleles.reverse()
    allele_counts.reverse()
    allele_codons.reverse()

    return snp_alleles, allele_counts, allele_codons


def make_codon_change_keys(allele_codons: List[str]) -> str:
    """
    Combine rooted codons into a string representing
    the direction of the codon change.
    This will be used as a key in the dictionary.
    """

    codon_change = "->".join(allele_codons)

    return codon_change


def make_mutational_change_keys(line: List[str]) -> str:
    """
    Combine rooted mutational context into a string representing
    the direction of mutational change. This will be used to match
    codon change directionwith mutational context on neutral SNPs.
    This will be used as a key in the dictionary.
    """

    # Get the mutation context from the line
    refcontext = line[13][2:5]
    altcontext = line[14][2:5]

    mutational_change = "->".join([refcontext, altcontext])

    return mutational_change
