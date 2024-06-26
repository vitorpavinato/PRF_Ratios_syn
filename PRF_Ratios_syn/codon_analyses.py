"""
Module for analysing codon changes.
"""

import json
import pickle
# from itertools import product  # For generating all possible synonymous codon changes
# from collections import defaultdict
# from typing import List
from pandas import DataFrame


# # Dictionary and list
# genetic_code = {
#     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
#     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
#     'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#     'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#     'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#     'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#     'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#     'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
# }


# # Fuction to generate synonymous pairs from the genetic code
# def generate_synonymous_pairs(genetic_code: dict) -> List[str]:
#     """
#     Function to generate synonymous pairs from the genetic code.
#     """
#     synonymous_pairs = []
#     for codon1, codon2 in product(genetic_code.keys(), repeat=2):
#         if codon1 != codon2 and genetic_code[codon1] == genetic_code[codon2]:
#             synonymous_pairs.append(f"{codon1}->{codon2}")
#     return synonymous_pairs


# # Instantiate the synonymous pairs from the genetic code
# # The list has 180 synonymous changes including stop codons
# # and changes involving more than one nucleotide
# synonymous_pairs = generate_synonymous_pairs(genetic_code)


# The list include 134 synonymous changes involving only one nucleotide
synonymous_1nt_pairs = ["TTT->TTC", "TTC->TTT", "CTT->CTC", "CTC->CTT", "CTT->CTA", "CTA->CTT",
                        "CTT->CTG", "CTG->CTT", "CTA->CTC", "CTC->CTA", "CTA->CTG", "CTG->CTA",
                        "CTC->CTG", "CTG->CTC", "CTG->TTG", "TTG->CTG", "CTA->TTA", "TTA->CTA",
                        "TTA->TTG", "TTG->TTA", "ATT->ATC", "ATC->ATT", "ATT->ATA", "ATA->ATT",
                        "ATA->ATC", "ATC->ATA", "GTT->GTC", "GTC->GTT", "GTT->GTA", "GTA->GTT",
                        "GTT->GTG", "GTG->GTT", "GTC->GTA", "GTA->GTC", "GTC->GTG", "GTG->GTC",
                        "GTA->GTG", "GTG->GTA", "TCT->TCC", "TCC->TCT", "TCT->TCA", "TCA->TCT",
                        "TCT->TCG", "TCG->TCT", "TCC->TCA", "TCA->TCC", "TCC->TCG", "TCG->TCC",
                        "TCA->TCG", "TCG->TCA", "AGT->AGC", "AGC->AGT", "CCT->CCC", "CCC->CCT",
                        "CCT->CCA", "CCA->CCT", "CCT->CCG", "CCG->CCT", "CCC->CCA", "CCA->CCC",
                        "CCC->CCG", "CCG->CCC", "CCA->CCG", "CCG->CCA", "ACT->ACC", "ACC->ACT",
                        "ACT->ACA", "ACA->ACT", "ACT->ACG", "ACG->ACT", "ACC->ACA", "ACA->ACC",
                        "ACC->ACG", "ACG->ACC", "ACA->ACG", "ACG->ACA", "GCT->GCC", "GCC->GCT",
                        "GCT->GCA", "GCA->GCT", "GCT->GCG", "GCG->GCT", "GCC->GCA", "GCA->GCC",
                        "GCC->GCG", "GCG->GCC", "GCA->GCG", "GCG->GCA", "TAT->TAC", "TAC->TAT",
                        "CAT->CAC", "CAC->CAT", "CAA->CAG", "CAG->CAA", "AAT->AAC", "AAC->AAT",
                        "AAA->AAG", "AAG->AAA", "GAT->GAC", "GAC->GAT", "GAA->GAG", "GAG->GAA",
                        "TGT->TGC", "TGC->TGT", "CGT->CGC", "CGC->CGT", "CGT->CGA", "CGA->CGT",
                        "CGT->CGG", "CGG->CGT", "CGC->CGA", "CGA->CGC", "CGC->CGG", "CGG->CGC",
                        "CGA->CGG", "CGG->CGA", "CGA->AGA", "AGA->CGA", "AGA->AGG", "AGG->AGA",
                        "CGG->AGG", "AGG->CGG", "GGT->GGC", "GGC->GGT", "GGT->GGA", "GGA->GGT",
                        "GGT->GGG", "GGG->GGT", "GGC->GGA", "GGA->GGC", "GGC->GGG", "GGG->GGC",
                        "GGA->GGG", "GGG->GGA"]


# # This code counts the number of nucleotide differences between two codons
# # Help to check if all changes in list 2 involve only one nucleotide
# def count_nucleotide_differences(codon1, codon2):
#     return sum(n1 != n2 for n1, n2 in zip(codon1, codon2))
#
# all_single_nucleotide = all(count_nucleotide_differences(change.split('->')[0], change.split('->')[1]) == 1 for change in synonymous_pairs)
#
# print(f"All changes in list 2 involve only one nucleotide: {all_single_nucleotide}")


def create_codon_change_sfs_dict(df: DataFrame, use_filter: bool) -> dict:
    """
    Function to create codon change dictionary of total counts
    and derived counts SFS. It filter out SNPs with exon-intron
    junctions annotation.
    """
    # Generate all possible synonymous codon changes
    synonymous_changes = synonymous_1nt_pairs

    # Create the nested dictionary structure
    # codon_dict = {change: defaultdict(lambda: defaultdict(int)) for change in synonymous_changes}
    codon_dict = {}
    for change in synonymous_changes:
        codon_dict[change] = {}

    if use_filter:
        df = df[df['custom_annotation'] != 'eij']

    # Fill the dictionary with data from the DataFrame
    for _, row in df.iterrows():
        codon_change = row['codon_change']

        # if codon_change in synonymous_changes and extra_annotation != 'eij':
        if codon_change in synonymous_changes:
            total_count = row['totalcount']
            alt_count = row['altcount']

            if total_count not in codon_dict[codon_change]:
                codon_dict[codon_change][total_count] = [0] * (total_count + 1)

            codon_dict[codon_change][total_count][alt_count] += 1

    return codon_dict


# Define the functions to save and load data
def save_data(data: dict, pickle_file: str, json_file: str):
    """
    Function to save data to pickle and json files.
    """
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    with open(json_file, 'w') as f:
        json.dump({k: dict(v) for k, v in data.items()}, f, indent=2)
