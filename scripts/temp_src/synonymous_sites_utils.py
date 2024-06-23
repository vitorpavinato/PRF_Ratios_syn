"""
Script to prepare synonymous SFSS to analyze with PRF-ratios.
"""

from typing import List


# Define the numerical keys for the synonymous pairs
rooted_synonymous_pairs_dict = {
    1: "TTT->TTC", 2: "TTC->TTT",
    3: "CTT->CTC", 4: "CTC->CTT",
    5: "CTT->CTA", 6: "CTA->CTT",
    7: "CTT->CTG", 8: "CTG->CTT",
    9: "CTA->CTC", 10: "CTC->CTA",
    11: "CTA->CTG", 12: "CTG->CTA",
    13: "CTC->CTG", 14: "CTG->CTC",
    15: "CTG->TTG", 16: "TTG->CTG",
    17: "CTA->TTA", 18: "TTA->CTA",
    19: "TTA->TTG", 20: "TTG->TTA",
    21: "ATT->ATC", 22: "ATC->ATT",
    23: "ATT->ATA", 24: "ATA->ATT",
    25: "ATA->ATC", 26: "ATC->ATA",
    27: "GTT->GTC", 28: "GTC->GTT",
    29: "GTT->GTA", 30: "GTA->GTT",
    31: "GTT->GTG", 32: "GTG->GTT",
    33: "GTC->GTA", 34: "GTA->GTC",
    35: "GTC->GTG", 36: "GTG->GTC",
    37: "GTA->GTG", 38: "GTG->GTA",
    39: "TCT->TCC", 40: "TCC->TCT",
    41: "TCT->TCA", 42: "TCA->TCT",
    43: "TCT->TCG", 44: "TCG->TCT",
    45: "TCC->TCA", 46: "TCA->TCC",
    47: "TCC->TCG", 48: "TCG->TCC",
    49: "TCA->TCG", 50: "TCG->TCA",
    51: "AGT->AGC", 52: "AGC->AGT",
    53: "CCT->CCC", 54: "CCC->CCT",
    55: "CCT->CCA", 56: "CCA->CCT",
    57: "CCT->CCG", 58: "CCG->CCT",
    59: "CCC->CCA", 60: "CCA->CCC",
    61: "CCC->CCG", 62: "CCG->CCC",
    63: "CCA->CCG", 64: "CCG->CCA",
    65: "ACT->ACC", 66: "ACC->ACT",
    67: "ACT->ACA", 68: "ACA->ACT",
    69: "ACT->ACG", 70: "ACG->ACT",
    71: "ACC->ACA", 72: "ACA->ACC",
    73: "ACC->ACG", 74: "ACG->ACC",
    75: "ACA->ACG", 76: "ACG->ACA",
    77: "GCT->GCC", 78: "GCC->GCT",
    79: "GCT->GCA", 80: "GCA->GCT",
    81: "GCT->GCG", 82: "GCG->GCT",
    83: "GCC->GCA", 84: "GCA->GCC",
    85: "GCC->GCG", 86: "GCG->GCC",
    87: "GCA->GCG", 88: "GCG->GCA",
    89: "TAT->TAC", 90: "TAC->TAT",
    91: "CAT->CAC", 92: "CAC->CAT",
    93: "CAA->CAG", 94: "CAG->CAA",
    95: "AAT->AAC", 96: "AAC->AAT",
    97: "AAA->AAG", 98: "AAG->AAA",
    99: "GAT->GAC", 100: "GAC->GAT",
    101: "GAA->GAG", 102: "GAG->GAA",
    103: "TGT->TGC", 104: "TGC->TGT",
    105: "CGT->CGC", 106: "CGC->CGT",
    107: "CGT->CGA", 108: "CGA->CGT",
    109: "CGT->CGG", 110: "CGG->CGT",
    111: "CGC->CGA", 112: "CGA->CGC",
    113: "CGC->CGG", 114: "CGG->CGC",
    115: "CGA->CGG", 116: "CGG->CGA",
    117: "CGA->AGA", 118: "AGA->CGA",
    119: "AGA->AGG", 120: "AGG->AGA",
    121: "CGG->AGG", 122: "AGG->CGG",
    123: "GGT->GGC", 124: "GGC->GGT",
    125: "GGT->GGA", 126: "GGA->GGT",
    127: "GGT->GGG", 128: "GGG->GGT",
    129: "GGC->GGA", 130: "GGA->GGC",
    131: "GGC->GGG", 132: "GGG->GGC",
    133: "GGA->GGG", 134: "GGG->GGA",
}

# Convert the numerical keys to string keys
rooted_synonymous_pairs_dict = {value: {} for value in rooted_synonymous_pairs_dict.values()}


def create_synonymous_codon_pair_sfs_dict(
    tsv_lines: List[List[str]],
    dict_to_fill: dict[int, str]
) -> dict[int, dict[int, list[int]]]:
    """
    Takes a list of lines from a .TSV file and creates a dictionary
    of sample sizes (total counts) for each synonymous pairs numeric key.
    """
    dict_to_fill = rooted_synonymous_pairs_dict

    # Make a copy of the dictionary
    dict_filled = dict_to_fill

    # Fill up codon pair dictionary with data
    for _, line in enumerate(tsv_lines_chr2L):
        # print(i, line)
        total_count = int(line[3])
        derived_count = int(line[5])
        codon_pair = line[11]

        # Fill up the rooted codon pair dictionary with data
        if total_count not in dict_filled[codon_pair]:
            dict_filled[codon_pair][total_count] = []

        dict_filled[codon_pair][total_count].append(derived_count)

    # Sort the list of derived counts for each codon pair total count:
    for codon_pair, total_counts in dict_filled.items():
        for total_count, derived_counts in total_counts.items():
            dict_filled[codon_pair][total_count] = sorted(derived_counts)

    return dict_filled


test_dict = create_synonymous_codon_pair_dict(tsv_lines_chr2L, rooted_synonymous_pairs_dict)

test_dict['CGT->CGA']


# Working on this part



codon_pairs_sfss_dict = {}


for codon_pair, total_counts in test_dict.items():
    for total_counts, derived_counts in total_counts.items():
        print(codon_pair, total_counts, derived_counts)
        if codon_pair not in codon_pairs_sfss_dict:
            codon_pairs_sfss_dict[codon_pair] = {}
            if total_counts not in codon_pairs_sfss_dict[codon_pair]:
                codon_pairs_sfss_dict[codon_pair][total_counts] = [0] * (total_counts + 1)
                for derived_count in derived_counts:
                    codon_pairs_sfss_dict[codon_pair][total_counts][derived_count] += 1
        

for codon_pair, total_counts in codon_pairs_sfss_dict.items():
    for total_counts, derived_counts in total_counts.items():
        print(codon_pair, total_counts, derived_counts)
        

test_dict['CGT->CGA']
codon_pairs_sfss_dict['CGT->CGA']


for codon_pair, total_counts in snp_codon_pairs_dict.items():
    for total_count, derived_counts in total_counts.items():
        if codon_pair not in codon_pairs_sfss_dict:
            codon_pairs_sfss_dict[codon_pair] = {}
            codon_pairs_sfss_dict[codon_pair][total_count] = [0] * (total_count + 1)
            for derived_count in derived_counts:
                codon_pairs_sfss_dict[codon_pair][total_count][derived_count] += 1


