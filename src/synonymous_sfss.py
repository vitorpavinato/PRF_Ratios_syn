"""
Script to prepare synonymous SFSS to analyze with PRF-ratios.
"""

synonymous_pairs_dict = {
    1: "TTT->TTC",
    2: "TTC->TTT",
    3: "CTT->CTC",
    4: "CTC->CTT",
    5: "CTT->CTA",
    6: "CTA->CTT",
    7: "CTT->CTG",
    8: "CTG->CTT",
    9: "CTA->CTC",
    10: "CTC->CTA",
    11: "CTA->CTG",
    12: "CTG->CTA",
    13: "CTC->CTG",
    14: "CTG->CTC",
    15: "CTG->TTG",
    16: "TTG->CTG",
    17: "CTA->TTA",
    18: "TTA->CTA",
    19: "TTA->TTG",
    20: "TTG->TTA",
    21: "ATT->ATC",
    22: "ATC->ATT",
    23: "ATT->ATA",
    24: "ATA->ATT",
    25: "ATA->ATC",
    26: "ATC->ATA",
    27: "GTT->GTC",
    28: "GTC->GTT",
    29: "GTT->GTA",
    30: "GTA->GTT",
    31: "GTT->GTG",
    32: "GTG->GTT",
    33: "GTC->GTA",
    34: "GTA->GTC",
    35: "GTC->GTG",
    36: "GTG->GTC",
    37: "GTA->GTG",
    38: "GTG->GTA",
    39: "TCT->TCC",
    40: "TCC->TCT",
    41: "TCT->TCA",
    42: "TCA->TCT",
    43: "TCT->TCG",
    44: "TCG->TCT",
    45: "TCC->TCA",
    46: "TCA->TCC",
    47: "TCC->TCG",
    48: "TCG->TCC",
    49: "TCA->TCG",
    50: "TCG->TCA",
    51: "AGT->AGC",
    52: "AGC->AGT",
    53: "CCT->CCC",
    54: "CCC->CCT",
    55: "CCT->CCA",
    56: "CCA->CCT",
    57: "CCT->CCG",
    58: "CCG->CCT",
    59: "CCC->CCA",
    60: "CCA->CCC",
    61: "CCC->CCG",
    62: "CCG->CCC",
    63: "CCA->CCG",
    64: "CCG->CCA",
    65: "ACT->ACC",
    66: "ACC->ACT",
    67: "ACT->ACA",
    68: "ACA->ACT",
    69: "ACT->ACG",
    70: "ACG->ACT",
    71: "ACC->ACA",
    72: "ACA->ACC",
    73: "ACC->ACG",
    74: "ACG->ACC",
    75: "ACA->ACG",
    76: "ACG->ACA",
    77: "GCT->GCC",
    78: "GCC->GCT",
    79: "GCT->GCA",
    80: "GCA->GCT",
    81: "GCT->GCG",
    82: "GCG->GCT",
    83: "GCC->GCA",
    84: "GCA->GCC",
    85: "GCC->GCG",
    86: "GCG->GCC",
    87: "GCA->GCG",
    88: "GCG->GCA",
    89: "TAT->TAC",
    90: "TAC->TAT",
    CAT->CAC
    CAC->CAT
    CAA->CAG
    CAG->CAA
    AAT->AAC
    AAC->AAT
    AAA->AAG
    AAG->AAA
    GAT->GAC
    GAC->GAT
    GAA->GAG
    GAG->GAA
    TGT->TGC
    TGC->TGT
    CGT->CGC
    CGC->CGT
    CGT->CGA
    CGA->CGT
    CGT->CGG
    CGG->CGT
    CGC->CGA
    CGA->CGC
    CGC->CGG
    CGG->CGC
    CGA->CGG
    CGG->CGA
    CGA->AGA
    AGA->CGA
    AGA->AGG
    AGG->AGA
    CGG->AGG
    AGG->CGG
    GGT->GGC
    GGC->GGT
    GGT->GGA
    GGA->GGT
    GGT->GGG
    GGG->GGT
    GGC->GGA
    GGA->GGC
    GGC->GGG
    GGG->GGC
    GGA->GGG
    GGG->GGA
}