"""
Codon change dictionaries
"""

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