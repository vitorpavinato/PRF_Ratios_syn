"""
Module for analysing codon changes.
"""

import json
import pickle
from itertools import product
from collections import defaultdict
from typing import List
from pandas import DataFrame
import numpy as np
import pandas as pd


# Define the genetic code
genetic_code = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


# Define the synonymous pairs
def generate_synonymous_pairs(genetic_code: dict) -> List[str]:
    """
    Function to generate synonymous pairs from the genetic code.
    """
    synonymous_pairs = []
    for codon1, codon2 in product(genetic_code.keys(), repeat=2):
        if codon1 != codon2 and genetic_code[codon1] == genetic_code[codon2]:
            synonymous_pairs.append(f"{codon1}->{codon2}")
    return synonymous_pairs


def create_codon_change_dict(df: DataFrame) -> dict:
    # Generate all possible synonymous codon changes
    synonymous_changes = generate_synonymous_pairs()
    
    # Create the nested dictionary structure
    codon_dict = {change: defaultdict(lambda: defaultdict(int)) for change in synonymous_changes}
    phylop_dict = {change: defaultdict(lambda: defaultdict(int)) for change in synonymous_changes}
    phastcons_dict = {change: defaultdict(lambda: defaultdict(int)) for change in synonymous_changes}

    
    # Fill the dictionary with data from the DataFrame
    for _, row in df.iterrows():
        codon_change = row['codon_change']
        if codon_change in synonymous_changes:
            total_count = row['totalcount']
            alt_count = row['altcount']
            extra_annotation = row['custom_annotation']
            phylop_value = row['phyloP']
            phastcons_value = row['phastCons']

            # Skip if the SNP is in an exon-intron junction
            if extra_annotation == 'eij':
                continue

            if (row['maineffect'] == 'SYNONYMOUS_CODING' and not pd.isna(phylop_value) and not pd.isna(row['phastCons'])):
                codon_dict[codon_change][total_count].append(alt_count)
                phylop_dict[codon_change].append(phylop_value)
                phastcons_dict[codon_change].append(phastcons_value)

    return codon_dict, phylop_dict, phastcons_dict


# Define the function to process the main dictionary 
def process_main_data(df: DataFrame) -> dict:
    """
    Function to process main data dictionary.
    """
    codon_data = defaultdict(lambda: defaultdict(list))
    for _, row in df.iterrows():
        codon_change = row['codon_change']
        total_count = row['totalcount']
        alt_count = row['altcount']
        custom_annotation = row['custom_annotation']
        phylop_value = row['phyloP']
        phastcons_value = row['phastCons']

        # Skip if the SNP is in an exon-intron junction
        if custom_annotation == 'eij':
            continue

        # Only include synonymous mutations with phyloP data
        if row['refaa'] == row['altaa'] and not pd.isna(phylop_value):
            codon_data[codon_change][total_count].append(alt_count)
            
            # Add phyloP and phastCons values to their respective distributions
            phylop_distribution[codon_change].append(phylop_value)
            if not pd.isna(phastcons_value):
                phastcons_distribution[codon_change].append(phastcons_value)
    
    return codon_data, phylop_distribution, phastcons_distribution


# Define the functions to save and load data
def save_data(data: dict, pickle_file: str, json_file: str):
    """
    Function to save data to pickle and json files.
    """
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    with open(json_file, 'w') as f:
        json.dump({k: dict(v) for k, v in data.items()}, f, indent=2)


def load_data(pickle_file: str) -> dict:
    """
    Function to load data from pickle file.
    """
    with open(pickle_file, 'rb') as f:
        return pickle.load(f)


# Define the function to downsample an SFS
def downsample_sfs(sfs: List[int], target_size: int) -> List[int]:
    """
    Function to downsample an SFS to a target size.
    """

    current_size = sum(sfs)
    if current_size <= target_size:
        return sfs
    p = np.array(sfs) / current_size
    return np.random.multinomial(target_size, p)


# Define the function to generate downsampled SFS
def generate_downsampled_sfs(codon_data: dict, sample_sizes: List[int]) -> dict:
    """
    Function to generate downsampled SFS.
    """
    downsampled_data = defaultdict(lambda: defaultdict(dict))
    for codon_change, size_data in codon_data.items():
        for size, counts in size_data.items():
            sfs = np.bincount(counts, minlength=size+1)
            for target_size in sample_sizes:
                if size >= target_size:
                    downsampled_sfs = downsample_sfs(sfs, target_size)
                    downsampled_data[codon_change][target_size][size] = downsampled_sfs.tolist()
    return downsampled_data


def count_codon_changes(df):
    return df['codon_change'].value_counts()


def analyze_conservation_scores(phylop_distribution, phastcons_distribution):
    phylop_stats = {}
    phastcons_stats = {}
    
    for codon_change, values in phylop_distribution.items():
        phylop_stats[codon_change] = {
            'mean': np.mean(values),
            'median': np.median(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values)
        }
    
    for codon_change, values in phastcons_distribution.items():
        phastcons_stats[codon_change] = {
            'mean': np.mean(values),
            'median': np.median(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values)
        }
    
    return phylop_stats, phastcons_stats