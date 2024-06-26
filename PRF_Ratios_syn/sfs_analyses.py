"""
Module containing functions for performing analyses of synonymous SFSs.
"""

import pickle
from typing import List
import numpy as np


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