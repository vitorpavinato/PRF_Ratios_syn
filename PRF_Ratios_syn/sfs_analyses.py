"""
Module containing functions for performing analyses of synonymous SFSs.
"""

import pickle
import json
from typing import List
import numpy as np
from scipy.stats import hypergeom


# Load data
def load_data(pickle_file: str) -> dict:
    """
    Function to load data from pickle file.
    """
    with open(pickle_file, 'rb') as f:
        return pickle.load(f)


# Downsample SFS
def downsample_sfs(
    original_sfs: list[int],
    original_size: int,
    sample_size: int
) -> list[int | float]:
    """
    Project a distribution an unfolded or folded site-frequency spectrum  
    to a sample size m < n.
    :param original_sfs: List of counts for each bin in the original SFS.
    Expects an SFS as a list of length sample size (n) + 1,
    with position i for i gene copies;
    :param popsize: Total number of gene copies in the original distribution.
    :param sampsize: Total number of gene copies in the
    downsampled distribution (m < n).
    :return: List of expected counts for each bin in the
    projected distribution.
    """

    # Raise error for an empty SFS
    if not original_sfs:
        raise ValueError("SFS is empty")

    # Create an empty sample SFS with sample_size + 1 bins
    sample_sfs = [0]*(sample_size+1)

    for pi, count in enumerate(original_sfs):
        for si in range(sample_size+1):
            prob = hypergeom.pmf(si, original_size, pi, sample_size)
            sample_sfs[si] += prob * count

    return sample_sfs


# Downsample codon change SFS in a dictionary
def downsample_codon_change_sfs_in_dict(
    codon_change_sfs_dict: dict, target_sample_sizes: List[int]
) -> dict:
    """
    Function to downsample codon change SFSs in a dictionary.
    You can specify one or more target sample sizes.
    """

    # Check if the dictionary is empty
    if not codon_change_sfs_dict:
        raise ValueError("SFS dictionary is empty")

    # Check if the target sample sizes are empty
    if not target_sample_sizes:
        raise ValueError("Target sample sizes list is empty")

    # Define the target sample sizes
    sample_sizes = target_sample_sizes

    # Create an empty dictionary to fill
    targeted_sizes_sfs_dict = {}

    for codon_change, size_data in codon_change_sfs_dict.items():
        targeted_sizes_sfs_dict[codon_change] = {}

        for sample_size in sample_sizes:
            # List of downsampled SFSs
            list_ds_sfs = []

            # Now donwsample only SFSs with key values higher than sample_size
            for nsize, nsfs in size_data.items():
                if nsize >= sample_size:
                    ds_sfs = downsample_sfs(nsfs, nsize, sample_size)
                    list_ds_sfs.append(ds_sfs)

            # Conver the list of SFS to an np.array
            sfs_array = np.array(list_ds_sfs)

            # Sum column-wise to get the final SFS
            tsfs = list(np.sum(sfs_array, 0))

            targeted_sizes_sfs_dict[codon_change][sample_size] = tsfs

    return targeted_sizes_sfs_dict


# Define the functions to save and load data
def save_data(data: dict, pickle_file: str, json_file: str):
    """
    Function to save data to pickle and json files.
    """
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    with open(json_file, 'w') as f:
        json.dump({k: dict(v) for k, v in data.items()}, f, indent=2)
