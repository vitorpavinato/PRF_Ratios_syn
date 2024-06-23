"""
Module for processing tsv files.
"""

from typing import List, Tuple
import pandas as pd
from pandas import DataFrame


# Initial data processing functions
# # Function to swap values in a column
def swap_values(row: DataFrame, swap_pairs: List[Tuple[str, str]]) -> DataFrame:
    """
    Function to swap values in a column.
    It is used to swap reference and alternative values in the .TSV files.
    It expects a list of tuples with the column names to swap.
    Swapping is performed in place, row by row, and in the whole dataframe
    using apply() function.
    """
    if row['aainfo'] == 'root_alt':
        for col1, col2 in swap_pairs:
            row[col1], row[col2] = row[col2], row[col1]
    return row


# # Function to create codon change string
def create_codon_change(row: DataFrame) -> str:
    """
    Function to create codon change string.
    """
    ref = row['refcodon']
    alt = row['altcodon']
    if pd.isna(ref) or pd.isna(alt):
        return 'NA'
    return f"{ref}->{alt}"


# # Function to incorporate swapping and codon change values in main table
def process_main_table(df: DataFrame, swap_pairs: List[Tuple[str, str]]):
    """
    Function to process main table.
    """
    df = df.apply(lambda row: swap_values(row, swap_pairs), axis=1)
    df['codon_change'] = df.apply(create_codon_change, axis=1)
    return df


# # Function to merge additioanal tables to main table
def merge_tables(
    main_df: DataFrame,
    custom_annotation_df: DataFrame,
    phylop_df: DataFrame,
    phastcons_df: DataFrame
) -> DataFrame:
    """
    Function to merge tables.
    """
    
    # Merge with custom_annotation table
    merged_df = pd.merge(main_df, custom_annotation_df, left_on=['chrom', 'pos'], right_on=['chrom', 'position'], how='left')

    # Merge with score1 table
    merged_df = pd.merge(merged_df, phylop_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phylop'))

    # Merge with score2 table
    merged_df = pd.merge(merged_df, phastcons_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phastcons'))

    # Step 4: Clean up column names
    merged_df = merged_df.drop(columns=['chromosome', 'start', 'step', 'span', 'position'])
    merged_df = merged_df.rename(columns={'score_phylop': 'phyloP', 'score_phastcons': 'phastCons'})

    # Fill NaN values in custom_annotation, score1, and score2 columns with 'NA'
    columns_to_fill = ['custom_annotation', 'phyloP', 'phastCons']
    merged_df[columns_to_fill] = merged_df[columns_to_fill].fillna('NA')

    return merged_df


def process_chromosome(main_file, custom_annotations_file, phylop_file, phastcons_file, swap_pairs):
    # Read files
    main_df = pd.read_csv(main_file, sep='\t')
    custom_annotations_df = pd.read_csv(custom_annotations_file, sep='\t')
    phylop_df = pd.read_csv(phylop_file, sep='\t')
    phastcons_df = pd.read_csv(phastcons_file, sep='\t')

    # Process main table
    processed_df = process_main_table(main_df, swap_pairs)

    # Merge all dataframes
    merged_df = merge_tables(processed_df, custom_annotations_df, phylop_df, phastcons_df)

    return merged_df


def process_all_chromosomes(chromosome_files, swap_pairs):
    all_data = []
    for files in chromosome_files:
        chromosome_data = process_chromosome(*files, swap_pairs)
        all_data.append(chromosome_data)
    
    # Combine all chromosome data
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df


def analyze_synonymous_mutations(df):
    """
    Function to analyze synonymous mutations.
    """
    synonymous = df[df['refaa'] == df['altaa']]
    return synonymous[['refaa', 'altaa', 'codon_change']]


