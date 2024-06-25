"""
Module for processing tsv files.
"""

from typing import List, Tuple
import pandas as pd
from pandas import DataFrame


# Data processing functions
def swap_values(row: DataFrame, swap_pairs: List[Tuple[str, str]]) -> DataFrame:
    """
    Function to swap values in a paired columns.
    It can be used to swap reference and alternative values in the .TSV files.
    It expects a list of tupleswith the column names (paired) to swap.
    Swapping is performed in place, row by row, and in the whole dataframe
    using apply() function.
    """
    if row['aainfo'] == 'root_alt':
        for col1, col2 in swap_pairs:
            row[col1], row[col2] = row[col2], row[col1]
    return row


def create_codon_change(row: DataFrame) -> str:
    """
    Function to create codon change string.
    """
    ref = row['refcodon']
    alt = row['altcodon']
    if pd.isna(ref) or pd.isna(alt):
        return 'NA'
    return f"{ref.upper()}->{alt.upper()}"


def process_main_table(df: DataFrame, swap_pairs: List[Tuple[str, str]]) -> DataFrame:
    """
    Function to process main table.
    """
    df = df.apply(lambda row: swap_values(row, swap_pairs), axis=1)
    df['codon_change'] = df.apply(create_codon_change, axis=1)
    return df


def merge_tables(
    main_df: DataFrame,
    extra_annotation_df: DataFrame,
    phylop_df: DataFrame,
    phastcons_df: DataFrame
) -> DataFrame:
    """
    Function to merge tables.
    It expects 4 tables: main_df, extra_annotation_df, phylop_df, phastcons_df.
    It returns a merged dataframe.
    """
    
    # Merge with custom_annotation table
    merged_df = pd.merge(main_df, extra_annotation_df, left_on=['chrom', 'pos'], right_on=['chrom', 'position'], how='left')

    # Merge with phyloP table
    merged_df = pd.merge(merged_df, phylop_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phylop'))

    # Merge with score2 table
    merged_df = pd.merge(merged_df, phastcons_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phastcons'))

    # Step 4: Clean up column names
    merged_df = merged_df.drop(columns=['position', 'chromosome', 'start', 'step', 'span', 'position_phylop','chromosome_phastcons', 'start_phastcons', 'step_phastcons', 'span_phastcons', 'position_phastcons'])
    merged_df = merged_df.rename(columns={'score': 'phyloP', 'score_phastcons': 'phastCons'})

    # Fill NaN values in custom_annotation, phyloP, and phastCons columns with 'NA'
    columns_to_fill = ['custom_annotation', 'phyloP', 'phastCons']
    merged_df[columns_to_fill] = merged_df[columns_to_fill].fillna('NA')

    return merged_df


def process_chromosome(
    main_table: str,
    extra_annotation_table: str,
    phylop_file, phastcons_file: str,
    swap_pairs: List[Tuple[str, str]]
) -> DataFrame:
    """
    Function to process a single chromosome.
    It expects 4 files: main_table, extra_annotation_table, phylop_file, phastcons_file.
    It returns a merged dataframe.
    """
    # Read files
    main_df = pd.read_table(main_table, low_memory=False, keep_default_na=True, na_values='NA')
    extra_annotation_df = pd.read_table(extra_annotation_table, keep_default_na=True, na_values='NA')
    phylop_df = pd.read_csv(phylop_file, sep=',')
    phastcons_df = pd.read_csv(phastcons_file, sep=',')

    # Process main table
    processed_df = process_main_table(main_df, swap_pairs)

    # Merge all dataframes
    merged_df = merge_tables(processed_df, extra_annotation_df, phylop_df, phastcons_df)

    return merged_df


def process_all_chromosomes(
    chromosome_files: List[Tuple[str, str, str, str]], 
    swap_pairs: List[Tuple[str, str]]
) -> DataFrame:
    """
    Function to process all chromosomes.
    """
    all_data = []
    for files in chromosome_files:
        chromosome_data = process_chromosome(*files, swap_pairs)
        all_data.append(chromosome_data)

    # Combine all chromosome data
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df
