"""
Scratch script for processing .TSV files.
"""

from typing import List, Tuple
from pandas import DataFrame
import pandas as pd


main_table = "data/dgrp2/tables/NC_Chr2L_tables.tsv"
extra_ann_table = "data/dgrp2/extra_ann_tables/NC_Chr2L_extra_ann.tsv"

phylop_table = "data/dgrp2/extra_ann_tables/dm6.phyloP27way_10klines.csv"
phastcons_table = "data/dgrp2/extra_ann_tables/dm6.27way.phastCons_10klines.csv"


# Define swap pairs
swap_pairs = [
    ('ref', 'alt'),
    ('refcount', 'altcount'),
    ('refcontext', 'altcontext'),
    ('refcontext_complrev', 'altcontext_complrev'),
    ('refcodon', 'altcodon'),
    ('refaa', 'altaa')
]


main_df = pd.read_table(main_table, low_memory=False, keep_default_na=True, na_values='NA')
extra_ann_df = pd.read_table(extra_ann_table, low_memory=False, keep_default_na=True, na_values='NA')
phylop_df = pd.read_csv(phylop_table, sep=',')
phastcons_df = pd.read_csv(phastcons_table, sep=',')

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


def create_codon_change(row: DataFrame) -> str:
    """
    Function to create codon change string.
    """
    ref = row['refcodon']
    alt = row['altcodon']
    if pd.isna(ref) or pd.isna(alt):
        return 'NA'
    return f"{ref.upper()}->{alt.upper()}"


def process_main_table(df: DataFrame, swap_pairs: List[Tuple[str, str]]):
    """
    Function to process main table.
    """
    df = df.apply(lambda row: swap_values(row, swap_pairs), axis=1)
    df['codon_change'] = df.apply(create_codon_change, axis=1)
    return df


main_df_swapped = process_main_table(main_df, swap_pairs)

main_df_swapped.to_csv('data/dgrp2/tables/main_df_swapped.tsv', sep='\t', index=False, na_rep='NA')


# Merge with custom_annotation table
merged_df = pd.merge(main_df_swapped, extra_ann_df, left_on=['chrom', 'pos'], right_on=['chrom', 'position'], how='left')

main_df_swapped.type()

print([col + ' : ' + str(main_df_swapped[col].dtype) for col in main_df_swapped.columns])

print([col + ' : ' + str(extra_ann_df[col].dtype) for col in extra_ann_df.columns])

for col in [0, 1, 2]:
    print(f"Column {col} unique values:")
    print(extra_ann_df.iloc[:, col].unique())
    print(f"Column {col} value counts:")
    print(extra_ann_df.iloc[:, col].value_counts())
    print("\n")