"""
Scratch script for processing .TSV files.
"""

from typing import List, Tuple
from pandas import DataFrame
import pandas as pd


main_table = "data/dgrp2/tables/NC_Chr2L_tables.tsv"
extra_ann_table = "data/dgrp2/extra_ann_tables/NC_Chr2L_extra_ann.tsv"

phylop_table = "data/dgrp2/extra_ann_tables/dm6.phyloP27way_chr2L.csv"
phastcons_table = "data/dgrp2/extra_ann_tables/dm6.27way.phastCons_chr2L.csv"


main_df = pd.read_table(main_table,
                        low_memory=False,
                        keep_default_na=True,
                        na_values='NA'
                        )


extra_ann_df = pd.read_table(extra_ann_table, keep_default_na=True, na_values='NA')

phylop_df = pd.read_csv(phylop_table, sep=',')
phastcons_df = pd.read_csv(phastcons_table, sep=',')


# Define swap pairs
swap_pairs = [
    ('ref', 'alt'),
    ('refcount', 'altcount'),
    ('refcontext', 'altcontext'),
    ('refcontext_complrev', 'altcontext_complrev'),
    ('refcodon', 'altcodon'),
    ('refaa', 'altaa')
]


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
main_df_swapped_merged = merged_df.to_csv('data/dgrp2/tables/main_df_swapped_merged.tsv', sep='\t', index=False, na_rep='NA')



# Merge with score1 table
merged_df_phylop = pd.merge(merged_df, phylop_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phylop'))
merged_phylop = merged_df_phylop.to_csv('data/dgrp2/tables/main_df_swapped_merged_phylop.tsv', sep='\t', index=False, na_rep='NA')



# Merge with score2 table
merged_phascons = pd.merge(merged_df_phylop, phastcons_df, left_on=['chrom', 'pos'], right_on=['chromosome', 'position'], how='left', suffixes=('', '_phastcons'))

merged_phascons = merged_phascons.to_csv('data/dgrp2/tables/main_df_swapped_merged_phylop_phascons.tsv', sep='\t', index=False, na_rep='NA')
