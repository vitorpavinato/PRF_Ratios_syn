"""
Main function to process .TSV files.
"""

import pandas as pd
from data_processing import process_main_table, process_main_data, save_data, generate_downsampled_sfs


def main():
    # Define swap pairs
    swap_pairs = [
        ('ref', 'alt'),
        ('refcount', 'altcount'),
        ('refcontext', 'altcontext'),
        ('refcontext_complrev', 'altcontext_complrev'),
        ('refcodon', 'altcodon'),
        ('refaa', 'altaa')
    ]

    # Read and process main table
    main_df = pd.read_csv('main_table.csv', sep='\t')
    processed_df = process_main_table(main_df, swap_pairs)

    # Process main data
    codon_data = process_main_data(processed_df)

    # Save data
    save_data(codon_data, 'codon_sfs_data.pickle', 'codon_sfs_data.json')

    # Generate downsampled SFS
    sample_sizes = [150, 100, 75, 50, 25]
    downsampled_data = generate_downsampled_sfs(codon_data, sample_sizes)

    # Save downsampled data
    save_data(downsampled_data, 'downsampled_sfs_data.pickle', 'downsampled_sfs_data.json')

    print("Processing complete. Data saved in pickle and JSON formats.")


if __name__ == "__main__":
    main()
