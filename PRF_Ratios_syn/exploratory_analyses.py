"""
Module for perfoming exploratory analyses.
"""

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pandas import DataFrame
from codon_changes_dict import synonymous_1nt_pairs


def make_groupby_table(
        df: DataFrame,
        custom_annotation_col: str,
        phylop_col: str,
        phastcons_col: str,
        use_filter: bool) -> DataFrame:
    """
    Make a table of scores by aggregating/grouping by codon change.
    For score parsed in scores_names, grouping will calculate the mean.
    """

    if use_filter:
        df = df[df[custom_annotation_col] != 'eij']

    # Create the table
    agg_table = df.groupby('codon_change').agg({
        phylop_col: 'mean',
        phastcons_col: 'mean'
    }).reset_index()

    # Add the count separately
    codon_counts = df['codon_change'].value_counts().reset_index()
    codon_counts.columns = ['codon_change', 'count']

    # Merge the count with the other statistics
    agg_table = agg_table.merge(codon_counts, on='codon_change')

    return agg_table


def plot_xvar_vs_yvars(
        dt: DataFrame,
        xvar_col: str,
        y1var_col: str,
        y2var_col: str,
        name: str) -> dict:
    """
    Calculate Spearman Correlation and return
    the plot for PhyloP and phastCons scores.
    """

    # Calculate Spearman's rank correlation
    correlation_y1var, p_value_y1var = stats.spearmanr(dt[xvar_col], dt[y1var_col])
    correlation_y2var, p_value_y2var = stats.spearmanr(dt[xvar_col], dt[y2var_col])

    # Print correlation results
    print(f"Spearman's rank correlation coefficient {y1var_col}: {correlation_y1var}")
    print(f"P-value {y1var_col}: {p_value_y1var}")
    print(f"Spearman's rank correlation coefficient {y2var_col}: {correlation_y2var}")
    print(f"P-value {{y2var_col}}: {p_value_y2var}")

    # Create scatter plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # y1var plot
    sns.regplot(x=xvar_col, y=y1var_col, data=dt, ax=ax1)
    ax1.set_xlabel(f'{xvar_col}')
    ax1.set_ylabel(f'{y1var_col}')
    ax1.set_title(f'{xvar_col} vs {y1var_col}')

    # Add correlation info to phyloP plot
    ax1.annotate(f'Spearman r = {correlation_y1var:.3f}\np-value = {p_value_y1var:.3e}',
                 xy=(0.05, 0.95), xycoords='axes fraction',
                 ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                 fontsize=12)

    # y1var plot
    sns.regplot(x=xvar_col, y=y2var_col, data=dt, ax=ax1)
    ax2.set_xlabel(f'{xvar_col}')
    ax2.set_ylabel(f'{y2var_col}')
    ax2.set_title(f'{xvar_col} vs {y2var_col}')

    # Add correlation info to phyloP plot
    ax2.annotate(f'Spearman r = {correlation_y2var:.3f}\np-value = {p_value_y2var:.3e}',
                 xy=(0.05, 0.95), xycoords='axes fraction',
                 ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                 fontsize=12)

    plt.tight_layout()

    # Generate timestamp for file names
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create directory if it doesn't exist
    os.makedirs('../exploratory', exist_ok=True)

    # Save the entire figure
    fig_path = f'../exploratory/{name}_{y1var_col}_{y2var_col}_{xvar_col}_{timestamp}.png'
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"Combined plot saved as {fig_path}")

    # Show the plot
    plt.show()

    # Close the figure to free up memory
    plt.close(fig)

    # Return correlation results
    return {
        'y1var': {'correlation': correlation_y1var, 'p_value': p_value_y1var},
        'y2var': {'correlation': correlation_y2var, 'p_value': p_value_y2var}
    }


def plot_scores(
        dt: DataFrame,
        score_name: str,
        name: str) -> None:
    """
    # Create a function to plot the scores.
    """
    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create the plot
    sns.barplot(x='codon_change', y=score_name, data=dt)

    # Set title and labels
    plt.title(f'Average {score_name} Score by Codon Change')
    plt.xlabel('Codon Change')
    plt.ylabel(f'Average {score_name} Score')

    # Rotate x-axis labels
    plt.xticks(rotation=90)

    # Adjust layout
    plt.tight_layout()

    # Generate timestamp
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create directory if it doesn't exist
    os.makedirs('../exploratory', exist_ok=True)

    # Save the figure
    file_path = f'../exploratory/{name}_{score_name}_plot_{timestamp}.png'
    fig.savefig(file_path, dpi=300, bbox_inches='tight')

    print(f"Plot saved as {file_path}")

    # Show the plot (optional, comment out if not needed)
    plt.show()

    # Close the figure to free up memory
    plt.close(fig)


def score_histogram(dt: DataFrame, score_col_name: str) -> None:
    """
    Make a histogram for a given score.
    """

    # Plot the histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(dt[score_col_name], kde=True)
    plt.xlabel(f'{score_col_name} Score')
    plt.ylabel('Count')
    plt.title(f'Distribution of {score_col_name} Scores for Synonymous Mutations')
    plt.show()


def create_codon_stats(
        df: DataFrame,
        use_filter=True,
        phylop_col='phyloP',
        phastcons_col='phastCons',
        custom_annotation_col='custom_annotation'):
    """
    Create a dictionary of codon change statistics from a DataFrame.

    Parameters:
    df (pandas.DataFrame): The input DataFrame containing codon change data.
    use_filter (bool): If True, filter out rows where custom_annotation is 'eij'.
    phylop_col (str): Name of the column containing phyloP scores.
    phastcons_col (str): Name of the column containing phastCons scores.
    custom_annotation_col (str): Name of the column containing custom annotations.

    Returns:
    dict: A dictionary containing statistics for each codon change.
    """

    # Check if required columns exist
    required_cols = ['codon_change', phylop_col, phastcons_col]
    
    if use_filter:
        required_cols.append(custom_annotation_col)

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

    # Generate all possible synonymous codon changes
    synonymous_changes = synonymous_1nt_pairs

    # Create the dictionary
    codon_stats = {change: {'phyloP': [], 'phastCons': [], 'count': 0} for change in synonymous_changes}

    if use_filter:
        df = df[df[custom_annotation_col] != 'eij']

    # Fill the dictionary with data from the DataFrame
    for _, row in df.iterrows():
        codon_change = row['codon_change']

        if codon_change in synonymous_changes:
            phylop = row[phylop_col]
            phastcons = row[phastcons_col]

            if not np.isnan(phylop) and not np.isnan(phastcons):
                codon_stats[codon_change]['phyloP'].append(phylop)
                codon_stats[codon_change]['phastCons'].append(phastcons)
                codon_stats[codon_change]['count'] += 1

    # Calculate means and clean up empty entries
    for change, scores in codon_stats.copy().items():
        if scores['count'] == 0:
            del codon_stats[change]
        else:
            scores['mean_phyloP'] = np.mean(scores['phyloP']) if scores['phyloP'] else np.nan
            scores['mean_phastCons'] = np.mean(scores['phastCons']) if scores['phastCons'] else np.nan
            del scores['phyloP']
            del scores['phastCons']

    return codon_stats


def create_codon_stats_dataframe(codon_stats):
    """
    Create a pandas DataFrame from the codon_stats dictionary.

    Parameters:
    codon_stats (dict): A dictionary containing statistics for each codon change.

    Returns:
    pandas.DataFrame: A DataFrame summarizing the codon change statistics.
    """
    data = []
    for codon_change, values in codon_stats.items():
        data.append({
            'codon_change': codon_change,
            'mean_phyloP': values['mean_phyloP'],
            'mean_phastCons': values['mean_phastCons'],
            'count': values['count']
        })

    df_codon_stats = pd.DataFrame(data)

    # # Sort the DataFrame by count in descending order
    # df_codon_stats = df_codon_stats.sort_values('count', ascending=False)

    # # Reset the index
    # df_codon_stats = df_codon_stats.reset_index(drop=True)

    return df_codon_stats


# Function to get the reverse codon change
def get_reverse(codon_change: str) -> str:
    """
    Get the reverse codon change.
    """
    return '->'.join(codon_change.split('->')[::-1])


# Create a dictionary to store the results
def create_mean_dict(df: DataFrame) -> dict:
    """
    Create a dictionary with the mean phyloP and phastCons scores for each codon change.
    """

    # Create an empty dictionary
    result = {}

    # Process each codon change
    for _, row in df.iterrows():
        codon_change = row['codon_change']
        reverse_change = get_reverse(codon_change)
        # directions = (codon_change, reverse_change)

        if codon_change not in result:
            # Find the reverse change in the DataFrame
            reverse_row = df[df['codon_change'] == reverse_change]

            if not reverse_row.empty:
                result[codon_change] = {
                    'mean_phyloP': np.mean([row['mean_phyloP'], reverse_row['mean_phyloP'].values[0]]),
                    'mean_phastCons': np.mean([row['mean_phastCons'], reverse_row['mean_phastCons'].values[0]]),
                    'mean_counts': np.mean([row['count'], reverse_row['count'].values[0]])
                }

    return result


def normalize_codon_pair(pair: str) -> str:
    """
    Normalize the pair column in DataFrames
    """
    codons = pair.replace('[', '').replace(']', '').split(',') if '[' in pair else pair.split(':')
    return '->'.join(sorted(codons))


def normalize_dataframe(df: DataFrame) -> dict:
    """
    Normalize the pair column in DataFrames
    """

    # Applyt the normalize_codon_pair function to the 'pair' column
    df['pair'] = df['pair'].apply(normalize_codon_pair)

    # Make a dictionary from the DataFrame
    df_to_dict = {}
    for _, row in df.iterrows():
        df_to_dict[row['pair']] = {
            'codon_rate': row['codon_rate'],
            'aa_rate': row['aa_rate']
        }
    
    return df_to_dict


def combine_dicts(dict1, dict2) -> DataFrame:
    """
    Combine two dictionaries
    """

    # Make a copy of dict1
    dict1 = dict1.copy()

    # Combine the dictionaries
    for key, value in dict2.items():
        if key in dict1:
            dict1[key].update(value)

    filtered_dict = {}

    for key, value in dict1.items():
        if len(value) > 3:
            filtered_dict[key] = value

    # Conver the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(filtered_dict, orient='index')

    # Reset the index to make the codon change a column
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'codon_change'}, inplace=True)

    return df
