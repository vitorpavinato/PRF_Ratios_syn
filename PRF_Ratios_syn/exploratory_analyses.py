"""
Module for perfoming exploratory analyses.
"""

import os
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pandas import DataFrame


def make_score_counts_data(data: DataFrame) -> DataFrame:
    """
    Make average data for ploting
    """
    stats_data = data.groupby('codon_change').agg({
        'phyloP': 'mean',
        'phastCons': 'mean'
    }).reset_index()

    # Add the count separately
    codon_counts = data['codon_change'].value_counts().reset_index()
    codon_counts.columns = ['codon_change', 'count']

    # Merge the count with the other statistics
    stats_data = stats_data.merge(codon_counts, on='codon_change')

    return stats_data


def plot_codon_change_frequency_vs_scores(data: DataFrame, name: str) -> dict:
    """
    Calculate Spearman Correlation and return
    the plot for PhyloP and phastCons scores.
    """
    
    # Calculate Spearman's rank correlation
    correlation_phylop, p_value_phylop = stats.spearmanr(data['count'], data['phyloP'])
    correlation_phastcons, p_value_phastcons = stats.spearmanr(data['count'], data['phastCons'])

    # Print correlation results
    print(f"Spearman's rank correlation coefficient (phyloP): {correlation_phylop}")
    print(f"P-value (phyloP): {p_value_phylop}")
    print(f"Spearman's rank correlation coefficient (phastCons): {correlation_phastcons}")
    print(f"P-value (phastCons): {p_value_phastcons}")

    # Create scatter plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # phyloP plot
    sns.regplot(x='count', y='phyloP', data=data, ax=ax1)
    ax1.set_xlabel('Frequency of Codon Change')
    ax1.set_ylabel('Average phyloP Score')
    ax1.set_title('Codon Change Frequency vs Average phyloP Score')

    # Add correlation info to phyloP plot
    ax1.annotate(f'Spearman r = {correlation_phylop:.3f}\np-value = {p_value_phylop:.3e}',
                 xy=(0.05, 0.95), xycoords='axes fraction',
                 ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                 fontsize=12)

    # phastCons plot
    sns.regplot(x='count', y='phastCons', data=data, ax=ax2)
    ax2.set_xlabel('Frequency of Codon Change')
    ax2.set_ylabel('Average phastCons Score')
    ax2.set_title('Codon Change Frequency vs Average phastCons Score')

    # Add correlation info to phastCons plot
    ax2.annotate(f'Spearman r = {correlation_phastcons:.3f}\np-value = {p_value_phastcons:.3e}',
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
    fig_path = f'../exploratory/{name}_phylop_phastcons_codon_density_{timestamp}.png'
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"Combined plot saved as {fig_path}")

    # # Save individual subplots
    # phylop_path = f'../exploratory/{name}_phylop_codon_density_{timestamp}.png'
    # phastcons_path = f'../exploratory/{name}_phastcons_codon_density_{timestamp}.png'
    
    # ax1.figure.savefig(phylop_path, dpi=300, bbox_inches='tight')
    # ax2.figure.savefig(phastcons_path, dpi=300, bbox_inches='tight')
    
    # print(f"PhyloP plot saved as {phylop_path}")
    # print(f"PhastCons plot saved as {phastcons_path}")

    # Show the plot
    plt.show()

    # Close the figure to free up memory
    plt.close(fig)

    # Return correlation results
    return {
        'phyloP': {'correlation': correlation_phylop, 'p_value': p_value_phylop},
        'phastCons': {'correlation': correlation_phastcons, 'p_value': p_value_phastcons}
    }


def plot_scores(data: DataFrame, score_type: str, name: str) -> None:
    """
    # Create a function to plot the scores.
    """
    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create the plot
    sns.barplot(x='codon_change', y=score_type, data=data)
    
    # Set title and labels
    plt.title(f'Average {score_type} Score by Codon Change')
    plt.xlabel('Codon Change')
    plt.ylabel(f'Average {score_type} Score')
    
    # Rotate x-axis labels
    plt.xticks(rotation=90)
    
    # Adjust layout
    plt.tight_layout()
    
    # Generate timestamp
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create directory if it doesn't exist
    os.makedirs('../exploratory', exist_ok=True)

    # Save the figure
    file_path = f'../exploratory/{name}_{score_type}_plot_{timestamp}.png'
    fig.savefig(file_path, dpi=300, bbox_inches='tight')

    print(f"Plot saved as {file_path}")

    # Show the plot (optional, comment out if not needed)
    plt.show()

    # Close the figure to free up memory
    plt.close(fig)


def score_histogram(data: DataFrame, score: str) -> None:
    plt.figure(figsize=(10, 6))
    sns.histplot(data[score], kde=True)
    plt.xlabel(f'{score} Score')
    plt.ylabel('Count')
    plt.title(f'Distribution of {score} Scores for Synonymous Mutations')
    plt.show()
