#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AncestryChecker Lite: A simplified version of AncestryChecker.
This program checks the ancestry of individuals in a genealogy assuming the reference
genome corresponds to the 664c founder line.

@author: Luis Javier Madrigal-Roca
@date: 2025-04-09
"""

import argparse
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
import os

# Setting up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='AncestryChecker Lite: Check ancestry assuming REF = 664c.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-r', '--relationships', required=True, help='Path to the relationship map file')
    parser.add_argument('-o', '--output', default='ancestry_results', help='Output directory')
    parser.add_argument('--region', default=None, help='Region of interest (e.g., "Chr_01" or "Chr_01:1000000-2000000")')
    return parser.parse_args()

def read_vcf(path):
    """Read a VCF file into a pandas DataFrame."""
    logger.info(f"Reading VCF file: {path}")
    
    # Parse VCF header to get column names
    with open(path, 'r') as vcf_file:
        header_line = None
        for line in vcf_file:
            if line.startswith('#CHROM'):
                header_line = line.strip()
                break
    
    if not header_line:
        raise ValueError("Could not find header line in VCF file.")
    
    # Parse header line to get column names
    columns = header_line.strip('#').split()
    
    # Read VCF data
    vcf_data = pd.read_csv(
        path, 
        comment='#', 
        sep='\t', 
        names=columns,
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 
                'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str}
    )
    
    # Rename #CHROM to CHROM
    vcf_data = vcf_data.rename(columns={'#CHROM': 'CHROM'})
    
    logger.info(f"Read {len(vcf_data)} variants from VCF file")
    return vcf_data

def read_relationship_map(path):
    """Read relationship map."""
    logger.info(f"Reading relationship map: {path}")
    
    # Read CSV with no header, and assign column names
    relationships = pd.read_csv(
        path, 
        sep=r'\s+',      # Use whitespace as separator
        header=None,     # No header row in the file
        names=['Sample', 'Generation', 'Founder1', 'Founder2'],  # Assign column names
        dtype=str
    )
    
    logger.info(f"Found {len(relationships)} samples in relationship map")
    return relationships

def extract_gt(value, format_str):
    """Extract genotype from VCF format string."""
    if pd.isna(value):
        return np.nan
        
    fields = value.split(':')
    format_fields = format_str.split(':')
    gt_idx = format_fields.index('GT') if 'GT' in format_fields else 0
    
    if gt_idx < len(fields):
        return fields[gt_idx]
    return np.nan

def filter_by_region(vcf_df, region):
    """Filter VCF data by region."""
    if not region:
        return vcf_df
    
    logger.info(f"Filtering by region: {region}")
    
    if ':' in region:
        # Format is "CHROM:START-END"
        chrom, pos_range = region.split(':')
        start, end = map(int, pos_range.split('-'))
        filtered = vcf_df[(vcf_df['CHROM'] == chrom) & 
                            (vcf_df['POS'] >= start) & 
                            (vcf_df['POS'] <= end)]
    else:
        # Just chromosome name
        filtered = vcf_df[vcf_df['CHROM'] == region]
    
    logger.info(f"Retained {len(filtered)} variants in region {region}")
    return filtered

def determine_ancestry(vcf_df, relationships, samples_to_analyze):
    """
    Determine ancestry of samples assuming REF = 664c.
    
    Args:
        vcf_df (pd.DataFrame): VCF data
        relationships (pd.DataFrame): Relationship map
        samples_to_analyze (list): List of sample names to analyze
        
    Returns:
        dict: Dictionary mapping sample names to ancestry results DataFrames
    """
    results = {}
    
    # Get format string (assuming consistent across variants)
    format_str = vcf_df['FORMAT'].iloc[0]
    
    # Extract sample columns
    format_col_idx = vcf_df.columns.get_loc('FORMAT')
    all_sample_cols = vcf_df.columns[format_col_idx + 1:]
    
    for sample in samples_to_analyze:
        if sample not in all_sample_cols:
            logger.warning(f"Sample {sample} not found in VCF, skipping...")
            continue
        
        logger.info(f"Analyzing ancestry for {sample}...")
        
        # Get parents from relationship map
        sample_info = relationships[relationships['Sample'] == sample]
        if len(sample_info) == 0:
            logger.warning(f"Sample {sample} not found in relationship map, skipping...")
            continue
            
        founder1 = sample_info['Founder1'].iloc[0]
        founder2 = sample_info['Founder2'].iloc[0]
        
        # Create result DataFrame
        result_df = pd.DataFrame({
            'CHROM': vcf_df['CHROM'],
            'POS': vcf_df['POS'],
            'REF': vcf_df['REF'],  # REF = 664c allele
            'ALT': vcf_df['ALT']
        })
        
        # Extract genotypes
        result_df['F2'] = vcf_df[sample].apply(lambda x: extract_gt(x, format_str))
        
        # We assume REF = 664c, so we only need to extract genotype for the other founder
        other_founder = founder1 if founder1 != '664c' else founder2
        if other_founder in all_sample_cols:
            result_df[other_founder] = vcf_df[other_founder].apply(lambda x: extract_gt(x, format_str))
        
        # Determine ancestry
        def assign_ancestry(row):
            if pd.isna(row['F2']):
                return 'Missing'
                
            if row['F2'] == './.' or row['F2'] == '.':
                return 'Missing'
                
            # Extract alleles
            f2_alleles = row['F2'].replace('|', '/').split('/')
            
            # Check if homozygous reference (same as 664c)
            if row['F2'] == '0/0' and row[other_founder] != '0/0':
                return '664c'
            
            if row['F2'] == '0/0' and row[other_founder] == '0/0':
                return 'Unresolved'
                
            # Check if homozygous alternate (same as other founder)
            if row['F2'] == '1/1' and other_founder in row and row[other_founder] == '1/1':
                return other_founder
                
            # Check if heterozygous (mixture)
            if '0' in f2_alleles and '1' in f2_alleles:
                return f'664c/{other_founder}'
                
            return 'Novel'
            
        result_df['Ancestry'] = result_df.apply(assign_ancestry, axis=1)
        results[sample] = result_df
        
    return results

def plot_chromosomal_painting(ancestry_df, sample_id, output_dir, color_map=None):
    """
    Create a chromosomal painting to visualize recombination events.
    
    This function plots the ancestry assignments along each chromosome,
    showing recombination breakpoints where ancestry changes.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create default color map if not provided
    if color_map is None:
        color_map = {
            '664c': '#FF9999',
            'Novel': '#66B2FF',
            'Unresolved': '#99FF99',
            'Missing': '#FFCC99'
        }
        
        # Add colors for other ancestry types
        potential_ancestry_types = set(ancestry_df['Ancestry'].unique())
        base_colors = ['#c2c2f0', '#ffb347', '#a5d6a7', '#ef9a9a', '#90caf9']
        color_index = 0
        
        for ancestry in potential_ancestry_types:
            if ancestry not in color_map:
                if color_index < len(base_colors):
                    color_map[ancestry] = base_colors[color_index]
                    color_index += 1
                else:
                    color_map[ancestry] = f"#{np.random.randint(0, 0xFFFFFF):06x}"
    
    # Get unique chromosomes and sort them
    chromosomes = sorted(ancestry_df['CHROM'].unique())
    num_chromosomes = len(chromosomes)
    
    # Create figure
    fig_height = max(8, 1 + num_chromosomes * 0.5)  # Adjust height based on number of chromosomes
    fig, axes = plt.subplots(num_chromosomes, 1, figsize=(12, fig_height), sharex=False)
    
    # If there's only one chromosome, axes won't be an array
    if num_chromosomes == 1:
        axes = [axes]
    
    # Set overall title
    fig.suptitle(f'Chromosomal Painting for {sample_id}', fontsize=16)
    
    # Create legend handles
    legend_handles = []
    for ancestry, color in color_map.items():
        if ancestry in ancestry_df['Ancestry'].values:
            legend_handles.append(plt.Rectangle((0, 0), 1, 1, color=color, label=ancestry))
    
    # Plot each chromosome
    for i, chromosome in enumerate(chromosomes):
        ax = axes[i]
        
        # Filter data for this chromosome
        chrom_data = ancestry_df[ancestry_df['CHROM'] == chromosome].copy()
        chrom_data = chrom_data.sort_values('POS')
        
        if len(chrom_data) == 0:
            ax.text(0.5, 0.5, f"No data for {chromosome}", 
                    horizontalalignment='center', verticalalignment='center')
            ax.set_ylabel(chromosome)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
            
        # Group by consecutive ancestry segments
        chrom_data['segment'] = (chrom_data['Ancestry'] != chrom_data['Ancestry'].shift()).cumsum()
        segments = chrom_data.groupby('segment')
        
        # Plot segments
        for _, segment in segments:
            start_pos = segment['POS'].min()
            end_pos = segment['POS'].max()
            ancestry = segment['Ancestry'].iloc[0]
            color = color_map.get(ancestry, '#DDDDDD')
            
            # Plot the segment
            ax.axvspan(start_pos, end_pos, color=color, alpha=0.7)
            
            # Add text label if segment is large enough
            segment_width = end_pos - start_pos
            chromosome_span = chrom_data['POS'].max() - chrom_data['POS'].min()
            if segment_width / chromosome_span > 0.05:  # Only label segments that are >5% of chromosome length
                mid_pos = (start_pos + end_pos) / 2
                ax.text(mid_pos, 0.5, ancestry, 
                       horizontalalignment='center', verticalalignment='center',
                       fontsize=8, fontweight='bold')
        
        # Set axis labels and limits
        ax.set_xlim(chrom_data['POS'].min(), chrom_data['POS'].max())
        ax.set_ylim(0, 1)
        ax.set_ylabel(chromosome)
        
        # Remove y-ticks
        ax.set_yticks([])
        
        # Format x-axis with Mbp scale
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x/1000000:.1f}'))
        
        # Only show x-label for the bottom subplot
        if i == num_chromosomes - 1:
            ax.set_xlabel('Position (Mbp)')
        
    # Add legend to the figure
    fig.legend(handles=legend_handles, loc='upper center', ncol=min(len(color_map), 4),
               bbox_to_anchor=(0.5, 0.05), frameon=False)
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Adjust for the legend
    plt.savefig(f"{output_dir}/{sample_id}_chromosomal_painting.png", dpi=300)
    plt.close()
    
    logger.info(f"Saved chromosomal painting for {sample_id}")

def plot_ancestry(ancestry_df, sample_id, output_dir):
    """Plot ancestry results."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Count ancestry by chromosome
    chrom_counts = {}
    for chrom in ancestry_df['CHROM'].unique():
        chrom_data = ancestry_df[ancestry_df['CHROM'] == chrom]
        counts = chrom_data['Ancestry'].value_counts().to_dict()
        chrom_counts[chrom] = counts
    
    # Overall counts
    total_counts = ancestry_df['Ancestry'].value_counts()
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Define a consistent color mapping for ancestry types
    color_map = {
        '664c': '#FF9999',
        'Novel': '#66B2FF',
        'Unresolved': '#99FF99',
        'Missing': '#FFCC99'
    }
    
    # Add colors for other ancestry types that might be present
    # (like founder names or mixed ancestry)
    potential_ancestry_types = set(ancestry_df['Ancestry'].unique())
    base_colors = ['#c2c2f0', '#ffb347', '#a5d6a7', '#ef9a9a', '#90caf9']
    color_index = 0
    
    for ancestry in potential_ancestry_types:
        if ancestry not in color_map:
            if color_index < len(base_colors):
                color_map[ancestry] = base_colors[color_index]
                color_index += 1
            else:
                # Generate a random color if we run out of predefined colors
                color_map[ancestry] = f"#{np.random.randint(0, 0xFFFFFF):06x}"
    
    # Plot overall proportions with consistent colors
    pie_colors = [color_map.get(ancestry, '#DDDDDD') for ancestry in total_counts.index]
    total_counts.plot.pie(
        ax=ax1, 
        autopct='%1.1f%%',
        startangle=90,
        colors=pie_colors
    )
    ax1.set_title(f'Overall Ancestry Proportions for {sample_id}')
    ax1.set_ylabel('')
    
    # Plot by chromosome with consistent colors
    chroms = sorted(chrom_counts.keys())
    ancestry_types = sorted(set().union(*[set(counts.keys()) for counts in chrom_counts.values()]))
    
    x = np.arange(len(chroms))
    width = 0.8 / len(ancestry_types)
    
    for i, ancestry in enumerate(ancestry_types):
        values = [chrom_counts[chrom].get(ancestry, 0) for chrom in chroms]
        ax2.bar(x + i * width, values, width, label=ancestry, color=color_map.get(ancestry, '#DDDDDD'))
    
    ax2.set_title(f'Ancestry Composition by Chromosome for {sample_id}')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Number of SNPs')
    ax2.set_xticks(x + width * (len(ancestry_types) - 1) / 2)
    ax2.set_xticklabels(chroms, rotation=45)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{sample_id}_ancestry.png")
    plt.close()
    
    logger.info(f"Saved ancestry plot for {sample_id}")
    
    # Return the color map so it can be reused for chromosome painting
    return color_map

def main():
    """Main function."""
    args = parse_arguments()
    
    # Read input files
    vcf_data = read_vcf(args.vcf)
    relationships = read_relationship_map(args.relationships)
    
    # Filter by region if specified
    if args.region:
        vcf_data = filter_by_region(vcf_data, args.region)
    
    # Get F2 samples to analyze
    f2_samples = relationships[relationships['Generation'] == 'F2']['Sample'].tolist()
    logger.info(f"Found {len(f2_samples)} F2 samples to analyze")
    
    # Determine ancestry
    ancestry_results = determine_ancestry(vcf_data, relationships, f2_samples)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Plot results and save summaries
    for sample, result_df in ancestry_results.items():
        logger.info(f"Processing results for {sample}...")
        
        # Save detailed results to CSV
        result_df.to_csv(f"{args.output}/{sample}_ancestry_details.csv", index=False)
        
        # Plot ancestry
        color_map = plot_ancestry(result_df, sample, args.output)
        
        # Plot chromosomal painting to visualize recombination
        plot_chromosomal_painting(result_df, sample, args.output, color_map)
        
        # Print summary statistics
        ancestry_counts = result_df['Ancestry'].value_counts()
        logger.info(f"Ancestry composition for {sample}:")
        for ancestry, count in ancestry_counts.items():
            percentage = 100 * count / len(result_df)
            logger.info(f"  {ancestry}: {count} SNPs ({percentage:.2f}%)")

    logger.info(f"Results saved to {args.output} directory")

if __name__ == "__main__":
    main()