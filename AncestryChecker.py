#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main program for AncestryChecker.
This program checks the ancestry of a given individual in a genealogy defined by the user.

@author: Luis Javier Madrigal-Roca

@date: 2025-04-01
@version: 1.0

"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main program for AncestryChecker.
This program checks the ancestry of a given individual in a genealogy defined by the user.

@author: Luis Javier Madrigal-Roca

@date: 2025-04-01
@version: 1.0

"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.vcf_reader_utilities import read_vcf
from src.data_tidyer import read_relationship_map, filter_vcf_data, identify_informative_snps

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='AncestryChecker: Check the ancestry of F2 individuals.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-r', '--relationships', required=True, help='Path to the relationship map file')
    parser.add_argument('-o', '--output', default='ancestry_results', help='Output directory')
    parser.add_argument('--no-missing', action='store_true', help='Do not allow missing data in samples')
    
    return parser.parse_args()

def determine_ancestry(vcf_data, relationships, sample_col):
    """
    Determine the ancestry of each variant in an F2 individual.
    
    Args:
        vcf_data (pd.DataFrame): Processed VCF data
        relationships (pd.DataFrame): Relationship map
        sample_col (str): Column name for the F2 sample
        
    Returns:
        pd.DataFrame: DataFrame with ancestry assignments
    """
    # Get the founders for this sample
    sample_id = sample_col.split('_')[0]  # Assuming the column name is in format "sample_GT"
    sample_info = relationships[relationships['Sample'] == sample_id]
    
    if len(sample_info) == 0:
        raise ValueError(f"Sample {sample_id} not found in relationship map")
    
    founder1 = sample_info['Founder1'].iloc[0]
    founder2 = sample_info['Founder2'].iloc[0]
    
    founder1_col = f"{founder1}_GT"
    founder2_col = f"{founder2}_GT"
    
    # Create a new DataFrame for ancestry results
    results = pd.DataFrame({
        'CHROM': vcf_data['CHROM'],
        'POS': vcf_data['POS'],
        'REF': vcf_data['REF'],
        'ALT': vcf_data['ALT'],
        'F2': vcf_data[sample_col],
        'Founder1': vcf_data[founder1_col],
        'Founder2': vcf_data[founder2_col]
    })
    
    # Determine ancestry
    def assign_ancestry(row):
        if pd.isna(row['F2']):
            return 'Missing'
        elif row['F2'] == row['Founder1'] and row['F2'] != row['Founder2']:
            return 'Founder1'
        elif row['F2'] == row['Founder2'] and row['F2'] != row['Founder1']:
            return 'Founder2'
        elif row['F2'] == row['Founder1'] and row['F2'] == row['Founder2']:
            return 'Both'  # Uninformative
        else:
            return 'Novel'  # New variant not in founders
    
    results['Ancestry'] = results.apply(assign_ancestry, axis=1)
    return results

def plot_ancestry(ancestry_data, sample_id, output_dir):
    """
    Create a visualization of ancestry across chromosomes.
    
    Args:
        ancestry_data (pd.DataFrame): Ancestry assignment data
        sample_id (str): ID of the F2 sample
        output_dir (str): Output directory
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get unique chromosomes
    chromosomes = ancestry_data['CHROM'].unique()
    
    plt.figure(figsize=(15, 10))
    
    # Color mapping
    colors = {'Founder1': 'blue', 'Founder2': 'red', 'Both': 'purple', 'Novel': 'green', 'Missing': 'gray'}
    
    for i, chrom in enumerate(chromosomes):
        chrom_data = ancestry_data[ancestry_data['CHROM'] == chrom]
        
        # Plot positions
        plt.subplot(len(chromosomes), 1, i+1)
        
        for ancestry, color in colors.items():
            subset = chrom_data[chrom_data['Ancestry'] == ancestry]
            plt.scatter(subset['POS'], [i] * len(subset), c=color, s=5, label=ancestry if i == 0 else "")
        
        plt.title(f'Chromosome {chrom}')
        plt.ylabel('Ancestry')
        plt.xlim(min(chrom_data['POS']), max(chrom_data['POS']))
        plt.yticks([])
    
    # Add legend to the first subplot
    plt.subplot(len(chromosomes), 1, 1)
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_id}_ancestry.png'))
    plt.close()
    
    # Also save the data to a CSV file
    ancestry_data.to_csv(os.path.join(output_dir, f'{sample_id}_ancestry.csv'), index=False)

def main():
    """Main function for the AncestryChecker program."""
    args = parse_arguments()
    
    # Read input files
    print(f"Reading VCF file: {args.vcf}")
    vcf_data = read_vcf(args.vcf)
    
    print(f"Reading relationship map: {args.relationships}")
    relationships = read_relationship_map(args.relationships)
    
    # Extract sample names
    format_col_idx = vcf_data.columns.get_loc('FORMAT')
    sample_cols = vcf_data.columns[format_col_idx + 1:]
    
    # Get all founder names from the relationship map
    founders = set(relationships['Founder1'].tolist() + relationships['Founder2'].tolist())
    
    # Get all F2 samples from the relationship map
    f2_samples = relationships[relationships['Generation'] == 'F2']['Sample'].tolist()
    
    print(f"Found {len(founders)} founders and {len(f2_samples)} F2 samples")
    
    # Filter VCF data
    print("Filtering VCF data...")
    filtered_vcf = filter_vcf_data(vcf_data, list(founders), f2_samples, not args.no_missing)
    
    # Identify informative SNPs
    print("Identifying informative SNPs...")
    founder_cols = [col for col in sample_cols if col in founders]
    informative_vcf = identify_informative_snps(filtered_vcf, founder_cols)
    
    print(f"Retained {len(informative_vcf)} informative SNPs out of {len(vcf_data)}")
    
    # Process each F2 sample
    for f2_id in f2_samples:
        f2_col = f"{f2_id}_GT"
        
        if f2_col in informative_vcf.columns:
            print(f"Analyzing ancestry for sample {f2_id}...")
            ancestry_results = determine_ancestry(informative_vcf, relationships, f2_col)
            
            print(f"Creating plots for sample {f2_id}...")
            plot_ancestry(ancestry_results, f2_id, args.output)
            
            # Print summary statistics
            ancestry_counts = ancestry_results['Ancestry'].value_counts()
            print(f"Ancestry composition for {f2_id}:")
            for ancestry, count in ancestry_counts.items():
                percentage = 100 * count / len(ancestry_results)
                print(f"  {ancestry}: {count} SNPs ({percentage:.2f}%)")
        else:
            print(f"Warning: Sample {f2_id} not found in VCF data, skipping...")
    
    print(f"Results saved to {args.output} directory")

if __name__ == "__main__":
    main()