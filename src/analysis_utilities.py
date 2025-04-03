"""
Module for analyzing VCF data and determining ancestry of F2 individuals.

This module provides functions to read VCF files, filter data, identify informative SNPs,
determine ancestry, and visualize results. It is designed to work with genomic data from
F2 individuals and their founders.

@module: AncestryChecker
@description: This module provides functions to analyze VCF data and determine ancestry of F2 individuals.
@date: 2025-04-01

@version: 1.0

"""

import pandas as pd
import os
import matplotlib.pyplot as plt

def normalize_genotype(gt):
        if pd.isna(gt) or gt in ['.', './.', '.', '.|.']:
            return None
            
        # Extract just the GT part if it contains colons
        if ':' in gt:
            gt = gt.split(':')[0]
            
        # Remove phasing and sort alleles for consistent comparison
        gt = gt.replace('|', '/').replace('\\', '/')  # Handle different separators
        
        # Handle missing alleles
        if '.' in gt:
            return None
            
        return gt

def assign_ancestry(row):
        f2 = row['F2_alleles']
        ref = row['RefFounder_alleles']
        other = row['OtherFounder_alleles']
        
        # Handle missing data
        if f2 is None:
            return 'Missing'
        if ref is None or other is None:
            return 'Missing_Founder'
            
        # For simplicity, sort the alleles to handle phasing differences
        f2.sort()
        ref.sort()
        other.sort()
        
        # Homozygous reference in F2 (0/0)
        if f2 == [0, 0]:
            if ref == [0, 0] and other != [0, 0]:
                return 'RefFounder'  # From 664c
            elif other == [0, 0] and ref != [0, 0]:
                return 'OtherFounder'  # From other founder
            elif ref == [0, 0] and other == [0, 0]:
                return 'Both'  # Uninformative
            else:
                return 'Novel'  # Unexpected combination
        
        # Heterozygous in F2 (0/1)
        elif sorted(f2) == [0, 1]:
            if ref == [0, 0] and other == [1, 1]:
                return 'Mixed'  # Clear case of mixed ancestry
            elif ref == [0, 1] and other == [0, 0]:
                return 'RefFounder'  # Heterozygosity from 664c
            elif ref == [0, 0] and other == [0, 1]:
                return 'OtherFounder'  # Heterozygosity from other founder
            elif ref == [0, 1] and other == [0, 1]:
                return 'Both'  # Could be from either
            else:
                return 'Novel'  # Unexpected combination
        
        # Homozygous alternate in F2 (1/1)
        elif f2 == [1, 1]:
            if ref != [1, 1] and other == [1, 1]:
                return 'OtherFounder'  # From other founder
            elif ref == [1, 1] and other != [1, 1]:
                return 'RefFounder'  # From 664c
            elif ref == [1, 1] and other == [1, 1]:
                return 'Both'  # Uninformative
            else:
                return 'Novel'  # Unexpected combination
        
        # Handle other cases (multiallelic variants or parsing errors)
        else:
            return 'Complex'

def determine_ancestry(vcf_data, relationships, sample_col, ref_founder='664c'):
    """
    Determine the ancestry of each variant in an F2 individual,
    accounting for the reference founder (664c) in the analysis.
    
    Args:
        vcf_data (pd.DataFrame): Processed VCF data
        relationships (pd.DataFrame): Relationship map
        sample_col (str): Column name for the F2 sample
        
    Returns:
        pd.DataFrame: DataFrame with ancestry assignments
    """
    # Get the founders for this sample
    sample_id = sample_col.split('_')[0]  # Extract sample ID from column name
    sample_info = relationships[relationships['Sample'] == sample_id]
    
    if len(sample_info) == 0:
        raise ValueError(f"Sample {sample_id} not found in relationship map")
    
    founder1 = sample_info['Founder1'].iloc[0]
    founder2 = sample_info['Founder2'].iloc[0]
    
    # Identify which founder is the reference (664c)
    if founder1 == ref_founder:
        ref_founder_col = f"{founder1}_GT"
        other_founder_col = f"{founder2}_GT"
        ref_founder_name = "Founder1"
        other_founder_name = "Founder2"
    elif founder2 == ref_founder:
        ref_founder_col = f"{founder2}_GT"
        other_founder_col = f"{founder1}_GT"
        ref_founder_name = "Founder2"
        other_founder_name = "Founder1"
    else:
        raise ValueError(f"Reference founder 664c not found in relationship for {sample_id}")
    
    # Create a new DataFrame for ancestry results
    results = pd.DataFrame({
        'CHROM': vcf_data['CHROM'],
        'POS': vcf_data['POS'],
        'REF': vcf_data['REF'],
        'ALT': vcf_data['ALT'],
        'F2': vcf_data[sample_col],
        'RefFounder': vcf_data[ref_founder_col],
        'OtherFounder': vcf_data[other_founder_col]
    })
    
    # Function to parse a genotype into a list of alleles (0s and 1s)
    def parse_genotype(gt):
        if pd.isna(gt) or gt == '.':
            return None
            
        # Extract just the GT part if it contains colons
        if ':' in gt:
            gt = gt.split(':')[0]
            
        # Replace separators and extract digits
        gt = gt.replace('|', '/').replace('\\', '/')
        
        # If any allele is missing (.), return None
        if '.' in gt:
            return None
            
        # Extract alleles as integers (0, 1, etc.)
        try:
            alleles = [int(a) for a in gt if a.isdigit()]
            return alleles if alleles else None
        except:
            return None
    
    # Parse genotypes for all samples
    results['F2_alleles'] = results['F2'].apply(parse_genotype)
    results['RefFounder_alleles'] = results['RefFounder'].apply(parse_genotype)
    results['OtherFounder_alleles'] = results['OtherFounder'].apply(parse_genotype)
    
    # Determine ancestry based on allele combinations
    results['Ancestry'] = results.apply(assign_ancestry, axis=1)
    
    # Replace generic ancestry labels with specific founder names
    ancestry_mapping = {
        'RefFounder': ref_founder,
        'OtherFounder': other_founder_name
    }
    
    results['Ancestry'] = results['Ancestry'].replace(ancestry_mapping)
    
    # Clean up temporary columns
    results = results.drop(['F2_alleles', 'RefFounder_alleles', 'OtherFounder_alleles'], axis=1)
    
    # Rename columns for clarity in the output
    results = results.rename(columns={
        'RefFounder': ref_founder,
        'OtherFounder': founder1 if founder2 == ref_founder else founder2
    })
    
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
    
    if len(chromosomes) == 0:
        print(f"Warning: No data to plot for sample {sample_id}")
        return
    
    # Adjust figure size based on number of chromosomes
    plt.figure(figsize=(15, max(5, 2 * len(chromosomes))))
    
    # Color mapping
    colors = {'Founder1': 'blue', 'Founder2': 'red', 'Both': 'purple', 'Novel': 'green', 'Missing': 'gray'}
    
    for i, chrom in enumerate(chromosomes):
        chrom_data = ancestry_data[ancestry_data['CHROM'] == chrom]
        
        # Skip if no data for this chromosome
        if len(chrom_data) == 0:
            continue
            
        # Plot positions
        plt.subplot(len(chromosomes), 1, i+1)
        
        for ancestry, color in colors.items():
            subset = chrom_data[chrom_data['Ancestry'] == ancestry]
            if len(subset) > 0:  # Only plot if we have data for this ancestry
                plt.scatter(subset['POS'], [0.5] * len(subset), c=color, s=5, label=ancestry if i == 0 else "")
        
        plt.title(f'Chromosome {chrom}')
        plt.ylabel('Ancestry')
        
        # Set x-axis limits based on data or specified range
        min_pos = min(chrom_data['POS'])
        max_pos = max(chrom_data['POS'])
        plt.xlim(min_pos, max_pos)
        plt.yticks([])
        
        # Add position markers
        plt.xlabel(f'Position (Range: {min_pos:,} - {max_pos:,})')
    
    # Add legend to the first subplot
    if len(chromosomes) > 0:
        plt.subplot(len(chromosomes), 1, 1)
        plt.legend(loc='upper right')
    
    plt.tight_layout()
    
    # Add region info to filename if analyzing a specific region
    filename = f'{sample_id}_ancestry'
    if len(chromosomes) == 1:
        chrom = chromosomes[0]
        min_pos = min(ancestry_data['POS'])
        max_pos = max(ancestry_data['POS'])
        filename += f'_{chrom}_{min_pos}_{max_pos}'
    
    # Save plot and data
    plt.savefig(os.path.join(output_dir, f'{filename}.png'))
    plt.close()
    
    # Also save the data to a CSV file
    ancestry_data.to_csv(os.path.join(output_dir, f'{filename}.csv'), index=False)
