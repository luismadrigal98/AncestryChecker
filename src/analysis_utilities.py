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
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

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

def determine_ancestry(vcf_data, relationships, sample_col, ref_founder=None):
    """
    Determine the ancestry of each variant in an F2 individual.
    
    Args:
        vcf_data (pd.DataFrame): Processed VCF data
        relationships (pd.DataFrame): Relationship map
        sample_col (str): Column name for the F2 sample
        ref_founder (str, optional): Name of reference founder. If None, 
                                    the first founder will be used.
        
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
    
    # If no reference founder is specified, use the first founder
    if ref_founder is None:
        ref_founder = founder1
        print(f"No reference founder specified. Using {founder1} as reference.")
    
    # Identify which founder is the reference
    if founder1 == ref_founder:
        ref_founder_col = f"{founder1}_GT"
        other_founder_col = f"{founder2}_GT"
        other_founder_name = founder2
    elif founder2 == ref_founder:
        ref_founder_col = f"{founder2}_GT"
        other_founder_col = f"{founder1}_GT"
        other_founder_name = founder1
    else:
        raise ValueError(f"Reference founder {ref_founder} not found in relationship for {sample_id}. "
                         f"Available founders: {founder1}, {founder2}")
    
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
                return 'RefFounder'
            elif other == [0, 0] and ref != [0, 0]:
                return 'OtherFounder'
            elif ref == [0, 0] and other == [0, 0]:
                return 'Both'
            else:
                return 'Novel'
        
        # Heterozygous in F2 (0/1)
        elif sorted(f2) == [0, 1]:
            if ref == [0, 0] and other == [1, 1]:
                return 'Mixed'
            elif ref == [0, 1] and other == [0, 0]:
                return 'RefFounder'
            elif ref == [0, 0] and other == [0, 1]:
                return 'OtherFounder'
            elif ref == [0, 1] and other == [0, 1]:
                return 'Both'
            else:
                return 'Novel'
        
        # Homozygous alternate in F2 (1/1)
        elif f2 == [1, 1]:
            if ref != [1, 1] and other == [1, 1]:
                return 'OtherFounder'
            elif ref == [1, 1] and other != [1, 1]:
                return 'RefFounder'
            elif ref == [1, 1] and other == [1, 1]:
                return 'Both'
            else:
                return 'Novel'
        
        # Handle other cases (multiallelic variants or parsing errors)
        else:
            return 'Complex'
    
    results['Ancestry'] = results.apply(assign_ancestry, axis=1)
    
    # Create human-readable ancestry labels
    ancestry_mapping = {
        'RefFounder': ref_founder,
        'OtherFounder': other_founder_name,
        'Both': 'Both',
        'Mixed': 'Mixed',
        'Novel': 'Novel',
        'Complex': 'Complex',
        'Missing': 'Missing',
        'Missing_Founder': 'Missing_Founder'
    }
    
    # Apply the mapping but keep descriptive labels for non-founder specific categories
    results['Ancestry'] = results['Ancestry'].map(lambda x: ancestry_mapping.get(x, x))
    
    # Clean up temporary columns
    results = results.drop(['F2_alleles', 'RefFounder_alleles', 'OtherFounder_alleles'], axis=1)
    
    # Rename columns for clarity in the output
    results = results.rename(columns={
        'RefFounder': ref_founder,
        'OtherFounder': other_founder_name
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
    
    # Save the complete data to a CSV file
    ancestry_data.to_csv(os.path.join(output_dir, f'{sample_id}_ancestry_full.csv'), index=False)
    
    # Get unique chromosomes
    chromosomes = ancestry_data['CHROM'].unique()
    
    if len(chromosomes) == 0:
        print(f"Warning: No data to plot for sample {sample_id}")
        return
    
    # Create a multi-page PDF for all chromosomes
    from matplotlib.backends.backend_pdf import PdfPages
    pdf_path = os.path.join(output_dir, f'{sample_id}_ancestry_plots.pdf')
    
    # Create a color map for all ancestry types
    unique_ancestry = ancestry_data['Ancestry'].unique()
    colors = {
        'Both': 'purple',
        'Mixed': 'green',
        'Novel': 'orange',
        'Missing': 'lightgray',
        'Complex': 'black',
        'Missing_Founder': 'darkgray'
    }
    
    # Add colors for founder-specific ancestry
    # Use different shades of blue/red to distinguish between different founders
    founder_names = [a for a in unique_ancestry if a not in colors]
    blues = plt.cm.Blues(np.linspace(0.6, 1.0, (len(founder_names) + 1) // 2))
    reds = plt.cm.Reds(np.linspace(0.6, 1.0, len(founder_names) // 2 + 1))
    
    for i, founder in enumerate(founder_names):
        if i % 2 == 0:  # Alternate between blues and reds
            colors[founder] = blues[i // 2]
        else:
            colors[founder] = reds[i // 2]
    
    with PdfPages(pdf_path) as pdf:
        # Create a summary page
        plt.figure(figsize=(8, 6))
        plt.axis('off')
        plt.text(0.1, 0.9, f"Ancestry Analysis for Sample: {sample_id}", fontsize=16)
        plt.text(0.1, 0.8, f"Total variants analyzed: {len(ancestry_data)}", fontsize=12)
        
        # Add ancestry distribution information
        ancestry_counts = ancestry_data['Ancestry'].value_counts()
        plt.text(0.1, 0.7, "Ancestry Distribution:", fontsize=14)
        y_pos = 0.65
        for ancestry, count in ancestry_counts.items():
            percentage = 100 * count / len(ancestry_data)
            color = colors.get(ancestry, 'black')
            if isinstance(color, (list, np.ndarray)):  # Handle RGB arrays
                plt.text(0.1, y_pos, f"{ancestry}: {count} variants ({percentage:.2f}%)", 
                        fontsize=10, color=color)
            else:
                plt.text(0.1, y_pos, f"{ancestry}: {count} variants ({percentage:.2f}%)", 
                        fontsize=10, color=color)
            y_pos -= 0.05
        
        pdf.savefig()
        plt.close()
        
        # Then create individual plots for each chromosome
        for chrom in chromosomes:
            chrom_data = ancestry_data[ancestry_data['CHROM'] == chrom]
            
            if len(chrom_data) == 0:
                continue
                
            plt.figure(figsize=(12, 6))
            
            # Plot positions for each ancestry type
            for ancestry in unique_ancestry:
                if ancestry in chrom_data['Ancestry'].values:
                    subset = chrom_data[chrom_data['Ancestry'] == ancestry]
                    color = colors.get(ancestry, 'gray')
                    plt.scatter(subset['POS'], [0.5] * len(subset), 
                                c=color, s=5, label=ancestry)
            
            # Format x-axis to show positions in Mb
            def format_mb(x, pos):
                return f'{x/1e6:.1f}'
                
            from matplotlib.ticker import FuncFormatter
            plt.gca().xaxis.set_major_formatter(FuncFormatter(format_mb))
            
            # Set plot title and labels
            plt.title(f'Chromosome {chrom} - Sample {sample_id}')
            plt.xlabel('Position (Mb)')
            plt.yticks([])
            
            # Add region markers
            min_pos = min(chrom_data['POS'])
            max_pos = max(chrom_data['POS'])
            plt.xlim(min_pos, max_pos)
            
            # Add legend OUTSIDE the plot area
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            # Save to PDF
            pdf.savefig(bbox_inches='tight')
            
            # Also save individual PNGs for easier viewing
            plt.savefig(os.path.join(output_dir, f'{sample_id}_chr{chrom}_ancestry.png'), 
                        bbox_inches='tight')
            plt.close()
        
        # Create a chromosome overview plot
        plt.figure(figsize=(15, 10))
        
        # Try to sort chromosomes numerically if possible
        try:
            chrom_numeric = [int(str(c).lower().replace('chr', '')) for c in chromosomes]
            sorted_indices = sorted(range(len(chrom_numeric)), key=lambda k: chrom_numeric[k])
            sorted_chroms = [chromosomes[i] for i in sorted_indices]
        except:
            # If not numeric, just use the original order
            sorted_chroms = chromosomes
            
        # Create a summary plot showing ancestry distribution across all chromosomes
        for i, chrom in enumerate(sorted_chroms):
            chrom_data = ancestry_data[ancestry_data['CHROM'] == chrom]
            
            if len(chrom_data) == 0:
                continue
            
            # Normalize positions to 0-1 range for consistent width
            min_pos = min(chrom_data['POS'])
            max_pos = max(chrom_data['POS'])
            range_pos = max_pos - min_pos
            
            # Plot each chromosome as a row
            y_pos = len(sorted_chroms) - i
            
            # Plot background for chromosome
            plt.plot([0, 1], [y_pos, y_pos], 'k-', lw=1)
            
            # Add markers for each variant, colored by ancestry
            for ancestry in unique_ancestry:
                subset = chrom_data[chrom_data['Ancestry'] == ancestry]
                if len(subset) > 0:
                    # Normalize positions
                    norm_pos = [(p - min_pos)/range_pos for p in subset['POS']]
                    plt.scatter(norm_pos, [y_pos] * len(subset), 
                                c=colors.get(ancestry, 'gray'), 
                                s=3, alpha=0.7)
        
        plt.yticks(range(1, len(sorted_chroms)+1), sorted_chroms)
        plt.title(f'Ancestry Overview - Sample {sample_id}')
        plt.xlabel('Normalized Position')
        plt.ylabel('Chromosome')
        
        # Create a custom legend
        legend_elements = []
        for ancestry in unique_ancestry:
            color = colors.get(ancestry, 'gray')
            if isinstance(color, (list, np.ndarray)):
                color = tuple(color)
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                            label=ancestry, markerfacecolor=color, markersize=8))
                           
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        pdf.savefig(bbox_inches='tight')
        
        # Save the overview as a separate image too
        plt.savefig(os.path.join(output_dir, f'{sample_id}_ancestry_overview.png'), 
                    bbox_inches='tight')
        plt.close()
    
    print(f"Created ancestry plots for {sample_id} in {pdf_path}")