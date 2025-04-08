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
import logging

# Setting up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

logger = logging.getLogger(__name__)

def normalize_genotype(gt):    
        """
        Normalizes a genotype string by handling missing values, extracting the genotype part,
        and ensuring consistent formatting for comparison.
        Args:
            gt (str): The genotype string to normalize. It may include phasing, colons, or missing values.
        Returns:
            str or None: The normalized genotype string with consistent formatting, or None if the input
            is missing or invalid.
        Notes:
            - Missing values (e.g., '.', './.', '.|.') are converted to None.
            - If the genotype contains additional information separated by colons, only the part before
                the first colon is retained.
            - Phasing indicators (e.g., '|', '\\') are replaced with '/' for consistency.
            - Genotypes containing missing alleles (e.g., '.') are also converted to None.
        
        """
        
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
    """
    Determines the ancestry of an individual based on allele information from F2 and founders.
    Uses an adaptive approach that doesn't assume specific allele combinations.
    
    Args:
        row (dict): Dictionary with F2_alleles, RefFounder_alleles, and AltFounder_alleles
        
    Returns:
        str: Ancestry classification
    """
    f2 = row['F2_alleles']
    ref = row['RefFounder_alleles']
    alt = row['AltFounder_alleles']
    
    # Handle missing data
    if f2 is None:
        return 'Missing'
    if ref is None or alt is None:
        return 'Missing_Founder'
        
    # For simplicity, sort the alleles to handle phasing differences
    f2 = sorted(f2)
    ref = sorted(ref)
    alt = sorted(alt)
    
    # Check if founders have different genotypes (informative variant)
    is_informative = ref != alt
    
    # Handle F2 homozygous cases
    if f2[0] == f2[1]:  # F2 is homozygous
        if f2 == ref and f2 != alt and is_informative:
            return 'RefFounder'  # F2 matches reference founder only
        elif f2 != ref and f2 == alt and is_informative:
            return 'AltFounder'  # F2 matches alternative founder only
        elif f2 == ref and f2 == alt:
            return 'Unresolved'  # Both founders have the same genotype
        else:
            return 'Novel'  # F2 doesn't match either founder
    
    # Handle F2 heterozygous cases
    elif f2[0] != f2[1]:  # F2 is heterozygous
        # Check if F2's heterozygosity could come from the two founders
        # This checks if one allele matches ref and the other matches alt
        ref_set = set(ref)
        alt_set = set(alt)
        f2_set = set(f2)
        
        # Perfect mixed ancestry: each allele from a different founder
        if len(ref_set) == 1 and len(alt_set) == 1 and ref_set != alt_set:
            if f2_set == ref_set.union(alt_set):
                return 'Mixed (het)'
            
        # Check if F2 is matching one founder exactly
        if f2 == ref and f2 != alt:
            return 'RefFounder'
        elif f2 != ref and f2 == alt:
            return 'AltFounder'
        elif f2 == ref and f2 == alt:
            return 'Unresolved'
        
        # Novel heterozygous combinations
        return 'Novel'
    
    # Handle other unexpected cases (should be rare)
    else:
        return 'Complex'

# Function to parse a genotype into a list of alleles (0s and 1s)
def parse_genotype(gt):
    """
    Parses a genotype string and extracts allele information.
    This function processes a genotype string (commonly found in VCF files) 
    and returns a list of alleles as integers. It handles missing data, 
    separators, and ensures the output is in a consistent format.
    Args:
        gt (str): The genotype string to parse. It may include alleles, 
                    separators ('|', '/', '\\'), and additional information 
                    separated by colons.
    Returns:
        list[int] or None: A list of alleles as integers if parsing is 
                            successful, or None if the input is invalid 
                            or contains missing data.
    Notes:
        - If the genotype string is missing (NaN or '.'), the function 
            returns None.
        - If the genotype string contains colons, only the part before 
            the first colon is considered.
        - Separators ('|', '\\') are replaced with '/' for uniformity.
        - If any allele is missing ('.'), the function returns None.
        - Alleles are extracted as integers. If no valid alleles are 
            found, the function returns None.

    """
    
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
        alt_founder_col = f"{founder2}_GT"
        other_founder_name = "Alternative"
    elif founder2 == ref_founder:
        ref_founder_col = f"{founder2}_GT"
        alt_founder_col = f"{founder1}_GT"
        other_founder_name = "Alternative"
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
        'AltFounder': vcf_data[alt_founder_col]
    })
    
    # Parse genotypes for all samples
    results['F2_alleles'] = results['F2'].apply(parse_genotype)
    results['RefFounder_alleles'] = results['RefFounder'].apply(parse_genotype)
    results['AltFounder_alleles'] = results['AltFounder'].apply(parse_genotype)
    
    # Determine ancestry based on allele combinations
    results['Ancestry'] = results.apply(assign_ancestry, axis=1)
    
    # Replace generic ancestry labels with specific founder names
    ancestry_mapping = {
        'RefFounder': ref_founder,
        'AltFounder': other_founder_name
    }
    
    results['Ancestry'] = results['Ancestry'].replace(ancestry_mapping)
    
    # Clean up temporary columns
    results = results.drop(['F2_alleles', 'RefFounder_alleles', 'AltFounder_alleles'], axis=1)
    
    # Rename columns for clarity in the output
    results = results.rename(columns={
        'RefFounder': ref_founder,
        'AltFounder': founder1 if founder2 == ref_founder else founder2
    })
    
    return results

def plot_ancestry(ancestry_results, sample_id, output_dir, window_size=1000000):
    """
    Create ancestry plots for a sample using stacked bar plots.
    
    Args:
        ancestry_results (pd.DataFrame): DataFrame with ancestry analysis results
        sample_id (str): ID of the sample to plot
        output_dir (str): Directory to save plots
        window_size (int): Size of the genomic window for aggregation (default: 1 Mb)
    """
    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{sample_id}_ancestry_plots.pdf")
    
    # Define colors for ancestry categories
    ancestry_colors = {
        'Missing': '#e0e0e0',       # Light gray
        'Unresolved': '#a9a9a9',    # Dark gray
        'Alternative': '#1f77b4',   # Blue
        'Novel': '#ff7f0e',         # Orange  
        '664c': '#2ca02c',          # Green (use the actual reference founder name)
        'Mixed (het)': '#9467bd'    # Purple
    }
    
    # Group data by chromosome and window
    ancestry_results['Window'] = (ancestry_results['POS'] // window_size) * window_size
    grouped = ancestry_results.groupby(['CHROM', 'Window', 'Ancestry']).size().reset_index(name='Count')
    
    # Normalize counts to proportions within each window
    grouped['Proportion'] = grouped.groupby(['CHROM', 'Window'])['Count'].transform(lambda x: x / x.sum())
    
    # Get unique chromosomes
    chromosomes = sorted(ancestry_results['CHROM'].unique())
    
    # Create PDF for plots
    with PdfPages(output_file) as pdf:
        for chrom in chromosomes:
            chrom_data = grouped[grouped['CHROM'] == chrom]
            if chrom_data.empty:
                continue
            
            # Pivot data for stacked bar plot
            pivot_data = chrom_data.pivot(index='Window', columns='Ancestry', values='Proportion').fillna(0)
            pivot_data = pivot_data[list(ancestry_colors.keys())]  # Ensure consistent order of categories
            
            # Create stacked bar plot
            fig, ax = plt.subplots(figsize=(12, 4))
            bottom = None
            for ancestry, color in ancestry_colors.items():
                if ancestry in pivot_data:
                    ax.bar(
                        pivot_data.index, pivot_data[ancestry], width=window_size, 
                        bottom=bottom, color=color, label=ancestry, align='edge'
                    )
                    bottom = pivot_data[ancestry] if bottom is None else bottom + pivot_data[ancestry]
            
            # Configure plot
            ax.set_title(f'Chromosome {chrom} - Sample {sample_id}')
            ax.set_xlabel('Genomic Position (bp)')
            ax.set_ylabel('Proportion')
            ax.set_xlim(0, ancestry_results[ancestry_results['CHROM'] == chrom]['POS'].max())
            ax.set_ylim(0, 1)
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False)
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            
            # Save plot to PDF
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()
    
    logger.info(f"Created ancestry plots for {sample_id} in {output_file}")