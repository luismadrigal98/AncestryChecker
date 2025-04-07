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

def plot_ancestry(ancestry_results, sample_id, output_dir):
    """
    Create improved ancestry plots for a sample.
    
    Args:
        ancestry_results (pd.DataFrame): DataFrame with ancestry analysis results
        sample_id (str): ID of the sample to plot
        output_dir (str): Directory to save plots
    """
    import os
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{sample_id}_ancestry_plots.pdf")
    
    # Get unique chromosomes and sort them
    chromosomes = sorted(ancestry_results['CHROM'].unique())
    
    # Define colors for ancestry categories - matching your existing names
    ancestry_colors = {
        'Missing': '#e0e0e0',       # Light gray
        'Unresolved': '#a9a9a9',    # Dark gray
        'Alternative': '#1f77b4',   # Blue
        'Novel': '#ff7f0e',         # Orange  
        '664c': '#2ca02c',          # Green (use the actual reference founder name)
        'Mixed (het)': '#9467bd'    # Purple
    }
    
    # Group data by chromosome for plotting
    with PdfPages(output_file) as pdf:
        # Create a summary plot with all chromosomes
        fig, axes = plt.subplots(len(chromosomes), 1, 
                                figsize=(10, max(3, len(chromosomes) * 0.7)))
        
        if len(chromosomes) == 1:
            axes = [axes]  # Make it iterable if there's only one chromosome
            
        for i, chrom in enumerate(chromosomes):
            chrom_data = ancestry_results[ancestry_results['CHROM'] == chrom]
            
            # Skip if no data for this chromosome
            if len(chrom_data) == 0:
                continue
                
            # Get axis for this chromosome
            ax = axes[i]
            
            # Add light background
            max_pos = chrom_data['POS'].max() if len(chrom_data) > 0 else 1000
            ax.axhspan(0.3, 0.7, color='#f8f8f8', alpha=0.7)
            
            # Plot each ancestry category
            for ancestry in ancestry_results['Ancestry'].unique():
                subset = chrom_data[chrom_data['Ancestry'] == ancestry]
                if len(subset) == 0:
                    continue
                
                # Get color, default to gray if not in our color map
                color = ancestry_colors.get(ancestry, '#cccccc')
                
                # Add slight y-jitter to avoid perfect overlaps
                y_pos = 0.5 + (np.random.rand(len(subset)) - 0.5) * 0.1
                
                # Plot with explicit color to avoid the warning
                ax.scatter(subset['POS'], y_pos, color=color, s=10, 
                           alpha=0.7, label=ancestry if i == 0 else "")
            
            # Configure axis
            ax.set_xlim(0, max_pos * 1.05)
            ax.set_ylim(0.2, 0.8)
            ax.set_yticks([])
            
            if i == len(chromosomes) - 1:
                ax.set_xlabel('Position (bp)')
            else:
                ax.set_xticklabels([])
                
            ax.set_ylabel(f'Chr {chrom}')
            ax.grid(True, linestyle='--', alpha=0.4)
        
        # Add single legend at the bottom
        handles = []
        labels = []
        for ancestry in sorted(ancestry_results['Ancestry'].unique()):
            color = ancestry_colors.get(ancestry, '#cccccc')
            handles.append(mpatches.Patch(color=color, label=ancestry))
            labels.append(ancestry)
            
        fig.legend(handles=handles, labels=labels, 
                  loc='lower center', ncol=min(3, len(labels)),
                  bbox_to_anchor=(0.5, 0.01), frameon=True)
        
        # Add title
        plt.suptitle(f'Ancestry Distribution for Sample {sample_id}', 
                    fontsize=14, y=0.98)
        
        # Adjust layout to make room for the legend
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)
        
        # Save to PDF
        pdf.savefig(fig)
        plt.close()
        
        # Generate individual chromosome plots for more detail
        for chrom in chromosomes:
            chrom_data = ancestry_results[ancestry_results['CHROM'] == chrom]
            if len(chrom_data) == 0:
                continue
                
            fig, ax = plt.subplots(figsize=(10, 3))
            max_pos = chrom_data['POS'].max()
            
            # Add background
            ax.axhspan(0.3, 0.7, color='#f8f8f8', alpha=0.7)
            
            # Plot each ancestry category
            for ancestry in ancestry_results['Ancestry'].unique():
                subset = chrom_data[chrom_data['Ancestry'] == ancestry]
                if len(subset) == 0:
                    continue
                    
                color = ancestry_colors.get(ancestry, '#cccccc')
                y_pos = 0.5 + (np.random.rand(len(subset)) - 0.5) * 0.1
                
                ax.scatter(subset['POS'], y_pos, color=color, s=15,
                          alpha=0.7, label=ancestry)
            
            # Configure axis
            ax.set_xlim(0, max_pos * 1.05)
            ax.set_ylim(0.2, 0.8)
            ax.set_yticks([])
            ax.set_xlabel('Position (bp)')
            ax.set_title(f'Chromosome {chrom} - Sample {sample_id}')
            ax.grid(True, linestyle='--', alpha=0.4)
            
            # Add legend
            ax.legend(loc='upper right', frameon=True)
            
            # Save and close
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()
    
    logger.info(f"Created ancestry plots for {sample_id} in {output_file}")
    
    # Optional: Create interactive HTML plot if plotly is available
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        html_file = os.path.join(output_dir, f"{sample_id}_ancestry_plots.html")
        
        fig = make_subplots(rows=len(chromosomes), cols=1,
                         subplot_titles=[f"Chromosome {chrom}" for chrom in chromosomes],
                         vertical_spacing=0.05)
        
        for i, chrom in enumerate(chromosomes):
            chrom_data = ancestry_results[ancestry_results['CHROM'] == chrom]
            
            # Plot each ancestry type
            for ancestry in ancestry_results['Ancestry'].unique():
                subset = chrom_data[chrom_data['Ancestry'] == ancestry]
                if len(subset) == 0:
                    continue
                
                color = ancestry_colors.get(ancestry, '#cccccc')
                
                fig.add_trace(
                    go.Scatter(
                        x=subset['POS'],
                        y=[0.5] * len(subset),
                        mode='markers',
                        name=ancestry,
                        marker=dict(color=color, size=6),
                        legendgroup=ancestry,
                        showlegend=(i == 0)  # Only show in legend once
                    ),
                    row=i+1, col=1
                )
            
            # Update axis properties
            fig.update_yaxes(range=[0, 1], showticklabels=False, 
                           title_text="", row=i+1, col=1)
            
            if i < len(chromosomes) - 1:
                fig.update_xaxes(showticklabels=False, row=i+1, col=1)
            else:
                fig.update_xaxes(title_text="Position (bp)", row=i+1, col=1)
        
        # Update layout
        fig.update_layout(
            title=f"Ancestry Distribution for Sample {sample_id}",
            height=max(500, 200 * len(chromosomes)),
            width=900,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            )
        )
        
        # Save HTML
        fig.write_html(html_file)
        logger.info(f"Created interactive HTML ancestry plot at {html_file}")
        
    except ImportError:
        logger.info("Plotly not available - skipping interactive plot")
    
    return output_file