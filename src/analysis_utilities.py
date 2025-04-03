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
