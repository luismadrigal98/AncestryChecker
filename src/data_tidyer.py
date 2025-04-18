"""
Data preprocessing utilities for AncestryChecker.
This module provides functions to clean and preprocess VCF data for ancestry analysis.

@module: data_tidyer
@description: Functions for cleaning and preprocessing genomic data
@date: 2025-04-01
@version: 1.0
@license: MIT
"""

import pandas as pd
import numpy as np
import logging

# Setting up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger = logging.getLogger(__name__)

# Extract just the chromosome number (removing any 'Chr_', 'chr', etc. prefixes)
def normalize_chrom(x):
    # Convert to lowercase string
    x = str(x).lower()
    # Remove common prefixes
    for prefix in ['chr_', 'chr', 'chromosome_', 'chromosome']:
        if x.startswith(prefix):
            x = x[len(prefix):]
    # Handle leading zeros (e.g., '01' -> '1')
    try:
        return str(int(x))
    except ValueError:
        # For non-numeric chromosomes (e.g., 'X', 'Y', 'MT')
        return x

def filter_founder_homozygous(vcf_df, founders):
    """
    This function will filter out variants that are not homozygous in all founders.
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        founders (list): List of founder sample names
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    # Find GT columns for founders (they have _GT suffix after filter_vcf_data)
    founder_gt_cols = [f"{founder}_GT" for founder in founders]
    
    # Check if any founders were found
    if not all(col in vcf_df.columns for col in founder_gt_cols):
        logger.warning(f"Not all founder GT columns found in the data")
        available = [col for col in founder_gt_cols if col in vcf_df.columns]
        logger.warning(f"Available founder GT columns: {', '.join(available)}")
        founder_gt_cols = available  # Use only available columns
    
    # Check if SNPs are homozygous in all founders (any homozygous genotype)
    def is_homozygous(row):
        genotypes = [row[col] for col in founder_gt_cols if not pd.isna(row[col])]
        if not genotypes:  # If all are missing
            return False
        
        return all(gt.split('/')[0] == gt.split('/')[-1] if '/' in gt else 
                gt.split('|')[0] == gt.split('|')[-1] if '|' in gt else False 
                for gt in genotypes)
        
    homozygous = vcf_df.apply(is_homozygous, axis=1)

    filtered = vcf_df[homozygous]
    logger.info(f"Homozygosity filter: {len(filtered)}/{len(vcf_df)} variants retained ({len(filtered)/len(vcf_df)*100:.1f}%)")

    return vcf_df[homozygous]

def filter_by_region(vcf_df, chrom=None, start_pos=None, end_pos=None):
    """
    Filter VCF data by genomic region.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        chrom (str): Chromosome to filter for (e.g., "1", "Chr1")
        start_pos (int): Start position for region filtering
        end_pos (int): End position for region filtering
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    filtered_df = vcf_df.copy()
    
    if chrom is not None:
        # Convert chrom to string if it's not already
        chrom = str(chrom)
                
        # Get normalized version of the target chromosome
        target_chrom = normalize_chrom(chrom)
        
        # Apply the same normalization to the chromosome column and compare
        chrom_matches = filtered_df['CHROM'].apply(normalize_chrom) == target_chrom
        filtered_df = filtered_df[chrom_matches]
        
        # Print debug info
        if len(filtered_df) == 0:
            logger.warning(f"Warning: No variants found for chromosome '{chrom}'")
            logger.info(f"Available chromosomes: {', '.join(filtered_df['CHROM'].unique())}")

    if start_pos is not None:
        filtered_df = filtered_df[filtered_df['POS'] >= start_pos]
        
    if end_pos is not None:
        filtered_df = filtered_df[filtered_df['POS'] <= end_pos]
    
    return filtered_df

def read_relationship_map(path):
    """
    Read and parse a relationship map file.
    
    Args:
        path (str): Path to the relationship map file
        
    Returns:
        pd.DataFrame: DataFrame containing the relationship information
    """
    # Explicitly set dtype=str for all columns to prevent numeric interpretation
    return pd.read_csv(path, sep='\t', header=None, 
                        names=['Sample', 'Generation', 'Founder1', 'Founder2'],
                        dtype=str)  # Force all columns to be read as strings

def extract_gt(value, vcf_df):
        """
        Extracts the genotype (GT) field from a VCF format string.
        This function parses a colon-delimited string from a VCF file and extracts 
        the genotype information based on the FORMAT specification.
        Parameters
        ----------
        value : str or float
            The VCF format string to parse, typically from a sample column.
            Can be NaN, in which case NaN is returned.
        vcf_df : pd.DataFrame
            DataFrame containing the VCF data, used to determine the index of the GT field.
        Returns
        -------
        str or numpy.nan
            The extracted genotype value if available, or numpy.nan if the input is NaN
            or if the GT field cannot be extracted.
        Notes
        -----
        The function depends on an external variable `filtered_vcf` that contains 
        the FORMAT column from the VCF file.
        """

        if pd.isna(value):
            return np.nan
        fields = value.split(':')
        format_fields = vcf_df['FORMAT'].iloc[0].split(':')
        gt_idx = format_fields.index('GT') if 'GT' in format_fields else 0
        return fields[gt_idx] if gt_idx < len(fields) else np.nan

def filter_vcf_data(vcf_df, founders, f2_samples, allow_missing=True):
    """
    Filter VCF data to keep only relevant SNPs for ancestry analysis.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        founders (list): List of founder sample names
        f2_samples (list): List of F2 sample names
        allow_missing (bool): Whether to allow missing data in F2 samples
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    # Extract genotype columns (all columns after FORMAT)
    format_col_idx = vcf_df.columns.get_loc('FORMAT')
    genotype_cols = vcf_df.columns[format_col_idx + 1:]
    
    # Find columns for founders and F2 samples
    founder_cols = [col for col in genotype_cols if col in founders]
    f2_cols = [col for col in genotype_cols if col in f2_samples]
    
    # Debug information
    missing_founders = [f for f in founders if f not in genotype_cols]
    missing_f2s = [f for f in f2_samples if f not in genotype_cols]
    
    if missing_founders:
        logger.warning(f"Warning: The following founders are not in the VCF: {', '.join(missing_founders)}")
        logger.info(f"Available samples in VCF: {', '.join(genotype_cols)}")
    
    if missing_f2s:
        logger.warning(f"Warning: The following F2 samples are not in the VCF: {', '.join(missing_f2s)}")
    
    # Ensure we have at least some founders and F2 samples
    if len(founder_cols) == 0:
        raise ValueError("No founders found in VCF file")
    if len(f2_cols) == 0:
        raise ValueError("No F2 samples found in VCF file")
    
    # Create a copy of the DataFrame to avoid modifying the original
    filtered_vcf = vcf_df.copy()
    
    # Apply the function to genotype columns
    for col in founder_cols + f2_cols:
        filtered_vcf[f'{col}_GT'] = filtered_vcf[col].apply(lambda x: extract_gt(x, filtered_vcf))
    
    # Filter rows with complete data for founders
    has_complete_founders = filtered_vcf[[f'{col}_GT' for col in founder_cols]].notna().all(axis=1)
    
    if allow_missing:
        # At least one F2 sample has data
        has_some_f2_data = filtered_vcf[[f'{col}_GT' for col in f2_cols]].notna().any(axis=1)
        filtered_vcf = filtered_vcf[has_complete_founders & has_some_f2_data]
    else:
        # Complete data for all samples
        has_complete_data = filtered_vcf[[f'{col}_GT' for col in founder_cols + f2_cols]].notna().all(axis=1)
        filtered_vcf = filtered_vcf[has_complete_data]
    
    return filtered_vcf

def identify_informative_snps(filtered_vcf, founder_cols):
    """
    Identify SNPs that are informative for ancestry determination.
    
    Args:
        filtered_vcf (pd.DataFrame): Filtered VCF data
        founder_cols (list): List of founder column names
        
    Returns:
        pd.DataFrame: VCF data with only informative SNPs
    """
    # SNPs are informative if they differ between founders
    gt_cols = [f'{col}_GT' for col in founder_cols]
    
    # Check if SNPs are different between founders
    def is_informative(row):
        genotypes = [row[col] for col in gt_cols if not pd.isna(row[col])]
        return len(set(genotypes)) > 1
    
    informative = filtered_vcf.apply(is_informative, axis=1)
    return filtered_vcf[informative]

def filter_by_allele_depth(vcf_df, sample_cols, founder_cols=None, min_depth=5, caller="freebayes"):
    """
    Filter VCF data by allele depth.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        founder_cols (list, optional): List of founder column names
        min_depth (int): Minimum depth threshold for filtering
        caller (str): Variant caller used (e.g., "freebayes", "bcftools", "gatk")
    
    Returns:
        pd.DataFrame: Filtered VCF data with sufficient read depth
    """
    if founder_cols is None:
        founder_cols = []
    
    all_cols = sample_cols + founder_cols
    filtered_df = vcf_df.copy()
    
    # Extract allele depth for each sample
    for col in all_cols:
        if col not in filtered_df.columns:
            logger.warning(f"Sample column {col} not found in VCF data")
            continue
            
        filtered_df[f'{col}_AD'] = filtered_df.apply(
            lambda row: extract_allele_depth(row[col], row['FORMAT'], caller), 
            axis=1
        )
    
    # Filter rows by minimum depth
    depth_cols = [f'{col}_AD' for col in all_cols if f'{col}_AD' in filtered_df.columns]
    if depth_cols:
        filtered_df = filtered_df[filtered_df[depth_cols].fillna(0).min(axis=1) >= min_depth]
        logger.info(f"Depth filter: {len(filtered_df)}/{len(vcf_df)} variants retained ({len(filtered_df)/len(vcf_df)*100:.1f}%)")
    else:
        logger.warning("No depth columns found for filtering")
    
    # Drop temporary depth columns
    filtered_df = filtered_df.drop(columns=[col for col in filtered_df.columns if col.endswith('_AD')], errors='ignore')
    
    return filtered_df

def extract_allele_depth(value, format_str, caller="freebayes"):
    """
    Extract the total allele depth from a VCF sample field.
    
    Args:
        value (str): Sample field value from VCF
        format_str (str): FORMAT field specifying the data structure
        caller (str): Variant caller used
    
    Returns:
        int or np.nan: Total allele depth or np.nan if not available
    """
    if pd.isna(value) or pd.isna(format_str):
        return np.nan
        
    format_fields = format_str.split(':')
    value_fields = value.split(':')
    
    # Handle different variant callers
    try:
        if caller == "freebayes":
            # In freebayes, can use DP field directly or sum RO+AO
            if 'DP' in format_fields and len(value_fields) > format_fields.index('DP'):
                dp_idx = format_fields.index('DP')
                return int(value_fields[dp_idx])
            elif all(field in format_fields for field in ['RO', 'AO']):
                ro_idx = format_fields.index('RO')
                ao_idx = format_fields.index('AO')
                if len(value_fields) > max(ro_idx, ao_idx):
                    return int(value_fields[ro_idx]) + int(value_fields[ao_idx])
        
        elif caller in ["bcftools", "gatk"]:
            # For bcftools/GATK, use AD field (comma-separated values)
            if 'AD' in format_fields:
                ad_idx = format_fields.index('AD')
                if len(value_fields) > ad_idx:
                    ad_values = value_fields[ad_idx].split(',')
                    return sum(int(ad) for ad in ad_values if ad.isdigit())
    except (ValueError, IndexError):
        pass
    
    return np.nan

def filter_biallelic_snps(vcf_df):
    """
    Filter VCF data to keep only biallelic SNPs.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        
    Returns:
        pd.DataFrame: Filtered VCF data with only biallelic SNPs
    """
    # Count comma-separated values in ALT field to determine number of alternate alleles
    is_biallelic = vcf_df['ALT'].apply(lambda x: ',' not in x)
    
    # Also filter out non-SNPs (where REF or ALT length > 1)
    is_snp = (vcf_df['REF'].str.len() == 1) & (vcf_df['ALT'].str.len() == 1)
    
    return vcf_df[is_biallelic & is_snp]

def calculate_maf(vcf_df, sample_cols, format_fields=['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL'], caller = "freebayes"):
    """
    Calculate minor allele frequency for each SNP.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        format_fields (list): List of format fields to extract from VCF data
        caller (str): Variant caller used (e.g., "freebayes", "gatk"). This is relevant because the vcf format may differ.
        
    Returns:
        pd.Series: Minor allele frequencies
    """
    # FORMAT: GT:DP:AD:RO:QR:AO:QA:GL

    # For each sample column, we need the observed counts of alleles, which are RO and AO
    # The frequency fo the minor allele is min{RO, AO} / (RO + AO)

    # Count alleles
    def get_allele_frequency(row, format_fields=format_fields):
        """Calculate MAF with error handling for truncated format strings"""
        for col in sample_cols:
            if pd.isna(row[col]):
                return pd.NA
                
            try:
                if caller == 'freebayes':
                
                    # Get the format parts
                    parts = row[col].split(":")
                    
                    # Make sure we have enough elements
                    if len(parts) <= max(format_fields.index('RO'), format_fields.index('AO')):
                        return pd.NA
                        
                    # Get the allele counts
                    RO_ix = format_fields.index('RO')
                    AO_ix = format_fields.index('AO')
                    RO = int(parts[RO_ix])
                    AO = int(parts[AO_ix])
                    
                    # Calculate the minor allele frequency
                    total = RO + AO
                    if total == 0:
                        return pd.NA
                    maf = min(RO, AO) / total
                    return maf
                
                elif caller == 'bcftools':
                    # Get format definition and verify AD field exists
                    format_str = row['FORMAT']
                    if 'AD' not in format_str.split(':'):
                        return pd.NA
                    
                    # Find position of AD field in the format string
                    format_parts = format_str.split(':')
                    ad_index = format_parts.index('AD')
                    
                    # Split sample data and get AD values
                    parts = row[col].split(':')
                    
                    # Check if we have enough parts to reach AD index
                    if len(parts) <= ad_index:
                        return pd.NA
                    
                    # Extract allele depths
                    ad_values = parts[ad_index].split(',')
                    if len(ad_values) < 2:
                        return pd.NA
                    
                    try:
                        RO = int(ad_values[0])  # Reference allele count
                        AO = int(ad_values[1])  # Alternate allele count
                        
                        # Calculate the minor allele frequency
                        total = RO + AO
                        if total == 0:
                            return pd.NA
                        maf = min(RO, AO) / total
                        return maf
                    except (ValueError, IndexError):
                        return pd.NA
                
                elif caller == 'gatk':
                    raise NotImplementedError("GATK format parsing not implemented yet.")
                else:
                    raise ValueError(f"Unknown variant caller: {caller}")
                
            except (IndexError, ValueError) as e:
                # Handle any parsing errors
                return pd.NA
        
        # In case no valid samples were found
        return pd.NA
    
    # Apply the function to each row
    return vcf_df.apply(get_allele_frequency, axis=1)

def filter_by_maf(vcf_df, sample_cols, min_maf=0.05, format_fields=['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL'], caller = "freebayes"):
    """
    Filter VCF data by minor allele frequency.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        min_maf (float): Minimum minor allele frequency threshold
        format_fields (list): List of format fields to extract from VCF data
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    maf_values = calculate_maf(vcf_df, sample_cols, format_fields, caller)
    return vcf_df[maf_values >= min_maf]

def filter_by_missing_rate(vcf_df, sample_cols, max_missing_rate=0.2):
    """
    Filter VCF data by missing rate.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        max_missing_rate (float): Maximum allowed rate of missing data
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    gt_cols = [f'{col}_GT' for col in sample_cols]
    missing_rates = vcf_df[gt_cols].isna().mean(axis=1)
    return vcf_df[missing_rates <= max_missing_rate]

def filter_by_qual(vcf_df, min_qual=30):
    """
    Filter VCF data by QUAL value.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        min_qual (float): Minimum QUAL value threshold
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    # Convert QUAL to float, handling any non-numeric values
    vcf_df['QUAL_numeric'] = pd.to_numeric(vcf_df['QUAL'], errors='coerce')
    filtered_df = vcf_df[vcf_df['QUAL_numeric'] >= min_qual]
    if 'QUAL_numeric' in filtered_df.columns:
        filtered_df = filtered_df.drop('QUAL_numeric', axis=1)
    return filtered_df