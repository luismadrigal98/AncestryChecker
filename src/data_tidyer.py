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
        print(f"Warning: The following founders are not in the VCF: {', '.join(missing_founders)}")
        print(f"Available samples in VCF: {', '.join(genotype_cols)}")
    
    if missing_f2s:
        print(f"Warning: The following F2 samples are not in the VCF: {', '.join(missing_f2s)}")
    
    # Ensure we have at least some founders and F2 samples
    if len(founder_cols) == 0:
        raise ValueError("No founders found in VCF file")
    if len(f2_cols) == 0:
        raise ValueError("No F2 samples found in VCF file")
    
    # Create a copy of the DataFrame to avoid modifying the original
    filtered_vcf = vcf_df.copy()
    
    # Extract genotype values (assuming FORMAT field contains 'GT')
    def extract_gt(value):
        if pd.isna(value):
            return np.nan
        fields = value.split(':')
        format_fields = filtered_vcf['FORMAT'].iloc[0].split(':')
        gt_idx = format_fields.index('GT') if 'GT' in format_fields else 0
        return fields[gt_idx] if gt_idx < len(fields) else np.nan
    
    # Apply the function to genotype columns
    for col in founder_cols + f2_cols:
        filtered_vcf[f'{col}_GT'] = filtered_vcf[col].apply(extract_gt)
    
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

def calculate_maf(vcf_df, sample_cols):
    """
    Calculate minor allele frequency for each SNP.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        
    Returns:
        pd.Series: Minor allele frequencies
    """
    # Extract genotype values
    gt_cols = [f'{col}_GT' for col in sample_cols]
    
    # Count alleles
    def count_alleles(row):
        allele_counts = {'0': 0, '1': 0}
        for col in gt_cols:
            if pd.isna(row[col]):
                continue
            
            # Handle different formats (0/0, 0|0, etc.)
            genotype = row[col].replace('|', '/').split('/')
            for allele in genotype:
                if allele in allele_counts:
                    allele_counts[allele] += 1
        
        # Calculate MAF
        total_alleles = sum(allele_counts.values())
        if total_alleles == 0:
            return 0
        
        min_count = min(allele_counts.values())
        return min_count / total_alleles if total_alleles > 0 else 0
    
    return vcf_df.apply(count_alleles, axis=1)

def filter_by_maf(vcf_df, sample_cols, min_maf=0.05):
    """
    Filter VCF data by minor allele frequency.
    
    Args:
        vcf_df (pd.DataFrame): DataFrame from read_vcf function
        sample_cols (list): List of sample column names
        min_maf (float): Minimum minor allele frequency threshold
        
    Returns:
        pd.DataFrame: Filtered VCF data
    """
    maf_values = calculate_maf(vcf_df, sample_cols)
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