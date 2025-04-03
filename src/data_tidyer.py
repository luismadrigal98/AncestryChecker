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
    return pd.read_csv(path, sep='\t', header=None, 
                        names=['Sample', 'Generation', 'Founder1', 'Founder2'])

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
    
    # Ensure all founders and at least one F2 sample are in the data
    if not all(founder in founder_cols for founder in founders):
        raise ValueError("Not all founders found in VCF file")
    if not any(f2 in f2_cols for f2 in f2_samples):
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