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
from src.analysis_utilities import determine_ancestry, plot_ancestry
from src.data_tidyer import (
    read_relationship_map, filter_vcf_data, identify_informative_snps,
    filter_biallelic_snps, filter_by_maf, filter_by_missing_rate, filter_by_qual,
    filter_by_region
)
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='AncestryChecker: Check the ancestry of F2 individuals.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-r', '--relationships', required=True, help='Path to the relationship map file')
    parser.add_argument('-o', '--output', default='ancestry_results', help='Output directory')
    parser.add_argument('--no-missing', action='store_true', help='Do not allow missing data in samples')
    
    # QC-related arguments
    parser.add_argument('--biallelic-only', action='store_true', 
                        help='Filter to keep only biallelic SNPs')
    parser.add_argument('--min-maf', type=float, default=0.0, 
                        help='Minimum minor allele frequency (default: 0.0 = no filtering)')
    parser.add_argument('--max-missing-rate', type=float, default=1.0, 
                        help='Maximum rate of missing data per SNP (default: 1.0 = no filtering)')
    parser.add_argument('--min-qual', type=float, default=0.0, 
                        help='Minimum QUAL value for variants (default: 0.0 = no filtering)')
    
    # Region targeting arguments
    parser.add_argument('-c', '--chrom', type=str, 
                        help='Target a specific chromosome (e.g., "1", "Chr1")')
    parser.add_argument('--start-pos', type=int, 
                        help='Start position for targeted analysis')
    parser.add_argument('--end-pos', type=int, 
                        help='End position for targeted analysis')
    
    return parser.parse_args()

def main():
    """Main function for the AncestryChecker program."""
    args = parse_arguments()
    
    # Read input files
    print(f"Reading VCF file: {args.vcf}")
    vcf_data = read_vcf(args.vcf)
    initial_count = len(vcf_data)
    
    # Apply region filtering if specified
    if args.chrom or args.start_pos or args.end_pos:
        print("Applying region filtering...")
        region_str = f"Chromosome {args.chrom if args.chrom else 'All'}"
        if args.start_pos:
            region_str += f", from position {args.start_pos}"
        if args.end_pos:
            region_str += f" to {args.end_pos}"
        print(f"Targeting: {region_str}")
        
        vcf_data = filter_by_region(vcf_data, args.chrom, args.start_pos, args.end_pos)
        print(f"Retained {len(vcf_data)} variants in target region")
        
        if len(vcf_data) == 0:
            print("Error: No variants found in specified region. Check chromosome name and positions.")
            return
    
    print(f"Reading relationship map: {args.relationships}")
    relationships = read_relationship_map(args.relationships)
    
    # Ensure all columns are strings
    relationships = relationships.astype(str)
    
    # Extract sample names
    format_col_idx = vcf_data.columns.get_loc('FORMAT')
    sample_cols = vcf_data.columns[format_col_idx + 1:]
    
    # Get all unique founder names from the relationship map
    founders = set(relationships['Founder1'].tolist() + relationships['Founder2'].tolist())
    print(f"Found {len(founders)} unique founders: {', '.join(founders)}")
    
    # Get all F2 samples from the relationship map
    f2_samples = relationships[relationships['Generation'] == 'F2']['Sample'].tolist()
    print(f"Found {len(f2_samples)} F2 samples: {', '.join(f2_samples)}")
    
    print(f"Initial SNP count: {initial_count}")
    
    # Apply QC filters
    if args.biallelic_only:
        print("Filtering to keep only biallelic SNPs...")
        vcf_data = filter_biallelic_snps(vcf_data)
        print(f"Retained {len(vcf_data)} biallelic SNPs")
    
    # Filter VCF data for ancestry analysis
    print("Filtering VCF data for founder/sample completeness...")
    filtered_vcf = filter_vcf_data(vcf_data, list(founders), f2_samples, not args.no_missing)
    print(f"Retained {len(filtered_vcf)} SNPs after missing data filtering")
    
    # Extract genotype columns
    all_samples = list(founders) + f2_samples
    
    if args.min_maf > 0:
        print(f"Filtering SNPs with MAF < {args.min_maf}...")
        filtered_vcf = filter_by_maf(filtered_vcf, all_samples, args.min_maf)
        print(f"Retained {len(filtered_vcf)} SNPs after MAF filtering")
    
    if args.max_missing_rate < 1.0:
        print(f"Filtering SNPs with missing rate > {args.max_missing_rate}...")
        filtered_vcf = filter_by_missing_rate(filtered_vcf, all_samples, args.max_missing_rate)
        print(f"Retained {len(filtered_vcf)} SNPs after missing rate filtering")
    
    if args.min_qual > 0:
        print(f"Filtering SNPs with QUAL < {args.min_qual}...")
        filtered_vcf = filter_by_qual(filtered_vcf, args.min_qual)
        print(f"Retained {len(filtered_vcf)} SNPs after QUAL filtering")
    
    # Identify informative SNPs
    print("Identifying informative SNPs...")
    founder_cols = [col for col in sample_cols if col in founders]
    informative_vcf = identify_informative_snps(filtered_vcf, founder_cols)
    
    print(f"Final count: {len(informative_vcf)} informative SNPs (from initial {initial_count})")
    
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