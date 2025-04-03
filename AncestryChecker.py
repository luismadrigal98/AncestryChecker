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
from src.analysis_utilities import determine_ancestry, plot_ancestry
from src.data_tidyer import read_relationship_map, filter_vcf_data, identify_informative_snps

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='AncestryChecker: Check the ancestry of F2 individuals.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-r', '--relationships', required=True, help='Path to the relationship map file')
    parser.add_argument('-o', '--output', default='ancestry_results', help='Output directory')
    parser.add_argument('--no-missing', action='store_true', help='Do not allow missing data in samples')
    
    return parser.parse_args()

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