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
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from src.vcf_reader_utilities import read_vcf
from src.analysis_utilities import determine_ancestry, plot_ancestry
from src.data_tidyer import (
    read_relationship_map, filter_vcf_data, identify_informative_snps,
    filter_biallelic_snps, filter_by_maf, filter_by_missing_rate, filter_by_qual,
    filter_by_region, filter_founder_homozygous
)

# Setting up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='AncestryChecker: Check the ancestry of F2 individuals.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-r', '--relationships', required=True, help='Path to the relationship map file')
    parser.add_argument('-o', '--output', default='ancestry_results', help='Output directory')
    parser.add_argument('--no-missing', action='store_true', help='Do not allow missing data in samples')
    parser.add_argument('--ref-founder', type=str, help='Specify the reference founder (default: first founder)')

    # QC-related arguments
    parser.add_argument('--biallelic-only', action='store_true', 
                        help='Filter to keep only biallelic SNPs')
    parser.add_argument('--min-maf', type=float, default=0.0, 
                        help='Minimum minor allele frequency (default: 0.0 = no filtering)')
    parser.add_argument('--max-missing-rate', type=float, default=1.0, 
                        help='Maximum rate of missing data per SNP (default: 1.0 = no filtering)')
    parser.add_argument('--min-qual', type=float, default=0.0, 
                        help='Minimum QUAL value for variants (default: 0.0 = no filtering)')
    parser.add_argument('--founders_homozygous', action='store_true',
                        help='Check if founders are homozygous for all variants')
    parser.add_argument('--retain-informative-only', action='store_true',
                        help='Retain only informative SNPs for ancestry analysis')
    parser.add_argument('--allele_depth_thr',
                        type=float, default=0.0,
                        help='Minimum allele depth threshold for filtering (default: 0.0 = no filtering). THIS IS ONLY VALID FOR BCFTOOLS VCFs.',
                        required=False)

    # Region targeting arguments
    parser.add_argument(
        '--ROI',
        default=None,
        help='Region of interest in the format "chrom:start-end" (e.g., "1:1000-2000") or "chrom" only for a specific chromosome'
        )
    
    # Additional arguments
    parser.add_argument(
        '--vcf_from_caller',
        type = str,
        default = 'freebayes',
        choices = ['freebayes', 'gatk', 'bcftools'],
        help = 'Variant caller used to generate the VCF file (default: freebayes)'
    )


    return parser.parse_args()

def main():
    """Main function for the AncestryChecker program."""
    args = parse_arguments()
    
    # Read input files
    logger.info(f"Reading VCF file: {args.vcf}")
    vcf_data = read_vcf(args.vcf)
    format_fields = vcf_data['FORMAT'].iloc[0].split(':')
    total_variants = len(vcf_data)  # Total before any filtering
    
    # Apply region filtering if specified
    if args.ROI is not None:
        logger.info("Applying region filtering...")
        only_chrom = False if ':' in args.ROI else True
        if only_chrom:
            args.chrom = args.ROI
            args.start_pos = None
            args.end_pos = None
        else:
            chrom, positions = args.ROI.split(':')
            start_pos, end_pos = map(int, positions.split('-'))
            args.chrom = chrom
            args.start_pos = start_pos
            args.end_pos = end_pos
        region_str = f"Chromosome {args.chrom}:{args.start_pos}-{args.end_pos}" if args.start_pos and args.end_pos else ("Chromosome " + args.chrom)
        logger.info(f"Targeting: {region_str}")
        
        vcf_data = filter_by_region(vcf_data, args.chrom, args.start_pos, args.end_pos)
        
        # After region filtering, set initial count for subsequent reporting
        current_len = len(vcf_data)
        
        logger.info(f"Retained {current_len} variants in target region (from total {total_variants})")
        
        if len(vcf_data) == 0:
            logger.error("Error: No variants found in specified region. Check chromosome name and positions.")
            return
    
    logger.info(f"Reading relationship map: {args.relationships}")
    relationships = read_relationship_map(args.relationships)
    
    # Extract sample names
    format_col_idx = vcf_data.columns.get_loc('FORMAT')
    sample_cols = vcf_data.columns[format_col_idx + 1:]
    
    # Get all unique founder names from the relationship map
    founders = set(relationships['Founder1'].tolist() + relationships['Founder2'].tolist())
    logger.info(f"Found {len(founders)} unique founders: {', '.join(founders)}")
    
    # Get all F2 samples from the relationship map
    f2_samples = relationships[relationships['Generation'] == 'F2']['Sample'].tolist()
    logger.info(f"Found {len(f2_samples)} F2 samples: {', '.join(f2_samples)}")
    
    logger.info(f"Initial variant count: {current_len}")
    
    # Apply QC filters
    if args.biallelic_only:
        logger.info("Filtering to keep only biallelic SNPs...")
        vcf_data = filter_biallelic_snps(vcf_data)
        current_len = len(vcf_data)
        logger.info(f"Retained {current_len} biallelic SNPs")
    
    if args.min_qual > 0:
        print(f"Filtering SNPs with QUAL < {args.min_qual}...")
        vcf_data = filter_by_qual(vcf_data, args.min_qual)
        print(f"Retained {len(vcf_data)} SNPs after QUAL filtering")

    if args.allele_depth_thr > 0:
        print(f"Filtering SNPs with allele depth < {args.allele_depth_thr}...")
        vcf_data = filter_by_qual(vcf_data, args.allele_depth_thr)
        print(f"Retained {len(vcf_data)} SNPs after allele depth filtering")

    # Filter according to MAF values
    if args.min_maf > 0:
        logger.info(f"Filtering SNPs with MAF < {args.min_maf}...")
        vcf_data = filter_by_maf(vcf_data, sample_cols, args.min_maf, format_fields, args.vcf_from_caller)
        current_len = len(vcf_data)
        logger.info(f"Retained {current_len} SNPs after MAF filtering")

    # Filter VCF data for ancestry analysis
    logger.info("Filtering VCF data for founder/sample completeness...")
    filtered_vcf = filter_vcf_data(vcf_data, list(founders), 
                                    f2_samples, not args.no_missing)
    current_len = len(filtered_vcf)
    logger.info(f"Retained {current_len} SNPs after basic data filtering")
    
    if args.founders_homozygous:
        logger.info("Checking if founders are homozygous for all variants...")
        filtered_vcf = filter_founder_homozygous(filtered_vcf, list(founders))
        current_len = len(filtered_vcf)
        logger.info(f"Retained {current_len} variants after checking founders' homozygosity")

    # Extract genotype columns
    all_samples = list(founders) + f2_samples
    
    if args.max_missing_rate < 1.0:
        logger.info(f"Filtering SNPs with missing rate > {args.max_missing_rate}...")
        filtered_vcf = filter_by_missing_rate(filtered_vcf, all_samples, args.max_missing_rate)
        current_len = len(filtered_vcf)
        logger.info(f"Retained {len(filtered_vcf)} SNPs after missing rate filtering")
    
    # Identify informative SNPs
    if args.retain_informative_only:
        logger.info("Identifying informative SNPs...")
        founder_cols = [col for col in sample_cols if col in founders]
        filtered_vcf = identify_informative_snps(filtered_vcf, founder_cols)
        logger.info(f"Final count: {len(filtered_vcf)} informative SNPs (from initial {total_variants})")
    
    # Process each F2 sample
    for f2_id in f2_samples:
        f2_col = f"{f2_id}_GT"
        
        if f2_col in filtered_vcf.columns:
            print(f"Analyzing ancestry for sample {f2_id}...")
            ancestry_results = determine_ancestry(filtered_vcf, relationships, f2_col, args.ref_founder)

            if 'Novel' in ancestry_results['Ancestry'].values:
                novel_count = (ancestry_results['Ancestry'] == 'Novel').sum()
                print(f"Found {novel_count} novel variants ({novel_count/len(ancestry_results)*100:.2f}%)")
                
                # Show a few examples of novel variants for debugging
                print("Examples of novel variants:")
                novel_examples = ancestry_results[ancestry_results['Ancestry'] == 'Novel'].head(3)
                for _, row in novel_examples.iterrows():
                    print(f"  CHROM={row['CHROM']}, POS={row['POS']}")
                    # Use dynamic column names instead of hardcoded ones
                    f2_col = 'F2'  # This should always be present
                    # Get founder columns (excluding standard columns)
                    founder_cols = [col for col in row.index if col not in 
                                ['CHROM', 'POS', 'REF', 'ALT', 'F2', 'Ancestry']]
                    if founder_cols:
                        print(f"    F2={row[f2_col]}, {founder_cols[0]}={row[founder_cols[0]]}, " +
                            f"{founder_cols[1]}={row[founder_cols[1]] if len(founder_cols) > 1 else 'N/A'}")
            
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