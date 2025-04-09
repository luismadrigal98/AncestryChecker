# Ancestry Checker

Ancestry Checker is a set of tools for analyzing the genetic ancestry of F2 individuals in a nested association mapping (NAM) design framework. The package provides two implementations:

1. **AncestryChecker**: A full-featured tool with extensive filtering options and customization.
2. **AncestryCheckerLite**: A simplified version that assumes the reference genome corresponds to the 664c founder line.

## Overview

This program analyzes the ancestry of F2 individuals by comparing their genotypes with those of the founder lines. The analysis assumes that the F2 genome is a mosaic of the founder genomes. The program performs quality checks on the VCF data, filtering variants to retain only those suitable for ancestry analysis.

## Requirements

- Python 3.6+
- Dependencies:
  - pandas
  - numpy
  - matplotlib
  - logging

## Installation

Clone the repository:

```bash
git clone https://github.com/username/AncestryChecker.git
cd AncestryChecker
```

No additional installation steps are required, as the scripts can be run directly.

## Input Files

### 1. VCF File

A standard VCF file containing genotype information for both founders and F2 samples.

### 2. Relationship Map

A tab-separated text file defining the relationships between F2 individuals and founder lines. The format is:

```
[Sample ID]  F2  [Founder1]  [Founder2]
```

Example:
```
17  F2  444c    664c
22  F2  767c    664c
604 F2  1192c   664c
675 F2  155c    664c
```

Each row represents one F2 sample and its two founder parents. The file has no header row.

## Usage

### AncestryChecker (Full Version)

```bash
python AncestryChecker.py -v [VCF_FILE] -r [RELATIONSHIP_MAP] -o [OUTPUT_DIR] [OPTIONS]
```

#### Options:

- `-v, --vcf`: Path to VCF file (required)
- `-r, --relationships`: Path to relationship map file (required)
- `-o, --output`: Output directory (default: 'ancestry_results')
- `--no-missing`: Do not allow missing data in samples
- `--ref-founder`: Specify the reference founder (default: first founder)

##### Quality Control Options:
- `--biallelic-only`: Filter to keep only biallelic SNPs
- `--min-maf`: Minimum minor allele frequency (default: 0.0 = no filtering)
- `--max-missing-rate`: Maximum rate of missing data per SNP (default: 1.0 = no filtering)
- `--min-qual`: Minimum QUAL value for variants (default: 0.0 = no filtering)
- `--founders_homozygous`: Check if founders are homozygous for all variants
- `--retain-informative-only`: Retain only informative SNPs for ancestry analysis
- `--allele_depth_thr`: Minimum allele depth threshold (requires bcftools VCF)

##### Region Selection:
- `--ROI`: Region of interest in format "chrom:start-end" or just "chrom" 

##### VCF Source:
- `--vcf_from_caller`: Variant caller used: 'freebayes' (default), 'gatk', or 'bcftools'

### AncestryCheckerLite (Simplified Version)

```bash
python AncestryCheckerLite.py -v [VCF_FILE] -r [RELATIONSHIP_MAP] -o [OUTPUT_DIR] [OPTIONS]
```

#### Options:
- `-v, --vcf`: Path to VCF file (required)
- `-r, --relationships`: Path to relationship map file (required)
- `-o, --output`: Output directory (default: 'ancestry_results')
- `--region`: Region of interest (e.g., "Chr_01" or "Chr_01:1000000-2000000")

## Output Files

For each F2 sample, the program produces:

1. **CSV File**: Detailed ancestry results containing variant positions and ancestry assignments
   - Format: `[output_dir]/[sample_id]_ancestry_details.csv`

2. **Plot Files**: 
   - For AncestryChecker: PDF files with chromosome-specific ancestry proportion plots
   - For AncestryCheckerLite: PNG files showing overall ancestry proportions and by-chromosome SNP counts

## Analysis Approach

The program assigns ancestry to each variant in an F2 individual by comparing its genotype with those of the founder lines:

- **AncestryChecker**: Uses a sophisticated approach that analyzes allele combinations to determine if a variant matches one founder, both founders, or represents a novel combination
  
- **AncestryCheckerLite**: Uses a simplified approach that assumes the reference allele corresponds to the 664c founder line

The ancestry categories include:
- Founder-specific ancestry (e.g., '664c', 'Alternative')
- Mixed ancestry between founders
- Novel variants not matching either founder
- Unresolved variants (when founders have identical genotypes)
- Missing data

## Examples

### Basic Usage (Full Version)
```bash
python AncestryChecker.py -v data/variants.vcf -r examples/relationship_map.txt -o results
```

### Basic Usage (Lite Version)
```bash
python AncestryCheckerLite.py -v data/variants.vcf -r examples/relationship_map.txt -o results
```

### Advanced Example with Filtering
```bash
python AncestryChecker.py -v data/variants.vcf -r examples/relationship_map.txt -o filtered_results \
  --biallelic-only --min-maf 0.05 --max-missing-rate 0.2 --min-qual 30 \
  --retain-informative-only --ROI Chr_01:1000000-2000000
```

## Author

Luis Javier Madrigal-Roca
