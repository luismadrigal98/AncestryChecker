"""
Utilities for reading and processing VCF files.
This module provides functions to read VCF files, extract relevant information,
and convert it into a format suitable for further analysis.
@module: vcf_reader_utilities
@description: This module provides utilities for reading and processing VCF files.
@date: 2025-04-01
@version: 1.0
@license: MIT

"""

import io
import pandas as pd

## Following function was adopted from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744.

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})