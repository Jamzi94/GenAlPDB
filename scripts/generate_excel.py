#!/usr/bin/env python3
"""
AlphaFold 3 Excel Results Generator

This script converts runtime summary CSV data into an Excel format with enhanced formatting.
It parses job names to extract metadata (UniProt ID, GPU, CPUs, memory, partition) and
creates a comprehensive Excel workbook with the benchmark results.

Prerequisites:
- Runtime summary CSV must exist in logs/runtime_summary.csv
- Run parse_slurm_time.py first to generate the summary

Output:
- benchmark_results.xlsx in the project root directory
"""

import pandas as pd
import pathlib

# Get project root directory (parent of scripts/)
root = pathlib.Path(__file__).resolve().parent.parent

# Path to the runtime summary CSV file
summary = root / 'logs' / 'runtime_summary.csv'

# Check if the summary file exists
if not summary.exists():
    print('Error: Runtime summary not found. Please run parse_slurm_time.py first.')
    exit()

# Load the runtime summary data
df = pd.read_csv(summary)

def split(jobname):
    """
    Extract metadata from job name.
    
    Job names are expected to follow the format:
    uniprot_gpu_cpus_mem_partition
    
    Args:
        jobname: Job name string to parse
        
    Returns:
        pandas.Series with extracted metadata fields
    """
    parts = jobname.split('_')
    return pd.Series(dict(
        uniprot=parts[0],
        gpu=parts[1],
        cpus=parts[2],
        mem=parts[3],
        partition=parts[4]
    ))

# Parse job names and add metadata columns to the dataframe
df = df.join(df['jobname'].apply(split))

# Export to Excel format
df.to_excel(root / 'benchmark_results.xlsx', index=False)

print('Excel file successfully generated: benchmark_results.xlsx')
