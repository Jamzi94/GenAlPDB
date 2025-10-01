#!/usr/bin/env python3
"""
AlphaFold 3 Runtime Parser

This script parses SLURM output files to extract runtime information for AlphaFold 3 jobs,
including pipeline time, inference time, and total execution time.

Key improvements:
- Adds robust error handling and validation
- Implements better pattern matching with more reliable regex
- Adds detailed logging for better debugging
- Creates a comprehensive summary with additional metrics
"""

import pathlib
import re
import pandas as pd
import logging
import sys
import os
from typing import Dict, List, Optional, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Improved regex patterns for better matching
PATTERN_PIPELINE = re.compile(r'Running data pipeline.*?took ([\d\.]+) seconds')
PATTERN_INFERENCE = re.compile(r'Running model inference.*?took ([\d\.]+) seconds')
PATTERN_TOTAL = re.compile(r'Running fold job.*?Fold job.*?done', re.DOTALL)
PATTERN_GPU_COUNT = re.compile(r'Running on (\d+) GPU')
PATTERN_MEMORY = re.compile(r'Memory usage: ([\d\.]+) GB')
PATTERN_TOKEN_COUNT = re.compile(r'Input has (\d+) tokens')

def parse_runtime_from_logs(logs_dir: Optional[str] = None) -> pd.DataFrame:
    """
    Parses SLURM output files to extract runtime information.
    
    Args:
        logs_dir: Directory containing SLURM output files
        
    Returns:
        DataFrame with parsed runtime information
    """
    # Get the root directory (parent of the script's directory)
    root = pathlib.Path(__file__).resolve().parent.parent
    
    # Set default logs directory if not provided
    if logs_dir is None:
        logs_dir = root / 'logs'
    else:
        logs_dir = pathlib.Path(logs_dir)
    
    # Validate logs directory
    if not logs_dir.exists():
        logger.error(f"Logs directory not found: {logs_dir}")
        raise FileNotFoundError(f"Logs directory not found: {logs_dir}")
    
    logger.info(f"Parsing runtime information from logs in: {logs_dir}")
    
    rows = []
    
    # Find all .out files in the logs directory
    out_files = list(logs_dir.glob('*.out'))
    if not out_files:
        logger.warning(f"No output files found in {logs_dir}")
        return pd.DataFrame()
    
    logger.info(f"Found {len(out_files)} output files to parse")
    
    for f in out_files:
        try:
            content = f.read_text()
            
            # Extract job name from filename
            jobname = f.stem
            
            # Extract runtime information using regex
            pipeline_match = PATTERN_PIPELINE.search(content)
            pipeline_s = float(pipeline_match.group(1)) if pipeline_match else None
            
            inference_match = PATTERN_INFERENCE.search(content)
            inference_s = float(inference_match.group(1)) if inference_match else None
            
            # Extract GPU count if available
            gpu_match = PATTERN_GPU_COUNT.search(content)
            gpu_count = int(gpu_match.group(1)) if gpu_match else None
            
            # Extract memory usage if available
            memory_match = PATTERN_MEMORY.search(content)
            memory_gb = float(memory_match.group(1)) if memory_match else None
            
            # Extract token count if available
            token_match = PATTERN_TOKEN_COUNT.search(content)
            token_count = int(token_match.group(1)) if token_match else None
            
            # Calculate total time (sum of pipeline and inference)
            # If either is missing, use the available one
            if pipeline_s is not None and inference_s is not None:
                total_s = pipeline_s + inference_s
            elif pipeline_s is not None:
                total_s = pipeline_s
            elif inference_s is not None:
                total_s = inference_s
            else:
                total_s = None
            
            # Only add rows with valid timing information
            if total_s is not None:
                rows.append({
                    'jobname': jobname,
                    'pipeline_s': pipeline_s,
                    'inference_s': inference_s,
                    'total_s': total_s,
                    'gpu_count': gpu_count,
                    'memory_gb': memory_gb,
                    'token_count': token_count
                })
                logger.debug(f"Parsed job {jobname}: pipeline={pipeline_s}s, inference={inference_s}s, total={total_s}s")
            else:
                logger.warning(f"No timing information found for job {jobname}")
                
        except Exception as e:
            logger.error(f"Error parsing file {f}: {e}")
    
    # Create DataFrame from parsed data
    if rows:
        df = pd.DataFrame(rows)
        logger.info(f"Successfully parsed {len(df)} jobs with timing information")
        
        # Add derived metrics if possible
        if 'gpu_count' in df.columns and 'total_s' in df.columns:
            df['gpu_seconds'] = df['gpu_count'] * df['total_s']
            
        if 'token_count' in df.columns and 'total_s' in df.columns:
            df['tokens_per_second'] = df['token_count'] / df['total_s']
        
        return df
    else:
        logger.warning("No valid timing information found in any log file")
        return pd.DataFrame()

def save_runtime_summary(df: pd.DataFrame, output_path: Optional[str] = None) -> str:
    """
    Saves runtime summary to CSV and Excel files.
    
    Args:
        df: DataFrame with runtime information
        output_path: Path where to save the summary (without extension)
        
    Returns:
        Path to the saved CSV file
    """
    if df.empty:
        logger.warning("No data to save")
        return ""
    
    # Get the root directory (parent of the script's directory)
    root = pathlib.Path(__file__).resolve().parent.parent
    
    # Set default output path if not provided
    if output_path is None:
        logs_dir = root / 'logs'
        output_path = logs_dir / 'runtime_summary'
    else:
        output_path = pathlib.Path(output_path)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(exist_ok=True)
    
    # Save as CSV
    csv_path = f"{output_path}.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Runtime summary saved to CSV: {csv_path}")
    
    # Save as Excel if possible
    try:
        excel_path = f"{output_path}.xlsx"
        df.to_excel(excel_path, index=False)
        logger.info(f"Runtime summary saved to Excel: {excel_path}")
    except Exception as e:
        logger.warning(f"Could not save to Excel format: {e}")
    
    return csv_path

def main():
    """Main function to run the runtime parsing workflow"""
    try:
        # Parse runtime information from logs
        df = parse_runtime_from_logs()
        
        if not df.empty:
            # Save runtime summary
            csv_path = save_runtime_summary(df)
            
            # Print summary statistics
            logger.info("\nRuntime Summary Statistics:")
            logger.info(f"Total jobs analyzed: {len(df)}")
            
            if 'pipeline_s' in df.columns:
                logger.info(f"Average pipeline time: {df['pipeline_s'].mean():.2f} seconds")
                
            if 'inference_s' in df.columns:
                logger.info(f"Average inference time: {df['inference_s'].mean():.2f} seconds")
                
            if 'total_s' in df.columns:
                logger.info(f"Average total time: {df['total_s'].mean():.2f} seconds")
                
            if 'gpu_count' in df.columns:
                logger.info(f"Average GPU count: {df['gpu_count'].mean():.2f}")
                
            if 'token_count' in df.columns:
                logger.info(f"Average token count: {df['token_count'].mean():.2f}")
                
            logger.info(f"Full details saved to: {csv_path}")
        else:
            logger.warning("No runtime data found to analyze")
            
        logger.info("Runtime parsing completed successfully")
        
    except Exception as e:
        logger.error(f"Runtime parsing failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

