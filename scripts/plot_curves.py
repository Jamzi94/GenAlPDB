
#!/usr/bin/env python3
"""
AlphaFold 3 Performance Visualization Tool

This script creates visualizations of AlphaFold 3 benchmark results,
plotting the relationship between GPU count and execution time.

Key improvements:
- Adds proper error handling and validation
- Implements better visualization with labels, titles, and grid
- Saves multiple plot formats for better compatibility
- Adds detailed logging for better debugging
"""

import pandas as pd
import matplotlib.pyplot as plt
import pathlib
import logging
import sys
import os
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def create_performance_plots(excel_path: Optional[str] = None, output_dir: Optional[str] = None):
    """
    Creates performance visualization plots from benchmark results.
    
    This function generates publication-quality plots showing the relationship
    between computational resources (GPU count) and performance metrics
    (execution time, memory usage). The plots help analyze scalability
    and resource efficiency of AlphaFold 3 predictions.
    
    Generated plots:
    - time_vs_gpu.png/pdf/svg: Execution time vs GPU count
    - memory_vs_gpu.png: Memory usage vs GPU count (if available)
    
    Args:
        excel_path: Path to the Excel file with benchmark results
        output_dir: Directory where plots will be saved
    """
    # Get the root directory (parent of the script's directory)
    root = pathlib.Path(__file__).resolve().parent.parent
    
    # Set default paths if not provided
    if excel_path is None:
        excel_path = root / 'benchmark_results.xlsx'
    else:
        excel_path = pathlib.Path(excel_path)
        
    if output_dir is None:
        output_dir = root / 'output'
    else:
        output_dir = pathlib.Path(output_dir)
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Validate input file
    if not excel_path.exists():
        logger.error(f"Excel file not found: {excel_path}")
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    
    # Load data
    try:
        logger.info(f"Loading benchmark data from: {excel_path}")
        df = pd.read_excel(excel_path)
        
        if df.empty:
            logger.warning("Excel file contains no data")
            return
            
        # Check required columns
        required_columns = ['gpu', 'elapsed_s']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Missing required columns: {missing_columns}")
            raise ValueError(f"Missing required columns: {missing_columns}")
            
        logger.info(f"Successfully loaded data with {len(df)} benchmark results")
    except Exception as e:
        logger.error(f"Failed to load Excel file: {e}")
        raise
    
    # Create time vs GPU plot
    try:
        logger.info("Creating Time vs GPU plot")
        plt.figure(figsize=(10, 6))
        
        # Create scatter plot with larger points and transparency for overlapping points
        plt.scatter(df['gpu'], df['elapsed_s'], s=100, alpha=0.7)
        
        # Add title and labels with appropriate font sizes
        plt.title('AlphaFold 3 Performance: Execution Time vs GPU Count', fontsize=14)
        plt.xlabel('Number of GPUs', fontsize=12)
        plt.ylabel('Execution Time (seconds)', fontsize=12)
        
        # Add grid for better readability
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Improve aesthetics with tight layout
        plt.tight_layout()
        
        # Save in multiple formats for different use cases
        plot_path_png = output_dir / 'time_vs_gpu.png'
        plot_path_pdf = output_dir / 'time_vs_gpu.pdf'
        plot_path_svg = output_dir / 'time_vs_gpu.svg'
        
        plt.savefig(plot_path_png, dpi=300)  # High-resolution PNG
        plt.savefig(plot_path_pdf)  # Vector format for publications
        plt.savefig(plot_path_svg)  # Vector format for web
        
        logger.info(f"Plots saved to: {plot_path_png}, {plot_path_pdf}, {plot_path_svg}")
        
        # Create additional plots if more columns are available
        if 'memory_gb' in df.columns:
            logger.info("Creating Memory vs GPU plot")
            plt.figure(figsize=(10, 6))
            plt.scatter(df['gpu'], df['memory_gb'], s=100, alpha=0.7)
            plt.title('AlphaFold 3 Performance: Memory Usage vs GPU Count', fontsize=14)
            plt.xlabel('Number of GPUs', fontsize=12)
            plt.ylabel('Memory Usage (GB)', fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            memory_plot_path = output_dir / 'memory_vs_gpu.png'
            plt.savefig(memory_plot_path, dpi=300)
            logger.info(f"Memory plot saved to: {memory_plot_path}")
    
    except Exception as e:
        logger.error(f"Failed to create plots: {e}")
        raise

def main():
    """Main function to run the visualization workflow"""
    try:
        create_performance_plots()
        logger.info("Performance visualization completed successfully")
    except Exception as e:
        logger.error(f"Performance visualization failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

