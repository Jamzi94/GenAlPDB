#!/usr/bin/env python3
"""
Remote Job Submit Script for AlphaFold 3 Workflow (SLURM Version)

This script handles the submission of AlphaFold 3 jobs to SLURM scheduler.
It processes JSON input files and submits them to the SLURM computing environment.

Key features:
- Validates input JSON files before submission
- Generates and submits SLURM job scripts
- Provides detailed logging of submission process
- Supports MSA detection and conditional job submission
- Manages job dependencies between MSA and inference jobs
"""

import os
import sys
import glob
import pathlib
import argparse
import logging
import subprocess
from typing import List, Dict, Optional, Tuple, Set

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def find_json_files(json_dir: pathlib.Path) -> List[pathlib.Path]:
    """
    Find all valid JSON input files in the specified directory.
    
    Args:
        json_dir: Path to the directory containing JSON files
        
    Returns:
        List of paths to valid JSON files
    """
    if not json_dir.exists():
        logger.error(f"JSON directory does not exist: {json_dir}")
        return []
    
    # Collect valid JSON files (excluding chain mapping files)
    valid_json_files = []
    for json_file in json_dir.glob("*.json"):
        if "_chain_mapping.json" not in json_file.name:
            valid_json_files.append(json_file)
    
    return valid_json_files

def extract_base_id(json_file: pathlib.Path) -> str:
    """
    Extract base ID from JSON filename (first two parts separated by _).
    
    Args:
        json_file: Path to the JSON file
        
    Returns:
        Base ID string
    """
    base_name = json_file.stem
    parts = base_name.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return base_name

def group_json_files_by_base_id(json_files: List[pathlib.Path]) -> Dict[str, List[pathlib.Path]]:
    """
    Group JSON files by their base ID.
    
    Args:
        json_files: List of JSON file paths
        
    Returns:
        Dictionary mapping base IDs to lists of JSON files
    """
    base_id_files = {}
    
    for json_file in json_files:
        base_id = extract_base_id(json_file)
        if base_id not in base_id_files:
            base_id_files[base_id] = []
        base_id_files[base_id].append(json_file)
    
    return base_id_files

def submit_slurm_job(root_dir: pathlib.Path) -> bool:
    """
    Submit AlphaFold 3 jobs to SLURM using the advanced submission script.
    
    Args:
        root_dir: Path to the project root directory
        
    Returns:
        True if submission was successful, False otherwise
    """
    try:
        # Check if slurm_job_submit.sh exists at the root
        slurm_script = root_dir / "slurm_job_submit.sh"
        if not slurm_script.exists():
            logger.error(f"SLURM job script not found: {slurm_script}")
            logger.error("Please ensure the script exists and is executable.")
            return False
        
        # Make sure the script is executable
        try:
            os.chmod(slurm_script, 0o755)
        except Exception as e:
            logger.warning(f"Could not set executable permission on {slurm_script}: {e}")
        
        # Execute the SLURM submission script
        logger.info(f"Executing SLURM job submission script: {slurm_script}")
        
        result = subprocess.run(
            ["bash", str(slurm_script), str(root_dir)],
            check=True,
            cwd=root_dir,
            capture_output=True,
            text=True
        )
        
        # Log the output
        if result.stdout:
            for line in result.stdout.strip().split('\n'):
                logger.info(f"SLURM: {line}")
        
        if result.stderr:
            for line in result.stderr.strip().split('\n'):
                logger.warning(f"SLURM stderr: {line}")
        
        logger.info("SLURM job submission completed successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"SLURM job submission failed with exit code {e.returncode}")
        if e.stdout:
            logger.error("Submission output:")
            for line in e.stdout.strip().split('\n'):
                logger.error(f"  {line}")
        if e.stderr:
            logger.error("Submission errors:")
            for line in e.stderr.strip().split('\n'):
                logger.error(f"  {line}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during SLURM job submission: {e}")
        return False

def main():
    """Main function to process and submit AlphaFold 3 jobs to SLURM."""
    parser = argparse.ArgumentParser(
        description='Submit AlphaFold 3 jobs to SLURM for JSON input files.'
    )
    parser.add_argument('root_dir', nargs='?', help='Path to the project root directory')
    
    args = parser.parse_args()
    
    # Get root directory
    if args.root_dir:
        root_dir = pathlib.Path(args.root_dir)
    else:
        # Default to current working directory
        root_dir = pathlib.Path.cwd()
    
    logger.info(f"Project root: {root_dir}")
    
    # Find JSON input files
    json_dir = root_dir / "af3_json"
    json_files = find_json_files(json_dir)
    
    if not json_files:
        logger.error("No valid JSON input files found. Cannot proceed.")
        sys.exit(1)
    
    logger.info(f"Found {len(json_files)} JSON input files to process.")
    
    # Group files by base ID
    base_id_files = group_json_files_by_base_id(json_files)
    logger.info(f"Found {len(base_id_files)} unique protein base IDs to process.")
    
    # Submit jobs using the SLURM submission script
    if submit_slurm_job(root_dir):
        logger.info("All jobs submitted successfully to SLURM.")
    else:
        logger.error("Job submission failed. Check the logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()

