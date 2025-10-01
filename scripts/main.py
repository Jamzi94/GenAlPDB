#!/usr/bin/env python3
"""
Main script for AlphaFold 3 benchmark workflow with SLURM integration and MSA fusion

This script manages the complete workflow:
1. Checks if JSON input files exist and are up-to-date
2. Generates JSON files from template if needed (with SMILES support)
3. Launches the SLURM job submission process via remote_job_submit.py

Key features:
- Automatic JSON generation with SMILES support
- Change detection and intelligent regeneration
- SLURM backend integration for HPC environments
- Comprehensive logging and error handling
"""

import subprocess
import pathlib
import sys
import os
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def get_file_modification_time(file_path: pathlib.Path) -> float:
    """
    Get the modification time of a file.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Modification time as a float (timestamp)
    """
    try:
        return file_path.stat().st_mtime
    except (OSError, FileNotFoundError):
        return 0.0

def check_json_files_exist(json_dir: pathlib.Path) -> tuple[bool, list[pathlib.Path]]:
    """
    Check if valid JSON input files exist in the JSON directory.
    
    Args:
        json_dir: Path to the JSON directory
        
    Returns:
        Tuple of (files_exist, list_of_json_files)
    """
    if not json_dir.exists():
        return False, []
    
    # Collect valid JSON files (excluding chain mapping files)
    valid_json_files = []
    for json_file in json_dir.glob("*.json"):
        if "_chain_mapping.json" not in json_file.name:
            valid_json_files.append(json_file)
    
    return len(valid_json_files) > 0, valid_json_files

def should_regenerate_json(root: pathlib.Path, json_files: list[pathlib.Path]) -> bool:
    """
    Determine if JSON files should be regenerated based on template modification time.
    
    This function implements intelligent change detection by comparing modification
    timestamps. JSON files are regenerated if:
    1. The template CSV has been modified more recently than any JSON file
    2. The SMILES database has been modified more recently than any JSON file
    
    Args:
        root: Root directory path
        json_files: List of existing JSON files
        
    Returns:
        True if JSON files should be regenerated, False otherwise
    """
    template_path = root / "config" / "af3_template.csv"
    
    if not template_path.exists():
        return True  # No template exists, cannot determine - regenerate
    
    if not json_files:
        return True  # No JSON files exist, need to generate
    
    # Get template modification time
    template_mtime = get_file_modification_time(template_path)
    
    # Check if any JSON file is older than the template
    for json_file in json_files:
        json_mtime = get_file_modification_time(json_file)
        if json_mtime < template_mtime:
            logger.info(f"Template is newer than {json_file.name}, regeneration needed")
            return True
    
    # Also check if SMILES database is newer (if it exists)
    smiles_db_path = root / "config" / "smiles_database.csv"
    if smiles_db_path.exists():
        smiles_mtime = get_file_modification_time(smiles_db_path)
        for json_file in json_files:
            json_mtime = get_file_modification_time(json_file)
            if json_mtime < smiles_mtime:
                logger.info(f"SMILES database is newer than {json_file.name}, regeneration needed")
                return True
    
    return False

def generate_json_files(root: pathlib.Path, force_regenerate: bool = False) -> bool:
    """
    Generate JSON files using the enhanced generator script with SMILES support.
    
    This function invokes the JSON generator script which:
    - Fetches protein sequences from UniProt
    - Resolves ligands and ions (with SMILES support if available)
    - Creates AlphaFold 3-compatible JSON input files
    
    Args:
        root: Root directory path
        force_regenerate: If True, regenerate even if files exist
        
    Returns:
        True if generation was successful, False otherwise
    """
    if force_regenerate:
        logger.info("Regenerating JSON files from template...")
    else:
        logger.info("Generating JSON files from template...")
    
    # Use the enhanced generator in scripts/
    generator_script = root / "scripts" / "generate_af3_json.py"
    
    if not generator_script.exists():
        logger.error(f"JSON generator script not found: {generator_script}")
        return False
    
    logger.info(f"Using JSON generator: {generator_script.name}")
    
    # Check if SMILES support is available
    smiles_support = (
        (root / "scripts" / "ligand_resolver.py").exists() and 
        (root / "config" / "smiles_database.csv").exists()
    )
    
    if smiles_support:
        logger.info("✓ SMILES support available - compound names can be used in templates")
        # Show number of available compounds
        try:
            with open(root / "config" / "smiles_database.csv", 'r') as f:
                compound_count = len(f.readlines()) - 1  # Subtract header
            logger.info(f"✓ {compound_count} pharmaceutical compounds available")
        except Exception as e:
            logger.warning(f"Could not read SMILES database: {e}")
    else:
        logger.info("Using traditional CCD-only mode")
    
    try:
        # Run the JSON generator
        result = subprocess.run(
            [sys.executable, str(generator_script)],
            cwd=root,
            check=True,
            capture_output=True,
            text=True
        )
        
        logger.info("✓ JSON generation completed successfully")
        
        # Parse and display useful information from the output
        if result.stdout:
            lines = result.stdout.strip().split('\n')
            for line in lines:
                if "JSON created:" in line or "Successfully" in line or "loaded" in line:
                    logger.info(f"  {line}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"✗ JSON generation failed with exit code {e.returncode}")
        if e.stdout:
            logger.error("Generator output:")
            for line in e.stdout.strip().split('\n'):
                logger.error(f"  {line}")
        if e.stderr:
            logger.error("Generator errors:")
            for line in e.stderr.strip().split('\n'):
                logger.error(f"  {line}")
        return False
    except Exception as e:
        logger.error(f"✗ Unexpected error during JSON generation: {e}")
        return False

def check_template_file(root: pathlib.Path) -> bool:
    """
    Check if the template CSV file exists and has data.
    
    Args:
        root: Root directory path
        
    Returns:
        True if template is valid, False otherwise
    """
    template_path = root / "config" / "af3_template.csv"
    
    if not template_path.exists():
        logger.error(f"✗ Template file not found: {template_path}")
        logger.error("Please create the template file with your protein and ligand specifications.")
        return False
    
    # Check if template has data (more than just header)
    try:
        with open(template_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) <= 1:
            logger.error("✗ Template file exists but contains no data rows.")
            logger.error("Please add protein entries to the template file.")
            return False
        
        logger.info(f"✓ Template file found with {len(lines) - 1} data rows")
        
        # Show template modification time
        mtime = datetime.fromtimestamp(get_file_modification_time(template_path))
        logger.info(f"  Last modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
        
        return True
        
    except Exception as e:
        logger.error(f"✗ Error reading template file: {e}")
        return False

def show_workflow_summary(root: pathlib.Path, json_files: list[pathlib.Path]):
    """
    Show a summary of the workflow status.
    
    Args:
        root: Root directory path
        json_files: List of JSON files
    """
    logger.info("=" * 50)
    logger.info("Workflow Summary")
    logger.info("=" * 50)
    
    # Template info
    template_path = root / "config" / "af3_template.csv"
    if template_path.exists():
        try:
            with open(template_path, 'r') as f:
                template_lines = len(f.readlines()) - 1
            logger.info(f"Template entries: {template_lines}")
        except Exception as e:
            logger.warning(f"Could not read template file: {e}")
            logger.info("Template: Present")
    
    # JSON files info
    logger.info(f"JSON input files: {len(json_files)}")
    if json_files:
        logger.info("Sample files:")
        for i, json_file in enumerate(json_files[:3]):
            logger.info(f"  - {json_file.name}")
        if len(json_files) > 3:
            logger.info(f"  ... and {len(json_files) - 3} more")
    
    # SMILES support info
    smiles_support = (
        (root / "scripts" / "ligand_resolver.py").exists() and 
        (root / "config" / "smiles_database.csv").exists()
    )
    if smiles_support:
        try:
            with open(root / "config" / "smiles_database.csv", 'r') as f:
                compound_count = len(f.readlines()) - 1
            logger.info(f"SMILES compounds: {compound_count} available")
        except Exception as e:
            logger.warning(f"Could not read SMILES database: {e}")
            logger.info("SMILES support: Available")
    else:
        logger.info("SMILES support: Not available (CCD-only mode)")
    
    logger.info("=" * 50)

def main():
    """Main workflow function."""
    logger.info("Starting AlphaFold 3 benchmark workflow...")
    
    # Get root directory (parent of scripts/)
    root = pathlib.Path(__file__).resolve().parent.parent
    json_dir = root / "af3_json"
    
    logger.info(f"Working directory: {root}")
    
    # Step 1: Check if template file exists and has data
    if not check_template_file(root):
        logger.error("Template file check failed. Cannot proceed.")
        sys.exit(1)
    
    # Step 2: Check if JSON files exist
    json_files_exist, json_files = check_json_files_exist(json_dir)
    
    # Step 3: Determine if JSON generation/regeneration is needed
    need_generation = False
    
    if not json_files_exist:
        logger.info("No JSON input files found.")
        need_generation = True
    elif should_regenerate_json(root, json_files):
        logger.info("JSON files are outdated.")
        need_generation = True
    else:
        logger.info("✓ JSON input files are up-to-date.")
    
    # Step 4: Generate/regenerate JSON files if needed
    if need_generation:
        if not generate_json_files(root, force_regenerate=json_files_exist):
            logger.error("JSON generation failed. Cannot proceed with benchmark.")
            sys.exit(1)
        
        # Refresh the file list
        json_files_exist, json_files = check_json_files_exist(json_dir)
        
        if not json_files_exist:
            logger.error("JSON generation completed but no valid files were created.")
            logger.error("Please check the template file and try again.")
            sys.exit(1)
        
        logger.info(f"✓ {len(json_files)} JSON files ready for processing.")
    
    # Step 5: Show workflow summary
    show_workflow_summary(root, json_files)
    
    # Step 6: Launch the SLURM job submission process
    logger.info("Launching SLURM job submission...")
    
    try:
        # Use the Python script for SLURM job submission
        remote_job_script = root / "remote_job_submit.py"
        
        if not remote_job_script.exists():
            logger.error(f"SLURM job submission script not found: {remote_job_script}")
            logger.error("Please ensure the script exists and is executable.")
            sys.exit(1)
        
        subprocess.run(
            [sys.executable, str(remote_job_script), str(root)],
            cwd=root,
            check=True
        )
        logger.info("✓ SLURM job submission completed successfully.")
        print('\nJobs submitted to SLURM! After jobs finish, run: parse, generate_excel, plot_curves.')
        
    except subprocess.CalledProcessError as e:
        logger.error(f"✗ SLURM job submission failed with exit code {e.returncode}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"✗ Unexpected error during SLURM job submission: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

