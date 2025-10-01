#!/usr/bin/env python3
"""
Submit Benchmark Script - Updated for SLURM Workflow

This script serves as the main entry point for the AlphaFold 3 benchmark workflow.
It calls the intelligent main.py which handles:
1. Automatic JSON generation (with SMILES support)
2. Change detection and regeneration
3. SLURM job submission with MSA fusion support

The main.py script manages the complete workflow automatically.
"""

import subprocess
import sys
import pathlib
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def main():
    """Launch the intelligent workflow via main.py"""
    # Get project root (parent of scripts/)
    root = pathlib.Path(__file__).resolve().parent.parent
    
    logger.info(f"Starting AlphaFold 3 benchmark workflow from {root}")
    
    # Call the intelligent main.py which handles everything
    try:
        subprocess.run(
            [sys.executable, str(root / 'scripts' / 'main.py')], 
            check=True,
            cwd=root
        )
        logger.info("Workflow completed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Workflow failed with exit code {e.returncode}")
        sys.exit(e.returncode)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

