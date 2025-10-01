#!/usr/bin/env python3
"""
AlphaFold 3 JSON Generator with SMILES Support (Reorganized Version)

This script generates input JSON files for AlphaFold 3 from a template CSV file.
It fetches protein sequences from UniProt, adds ligands and ions as specified,
and ensures all chain IDs are compatible with AlphaFold 3 requirements.

Reorganized structure:
- SMILES database in config/smiles_database.csv
- Ligand resolver in scripts/ligand_resolver.py
- Template in config/af3_template.csv

Key improvements over the original version:
- Supports both CCD codes and SMILES through intelligent ligand resolution
- Maintains full backward compatibility with existing templates
- Uses single uppercase letters for chain IDs (A-Z) as required by AlphaFold 3
- Implements robust error handling and validation
- Follows best practices from AlphaFold 3 documentation
"""

import os
import csv
import json
import string
from typing import List, Dict, Optional, Tuple
import pandas as pd
import requests
from collections import Counter
import shutil
import logging
import sys
import pathlib

# Import the ligand resolver (adjusted for new structure)
try:
    from ligand_resolver import get_resolver, LigandInfo
except ImportError:
    # Fallback if ligand resolver is not available
    logger = logging.getLogger(__name__)
    logger.warning("Ligand resolver not available - only CCD codes will be supported")
    get_resolver = None
    LigandInfo = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# ------------------------------------------------------------------------
# Configuration (adjusted for new structure)
def get_paths():
    """Get paths relative to the project root."""
    # Get root directory (parent of scripts/)
    root = pathlib.Path(__file__).resolve().parent.parent
    
    return {
        'TEMPLATE_PATH': root / "config" / "af3_template.csv",
        'OUTPUT_DIR': root / "af3_json",
        'SMILES_DATABASE_PATH': root / "config" / "smiles_database.csv",
        'ROOT': root
    }

RESULTS_SUFFIX = "_results.csv"

# ------------------------------------------------------------------------
def setup_output_directory(output_dir: pathlib.Path) -> None:
    """
    Sets up the output directory, creating it if it doesn't exist
    and optionally clearing its contents after user confirmation.
    
    Args:
        output_dir: Path to the output directory
    """
    output_dir.mkdir(exist_ok=True)
    
    if output_dir.exists() and any(output_dir.iterdir()):
        choice = input(f"The directory '{output_dir}' is not empty. Do you want to delete its contents? (y/N): ")
        if choice.strip().lower() == 'y':
            for entry in output_dir.iterdir():
                try:
                    if entry.is_dir():
                        shutil.rmtree(entry)
                    else:
                        entry.unlink()
                except Exception as e:
                    logger.error(f"Error deleting {entry}: {e}")
            logger.info(f"Contents of '{output_dir}' have been deleted.")
        else:
            logger.info(f"Contents of '{output_dir}' are kept.")

# ------------------------------------------------------------------------
def fetch_uniprot_json(uniprot_id: str) -> Dict:
    """
    Fetches protein data from UniProt API.
    
    Args:
        uniprot_id: UniProt accession ID
        
    Returns:
        Dictionary containing the protein data
        
    Raises:
        requests.exceptions.RequestException: If the API request fails
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch UniProt data for {uniprot_id}: {e}")
        raise

# ------------------------------------------------------------------------
def parse_uniprot_data(data: Dict) -> Dict:
    """
    Extracts gene name and sequence from UniProt data.
    
    Args:
        data: UniProt JSON data
        
    Returns:
        Dictionary with gene name and sequence
    """
    gene_name = data.get("genes", [{}])[0].get("geneName", {}).get("value", "")
    sequence = data.get("sequence", {}).get("value", "")
    return {"gene": gene_name, "sequence": sequence}

# ------------------------------------------------------------------------
def parse_compound_field(field: str) -> List[str]:
    """
    Parses a compound string like "ATP:2;MG:2" into a list: ["ATP","ATP","MG","MG"]
    
    Args:
        field: Compound string to parse
        
    Returns:
        List of individual compound codes
    """
    items: List[str] = []
    if pd.isna(field) or not field:
        return items
    
    for entry in str(field).split(";"):
        if not entry.strip():
            continue
        try:
            name, count = entry.split(":")
            items.extend([name.strip()] * int(count.strip()))
        except ValueError:
            items.append(entry.strip())
    return items

# ------------------------------------------------------------------------
def create_ligand_sequences_fallback(ligand_list: List[str], ion_list: List[str], 
                                   available_chain_ids: List[str]) -> Tuple[List[Dict], List[str]]:
    """
    Create ligand sequences using traditional CCD codes (fallback when SMILES not available).
    
    Args:
        ligand_list: List of ligand identifiers
        ion_list: List of ion identifiers
        available_chain_ids: List of available chain IDs
        
    Returns:
        Tuple of (sequences list, updated available_chain_ids)
    """
    sequences = []
    
    # Process ligands
    lig_counts = Counter(ligand_list)
    for lig, count in sorted(lig_counts.items()):
        for i in range(count):
            if not available_chain_ids:
                logger.warning(f"Warning: Ran out of available chain IDs. Skipping remaining ligands/ions.")
                break
                
            chain_id = available_chain_ids.pop(0)
            
            sequences.append({
                "ligand": {
                    "id": chain_id,
                    "ccdCodes": [lig.upper()]
                }
            })

    # Process ions
    ion_counts = Counter(ion_list)
    for ion, count in sorted(ion_counts.items()):
        for i in range(count):
            if not available_chain_ids:
                logger.warning(f"Warning: Ran out of available chain IDs. Skipping remaining ions.")
                break
                
            chain_id = available_chain_ids.pop(0)
            
            sequences.append({
                "ligand": {
                    "id": chain_id,
                    "ccdCodes": [ion.upper()]
                }
            })
    
    return sequences, available_chain_ids

# ------------------------------------------------------------------------
def create_ligand_sequences(ligand_list: List[str], ion_list: List[str], 
                          available_chain_ids: List[str], smiles_database_path: pathlib.Path) -> Tuple[List[Dict], List[str]]:
    """
    Create ligand sequences using the ligand resolver (with SMILES support).
    
    Args:
        ligand_list: List of ligand identifiers
        ion_list: List of ion identifiers
        available_chain_ids: List of available chain IDs
        smiles_database_path: Path to the SMILES database
        
    Returns:
        Tuple of (sequences list, updated available_chain_ids)
    """
    sequences = []
    
    # Check if SMILES support is available
    if get_resolver is None or not smiles_database_path.exists():
        logger.info("SMILES support not available, using traditional CCD-only mode")
        return create_ligand_sequences_fallback(ligand_list, ion_list, available_chain_ids)
    
    resolver = get_resolver(str(smiles_database_path))
    
    # Process ligands
    lig_counts = Counter(ligand_list)
    for lig, count in sorted(lig_counts.items()):
        for i in range(count):
            if not available_chain_ids:
                logger.warning(f"Warning: Ran out of available chain IDs. Skipping remaining ligands/ions.")
                break
                
            chain_id = available_chain_ids.pop(0)
            
            # Resolve the ligand
            ligand_info = resolver.resolve_ligand(lig)
            ligand_json = resolver.create_ligand_json(ligand_info, chain_id)
            sequences.append(ligand_json)
            
            # Log the resolution for debugging
            logger.debug(f"Ligand '{lig}' resolved as {ligand_info.ligand_type}: {ligand_info.value}")

    # Process ions (always treated as CCD codes)
    ion_counts = Counter(ion_list)
    for ion, count in sorted(ion_counts.items()):
        for i in range(count):
            if not available_chain_ids:
                logger.warning(f"Warning: Ran out of available chain IDs. Skipping remaining ions.")
                break
                
            chain_id = available_chain_ids.pop(0)
            
            # Ions are always treated as CCD codes
            sequences.append({
                "ligand": {
                    "id": chain_id,
                    "ccdCodes": [ion.upper()]
                }
            })
    
    return sequences, available_chain_ids

# ------------------------------------------------------------------------
def generate_af3_json_files() -> None:
    """
    Generates AlphaFold 3 JSON input files from a template CSV with SMILES support.
    """
    # Get paths for the new structure
    paths = get_paths()
    template_path = paths['TEMPLATE_PATH']
    output_dir = paths['OUTPUT_DIR']
    smiles_database_path = paths['SMILES_DATABASE_PATH']
    
    setup_output_directory(output_dir)

    # Initialize the ligand resolver if available
    if get_resolver is not None and smiles_database_path.exists():
        resolver = get_resolver(str(smiles_database_path))
        stats = resolver.get_database_stats()
        
        if stats['database_loaded']:
            logger.info(f"SMILES database loaded with {stats['total_compounds']} compounds")
            logger.info("Both CCD codes and compound names are supported")
        else:
            logger.warning("SMILES database not loaded - only CCD codes will be supported")
    else:
        logger.info("SMILES support not available - only CCD codes will be supported")

    # Create template if missing
    if not template_path.exists():
        template_path.parent.mkdir(exist_ok=True)
        with open(template_path, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["uniprot_id", "ligands", "ions"])
        logger.info(f"Template CSV created: {template_path}. Please fill it and rerun.")
        return

    # Read input CSV
    try:
        df = pd.read_csv(template_path)
        logger.info(f"Successfully loaded template from {template_path}")
    except Exception as e:
        logger.error(f"Failed to read template CSV: {e}")
        return
        
    records: List[Dict] = []

    for _, row in df.iterrows():
        uniprot_id = str(row.get("uniprot_id", "")).strip()
        if not uniprot_id:
            logger.warning("Skipping row with empty UniProt ID")
            continue
            
        ligand_list = parse_compound_field(row.get("ligands", ""))
        ion_list = parse_compound_field(row.get("ions", ""))

        # Fetch UniProt data
        try:
            data = fetch_uniprot_json(uniprot_id)
            info = parse_uniprot_data(data)
            gene = info.get("gene") or "UnknownGene"
            sequence = info.get("sequence")
            if not sequence:
                raise ValueError("Sequence is empty")
                
            logger.info(f"Successfully fetched data for {uniprot_id} ({gene})")
        except Exception as e:
            logger.error(f"Error for {uniprot_id}: {e}")
            records.append({
                "uniprot_id": uniprot_id,
                "gene": None,
                "status": f"Error: {e}",
                "ligands": ";".join(ligand_list),
                "ions": ";".join(ion_list),
                "json_file": ""
            })
            continue

        # Build sequences list with valid chain IDs (single uppercase letters)
        # Start with 'A' for the protein, then use subsequent letters for ligands and ions
        sequences = [{"protein": {"id": "A", "sequence": sequence}}]
        
        # Get available chain IDs (uppercase letters starting from B)
        available_chain_ids = list(string.ascii_uppercase[1:])  # B, C, D, ...
        
        # Create ligand sequences using the resolver
        ligand_sequences, remaining_chain_ids = create_ligand_sequences(
            ligand_list, ion_list, available_chain_ids, smiles_database_path
        )
        sequences.extend(ligand_sequences)

        # Create JSON object
        json_obj = {
            "name": gene,
            "sequences": sequences,
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 2  # Using version 2 to support external MSA paths
        }

        # Build filename with quantities (preserving original naming logic)
        lig_counts = Counter(ligand_list)
        ion_counts = Counter(ion_list)
        lig_suffix = "".join(f"_{lig}{count}" for lig, count in sorted(lig_counts.items()))
        ion_suffix = "".join(f"_{ion}{count}" for ion, count in sorted(ion_counts.items()))
        filename = f"{gene}_{uniprot_id}{lig_suffix}{ion_suffix}.json"
        filepath = output_dir / filename

        # Write JSON file
        try:
            with open(filepath, 'w') as jf:
                json.dump(json_obj, jf, indent=2)
            logger.info(f"JSON created: {filename}")
        except Exception as e:
            logger.error(f"Failed to write JSON file {filepath}: {e}")
            continue

        records.append({
            "uniprot_id": uniprot_id,
            "gene": gene,
            "status": "OK",
            "ligands": ";".join(ligand_list),
            "ions": ";".join(ion_list),
            "json_file": filename
        })

    # Write summary CSV
    if records:
        summary_path = template_path.with_suffix('').with_suffix('.csv').with_name(
            template_path.stem + RESULTS_SUFFIX
        )
        try:
            pd.DataFrame(records).to_csv(summary_path, index=False)
            logger.info(f"Summary of results written to: {summary_path}")
        except Exception as e:
            logger.error(f"Failed to write summary CSV: {e}")
    else:
        logger.warning("No records were processed successfully")

# ------------------------------------------------------------------------

if __name__ == "__main__":
    generate_af3_json_files()

# ------------------------------------------------------------------------
# Usage Examples:
#
# Template CSV with mixed CCD codes and compound names:
# uniprot_id,ligands,ions
# P08183,ATP:2;Diclofenac:1,MG:2
# O15245,Ketoprofen:1;NAD:1,CA:1
#
# In this example:
# - ATP, NAD are CCD codes (will use ccdCodes format)
# - Diclofenac, Ketoprofen are compound names (will use SMILES format if available)
# - MG, CA are ion CCD codes (will use ccdCodes format)
#
# The script automatically detects the type and generates appropriate JSON structures.
