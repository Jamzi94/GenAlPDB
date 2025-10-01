
#!/usr/bin/env python3
"""
AlphaFold 3 JSON Fusion Tool

This script intelligently merges MSA data from AlphaFold 3 data pipeline into original JSON files
while ensuring chain IDs remain compatible with AlphaFold 3 requirements.

Key improvements:
- Validates and fixes chain IDs to be single uppercase letters (A-Z) as required by AlphaFold 3
- Implements robust error handling and validation
- Uses deep recursive merging to preserve all data
- Follows best practices from AlphaFold 3 documentation
"""

import json
import sys
import os
import argparse
import string
import logging
from datetime import datetime
from typing import Dict, Tuple, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def deep_merge_protein_data(original_protein: Dict[str, Any], msa_protein: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merges all fields from MSA protein into original protein.
    Preserves original fields if they don't exist in MSA protein.
    
    This function performs a deep merge of protein data structures, combining
    MSA (Multiple Sequence Alignment) data with the original protein data.
    The MSA data typically contains additional fields like unpaired_msa,
    paired_msa, and templates that enhance the prediction quality.
    
    Args:
        original_protein: Protein data from original JSON (may include ligands)
        msa_protein: Protein data from MSA JSON with additional fields
        
    Returns:
        Merged protein data combining both sources
    """
    # Start with a copy of the original protein data
    merged = original_protein.copy()
    
    # Add or update fields from MSA protein
    for key, value in msa_protein.items():
        # If the field exists in both and both are dictionaries, merge recursively
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = deep_merge_protein_data(merged[key], value)
        # Otherwise, use the MSA value (overwrite existing or add new)
        else:
            merged[key] = value
            
    return merged

def fix_chain_ids(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fixes chain IDs to be compatible with AlphaFold 3 (single uppercase letters).
    
    AlphaFold 3 requires chain IDs to be single uppercase letters (A-Z).
    This function validates and normalizes chain IDs in the JSON data structure,
    ensuring compliance with AlphaFold 3 requirements.
    
    Args:
        data: The JSON data to fix
        
    Returns:
        The fixed JSON data with valid chain IDs (A-Z)
        
    Raises:
        ValueError: If there are more than 26 chains (exceeds A-Z limit)
    """
    if 'sequences' not in data:
        return data
    
    # Get available chain IDs (uppercase letters A through Z)
    available_chain_ids = list(string.ascii_uppercase)  # A, B, C, D, ..., Z
    
    # Create a mapping to track old ID to new ID conversions
    id_mapping = {}
    
    # First pass: collect all existing IDs and assign new valid IDs
    for seq in data['sequences']:
        if 'protein' in seq:
            old_id = seq['protein']['id']
            # Protein chains usually get 'A' by convention
            new_id = 'A'
            if new_id in id_mapping.values():
                # If 'A' is already taken, get the next available letter
                if available_chain_ids:
                    new_id = available_chain_ids.pop(0)
                else:
                    raise ValueError("Ran out of available chain IDs (max 26 chains: A-Z)")
            id_mapping[old_id] = new_id
            if new_id in available_chain_ids:
                available_chain_ids.remove(new_id)
        elif 'ligand' in seq:
            old_id = seq['ligand']['id']
            if old_id not in id_mapping:
                if available_chain_ids:
                    new_id = available_chain_ids.pop(0)
                    id_mapping[old_id] = new_id
                else:
                    raise ValueError("Ran out of available chain IDs (max 26 chains: A-Z)")
    
    # Second pass: update all IDs in the sequences
    for seq in data['sequences']:
        if 'protein' in seq:
            old_id = seq['protein']['id']
            seq['protein']['id'] = id_mapping[old_id]
        elif 'ligand' in seq:
            old_id = seq['ligand']['id']
            seq['ligand']['id'] = id_mapping[old_id]
    
    return data

def is_valid_json_file(file_path: str, required_keys: List[str] = None) -> bool:
    """
    Validates if a file is a valid JSON file containing required keys.
    
    Args:
        file_path: Path to the JSON file
        required_keys: List of keys that must be present in the JSON
        
    Returns:
        True if valid, False otherwise
    """
    if not os.path.isfile(file_path):
        logger.error(f"File does not exist: {file_path}")
        return False
        
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            
        if required_keys:
            for key in required_keys:
                if key not in data:
                    logger.error(f"Required key '{key}' not found in {file_path}")
                    return False
        return True
    except json.JSONDecodeError:
        logger.error(f"Invalid JSON format in {file_path}")
        return False
    except Exception as e:
        logger.error(f"Error validating {file_path}: {e}")
        return False

def fusion_json_msa(json_original_path: str, json_msa_path: str) -> bool:
    """
    Intelligently merges all data from MSA JSON into the original JSON,
    while ensuring chain IDs are compatible with AlphaFold 3.
    
    This function combines the original JSON (which may contain ligands and ions)
    with MSA data generated by the AlphaFold 3 data pipeline. The MSA data includes
    multiple sequence alignments that improve prediction quality. The merged file
    is written back to the original path, preserving ligand information.
    
    Args:
        json_original_path: Path to the original JSON file (may contain ligands)
        json_msa_path: Path to the MSA JSON file (contains MSA fields)
        
    Returns:
        True if successful, False otherwise
        
    Raises:
        ValueError: If input files are invalid or missing required sections
    """
    # Validate input files exist and have required structure
    if not is_valid_json_file(json_original_path, ['sequences']):
        raise ValueError(f"Original JSON file is invalid or missing required sections: {json_original_path}")
    
    if not is_valid_json_file(json_msa_path, ['sequences']):
        raise ValueError(f"MSA JSON file is invalid or missing required sections: {json_msa_path}")
    
    # Skip if the file is a mapping file (not a real input JSON)
    if "_chain_mapping.json" in json_original_path:
        logger.warning(f"Skipping chain mapping file (not a valid input JSON): {json_original_path}")
        return False
    
    # Load JSON files
    with open(json_original_path, 'r') as f_orig:
        data_original = json.load(f_orig)
    
    with open(json_msa_path, 'r') as f_msa:
        data_msa = json.load(f_msa)
    
    # Fix chain IDs in original data to ensure AlphaFold 3 compatibility
    data_original = fix_chain_ids(data_original)
    
    # Find protein sequence in the original file
    protein_original = None
    protein_index = -1
    
    for i, seq in enumerate(data_original['sequences']):
        if 'protein' in seq:
            protein_original = seq['protein']
            protein_index = i
            break
    
    if protein_original is None:
        raise ValueError(f"No protein sequence found in original JSON file: {json_original_path}")
    
    # Find protein sequence with MSA data in the MSA file
    protein_msa = None
    
    for seq in data_msa['sequences']:
        if 'protein' in seq:
            protein_msa = seq['protein']
            break
    
    if protein_msa is None:
        raise ValueError(f"No protein sequence found in MSA file: {json_msa_path}")
    
    # Merge all protein data intelligently (combines MSA fields with original)
    merged_protein = deep_merge_protein_data(protein_original, protein_msa)
    
    # Replace the original protein data with the merged data
    data_original['sequences'][protein_index]['protein'] = merged_protein
    
    # Update top-level fields from MSA if they don't override original fields
    # (except for sequences which we've already handled)
    for key, value in data_msa.items():
        if key != 'sequences':
            # If the field exists in both, prefer MSA version as it's likely more up-to-date
            # (e.g., version number for AlphaFold3)
            data_original[key] = value
    
    # Ensure version is set to at least 2 (required for MSA data)
    if 'version' in data_original and data_original['version'] < 2:
        data_original['version'] = 2
    
    # Write the merged file directly back to the original path
    with open(json_original_path, 'w') as f_out:
        json.dump(data_original, f_out, indent=2)
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Intelligently merges MSA data into original JSON files with AlphaFold 3 compatible chain IDs.'
    )
    parser.add_argument('json_original', help='Path to the original JSON file (potentially containing ligands)')
    parser.add_argument('json_msa', help='Path to the JSON file generated by MSA (containing additional fields)')
    
    args = parser.parse_args()
    
    try:
        success = fusion_json_msa(args.json_original, args.json_msa)
        if success:
            logger.info(f"Complete fusion successful! Original file modified: {args.json_original}")
    except Exception as e:
        logger.error(f"Error during fusion: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

