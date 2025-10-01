#!/usr/bin/env python3
"""
Ligand Resolver Module for AlphaFold 3 JSON Generator (Reorganized Version)

This module provides intelligent ligand resolution capabilities, automatically
distinguishing between CCD codes and compound names, and resolving compound
names to their corresponding SMILES strings.

Reorganized structure:
- SMILES database in config/smiles_database.csv
- Module in scripts/ligand_resolver.py

Key features:
- Automatic detection of CCD codes vs compound names
- SMILES database lookup for pharmaceutical compounds
- Fallback to CCD codes for unknown compounds
- Full compatibility with existing workflows
"""

import os
import re
import pandas as pd
import logging
import pathlib
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class LigandInfo:
    """Information about a resolved ligand."""
    identifier: str
    ligand_type: str  # 'ccd' or 'smiles'
    value: str  # CCD code or SMILES string
    original_name: str

class LigandResolver:
    """
    Intelligent ligand resolver that distinguishes between CCD codes and compound names,
    and resolves compound names to SMILES strings.
    """
    
    def __init__(self, smiles_database_path: str = None):
        """
        Initialize the ligand resolver.
        
        Args:
            smiles_database_path: Path to the SMILES database CSV file
                                If None, will use config/smiles_database.csv relative to project root
        """
        if smiles_database_path is None:
            # Get project root (parent of scripts/)
            root = pathlib.Path(__file__).resolve().parent.parent
            smiles_database_path = str(root / "config" / "smiles_database.csv")
        
        self.smiles_database_path = smiles_database_path
        self.smiles_db: Dict[str, str] = {}
        self.compound_names: set = set()
        self._load_smiles_database()
        
        # Common CCD codes for quick recognition
        self.common_ccd_codes = {
            'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'NAD', 'NADH', 'FAD', 'FMN',
            'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'K', 'NA', 'CL', 'SO4', 'PO4',
            'HOH', 'WAT', 'GOL', 'EDO', 'PEG', 'DMS', 'ACE', 'NH4'
        }
    
    def _load_smiles_database(self) -> None:
        """Load the SMILES database from CSV file."""
        if not os.path.exists(self.smiles_database_path):
            logger.warning(f"SMILES database not found at {self.smiles_database_path}")
            logger.warning("Only CCD codes will be supported.")
            return
        
        try:
            df = pd.read_csv(self.smiles_database_path)
            
            # Validate required columns
            required_columns = ['compound_name', 'smiles']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(f"Missing required columns in SMILES database: {missing_columns}")
            
            # Build the lookup dictionary (case-insensitive)
            for _, row in df.iterrows():
                compound_name = str(row['compound_name']).strip()
                smiles = str(row['smiles']).strip()
                
                if compound_name and smiles and smiles != 'nan':
                    # Store both original case and lowercase for flexible lookup
                    self.smiles_db[compound_name] = smiles
                    self.smiles_db[compound_name.lower()] = smiles
                    self.compound_names.add(compound_name.lower())
            
            logger.info(f"Loaded {len(self.compound_names)} compounds from SMILES database")
            
        except Exception as e:
            logger.error(f"Error loading SMILES database: {e}")
            logger.warning("Only CCD codes will be supported.")
    
    def _is_likely_ccd_code(self, identifier: str) -> bool:
        """
        Determine if an identifier is likely a CCD code based on patterns.
        
        Args:
            identifier: The identifier to check
            
        Returns:
            True if likely a CCD code, False otherwise
        """
        identifier = identifier.strip()
        
        # Check if it's a known common CCD code
        if identifier.upper() in self.common_ccd_codes:
            return True
        
        # CCD code patterns:
        # - Short (typically 2-6 characters)
        # - Often all uppercase
        # - May contain numbers
        # - May contain hyphens (e.g., LIG-1337)
        
        # Very short identifiers are likely CCD codes
        if len(identifier) <= 3:
            return True
        
        # All uppercase with possible numbers/hyphens
        if re.match(r'^[A-Z0-9\-]+$', identifier):
            return True
        
        # Pattern like LIG-#### or similar
        if re.match(r'^[A-Z]{2,4}-\d+$', identifier):
            return True
        
        return False
    
    def _is_known_compound(self, identifier: str) -> bool:
        """
        Check if an identifier is a known compound in the SMILES database.
        
        Args:
            identifier: The identifier to check
            
        Returns:
            True if it's a known compound, False otherwise
        """
        return identifier.lower() in self.compound_names
    
    def resolve_ligand(self, identifier: str) -> LigandInfo:
        """
        Resolve a ligand identifier to either a CCD code or SMILES string.
        
        Args:
            identifier: The ligand identifier to resolve
            
        Returns:
            LigandInfo object with resolution details
        """
        identifier = identifier.strip()
        
        # First, check if it's a known compound in our SMILES database
        if self._is_known_compound(identifier):
            smiles = self.smiles_db[identifier.lower()]
            logger.debug(f"Resolved '{identifier}' to SMILES: {smiles}")
            return LigandInfo(
                identifier=identifier,
                ligand_type='smiles',
                value=smiles,
                original_name=identifier
            )
        
        # If not in SMILES database, check if it looks like a CCD code
        if self._is_likely_ccd_code(identifier):
            logger.debug(f"Treating '{identifier}' as CCD code")
            return LigandInfo(
                identifier=identifier,
                ligand_type='ccd',
                value=identifier.upper(),  # CCD codes are typically uppercase
                original_name=identifier
            )
        
        # Fallback: treat as CCD code (for compatibility)
        logger.warning(f"Unknown identifier '{identifier}', treating as CCD code")
        return LigandInfo(
            identifier=identifier,
            ligand_type='ccd',
            value=identifier,
            original_name=identifier
        )
    
    def resolve_ligand_list(self, ligand_list: List[str]) -> List[LigandInfo]:
        """
        Resolve a list of ligand identifiers.
        
        Args:
            ligand_list: List of ligand identifiers
            
        Returns:
            List of LigandInfo objects
        """
        return [self.resolve_ligand(ligand) for ligand in ligand_list]
    
    def create_ligand_json(self, ligand_info: LigandInfo, chain_id: str) -> Dict:
        """
        Create the JSON structure for a ligand according to AlphaFold3 format.
        
        Args:
            ligand_info: Resolved ligand information
            chain_id: Chain identifier for the ligand
            
        Returns:
            Dictionary representing the ligand in AlphaFold3 JSON format
        """
        if ligand_info.ligand_type == 'smiles':
            return {
                "ligand": {
                    "id": chain_id,
                    "smiles": ligand_info.value
                }
            }
        else:  # CCD code
            return {
                "ligand": {
                    "id": chain_id,
                    "ccdCodes": [ligand_info.value]
                }
            }
    
    def get_database_stats(self) -> Dict[str, int]:
        """
        Get statistics about the loaded SMILES database.
        
        Returns:
            Dictionary with database statistics
        """
        return {
            "total_compounds": len(self.compound_names),
            "database_loaded": len(self.smiles_db) > 0
        }
    
    def list_available_compounds(self, limit: int = 20) -> List[str]:
        """
        List available compounds in the database.
        
        Args:
            limit: Maximum number of compounds to return
            
        Returns:
            List of compound names
        """
        compounds = sorted(list(self.compound_names))
        return compounds[:limit]

# Global resolver instance (initialized when module is imported)
_resolver_instance: Optional[LigandResolver] = None

def get_resolver(smiles_database_path: str = None) -> LigandResolver:
    """
    Get the global ligand resolver instance.
    
    Args:
        smiles_database_path: Path to the SMILES database CSV file
                            If None, will use config/smiles_database.csv
        
    Returns:
        LigandResolver instance
    """
    global _resolver_instance
    if _resolver_instance is None:
        _resolver_instance = LigandResolver(smiles_database_path)
    return _resolver_instance

def resolve_ligand(identifier: str) -> LigandInfo:
    """
    Convenience function to resolve a single ligand.
    
    Args:
        identifier: The ligand identifier to resolve
        
    Returns:
        LigandInfo object with resolution details
    """
    resolver = get_resolver()
    return resolver.resolve_ligand(identifier)

def resolve_ligand_list(ligand_list: List[str]) -> List[LigandInfo]:
    """
    Convenience function to resolve a list of ligands.
    
    Args:
        ligand_list: List of ligand identifiers
        
    Returns:
        List of LigandInfo objects
    """
    resolver = get_resolver()
    return resolver.resolve_ligand_list(ligand_list)

if __name__ == "__main__":
    # Test the resolver
    resolver = LigandResolver()
    
    # Test some examples
    test_ligands = ["ATP", "Metformin", "MG", "Ketoprofen", "NAD", "UnknownCompound"]
    
    print("Ligand Resolution Test:")
    print("=" * 50)
    
    for ligand in test_ligands:
        info = resolver.resolve_ligand(ligand)
        print(f"'{ligand}' -> {info.ligand_type}: {info.value}")
    
    print("\nDatabase Statistics:")
    stats = resolver.get_database_stats()
    print(f"Total compounds: {stats['total_compounds']}")
    print(f"Database loaded: {stats['database_loaded']}")
    
    if stats['database_loaded']:
        print("\nFirst 10 available compounds:")
        compounds = resolver.list_available_compounds(10)
        for compound in compounds:
            print(f"  - {compound}")


