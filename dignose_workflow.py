#!/usr/bin/env python3
"""
AlphaFold 3 Workflow Diagnostic Tool

This script performs comprehensive diagnostics on the AlphaFold 3 workflow setup,
identifying common issues and providing solutions.

Key features:
- Validates directory structure and required files
- Checks configuration files and templates
- Verifies script dependencies and consistency
- Analyzes JSON file status and validity
- Provides actionable recommendations for fixing issues
"""

import pathlib
import sys
import os
import json
import logging
from typing import List, Dict, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def check_file_exists(file_path: pathlib.Path, description: str) -> bool:
    """
    Check if a file exists and report status.
    
    Args:
        file_path: Path to the file to check
        description: Description of the file for reporting
        
    Returns:
        True if file exists, False otherwise
    """
    if file_path.exists():
        print(f"✓ {description}: {file_path}")
        return True
    else:
        print(f"✗ {description}: {file_path} (MISSING)")
        return False

def check_json_validity(json_dir: pathlib.Path) -> Tuple[int, int]:
    """
    Check the validity of JSON files in the specified directory.
    
    This function scans the JSON directory and validates each file to ensure
    it contains the required structure for AlphaFold 3 input. Valid files
    must have 'sequences' and 'name' fields.
    
    Args:
        json_dir: Path to the directory containing JSON files
        
    Returns:
        Tuple of (total_files, valid_files)
    """
    if not json_dir.exists():
        return 0, 0
    
    total_files = 0
    valid_files = 0
    
    for json_file in json_dir.glob("*.json"):
        # Skip auxiliary files like chain mapping files
        if "_chain_mapping.json" in json_file.name:
            continue
            
        total_files += 1
        
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Check for required keys in AlphaFold 3 JSON format
            if 'sequences' in data and 'name' in data:
                valid_files += 1
        except:
            # File is invalid JSON or missing required fields
            pass
    
    return total_files, valid_files

def check_script_consistency(scripts_dir: pathlib.Path) -> Dict[str, bool]:
    """
    Check the consistency of script files.
    
    This function validates that workflow scripts are correctly configured
    and reference the appropriate files. It checks for common misconfigurations
    such as calling .sh files instead of .py files.
    
    Args:
        scripts_dir: Path to the scripts directory
        
    Returns:
        Dictionary of script checks and their results
    """
    results = {}
    
    # Check if main.py calls remote_job_submit.py (not .sh)
    main_py = scripts_dir / "main.py"
    if main_py.exists():
        try:
            with open(main_py, 'r') as f:
                content = f.read()
            
            if "remote_job_submit.py" in content:
                results["main_py_calls_correct_script"] = True
            elif "remote_job_submit.sh" in content:
                results["main_py_calls_correct_script"] = False
            else:
                results["main_py_calls_correct_script"] = None
        except:
            results["main_py_calls_correct_script"] = None
    
    # Check if submit_benchmark.py calls main.py
    submit_py = scripts_dir / "submit_benchmark.py"
    if submit_py.exists():
        try:
            with open(submit_py, 'r') as f:
                content = f.read()
            
            if "main.py" in content:
                results["submit_py_calls_main"] = True
            else:
                results["submit_py_calls_main"] = False
        except:
            results["submit_py_calls_main"] = None
    
    return results

def check_workflow_setup():
    """
    Check the current workflow setup and provide diagnostic information.
    
    This function performs a comprehensive diagnostic of the AlphaFold 3 workflow,
    checking directory structure, configuration files, scripts, and identifying
    common issues that may prevent the workflow from running correctly.
    
    It provides actionable recommendations for fixing any detected problems.
    """
    print("=" * 60)
    print("AlphaFold 3 Workflow Diagnostic")
    print("=" * 60)
    
    # Get project root
    if len(sys.argv) > 1:
        root = pathlib.Path(sys.argv[1]).resolve()
    else:
        root = pathlib.Path.cwd()
    
    print(f"Project root: {root}")
    print()
    
    # Check directory structure
    print("1. Directory Structure:")
    print("-" * 30)
    
    required_dirs = [
        (root / "config", "Config directory"),
        (root / "scripts", "Scripts directory"),
        (root / "af3_json", "JSON output directory")
    ]
    
    for dir_path, desc in required_dirs:
        check_file_exists(dir_path, desc)
    
    print()
    
    # Check configuration files
    print("2. Configuration Files:")
    print("-" * 30)
    
    config_files = [
        (root / "config" / "af3_template.csv", "Template CSV"),
        (root / "config" / "benchmark_config.csv", "Benchmark config"),
        (root / "config" / "smiles_database.csv", "SMILES database (optional)")
    ]
    
    for file_path, desc in config_files:
        check_file_exists(file_path, desc)
    
    print()
    
    # Check script files
    print("3. Script Files:")
    print("-" * 30)
    
    script_files = [
        (root / "scripts" / "main.py", "Main workflow script"),
        (root / "scripts" / "generate_af3_json.py", "JSON generator"),
        (root / "scripts" / "ligand_resolver.py", "Ligand resolver (optional)"),
        (root / "scripts" / "submit_benchmark.py", "Submit benchmark script"),
        (root / "scripts" / "fusion_json.py", "Fusion script"),
        (root / "remote_job_submit.py", "Remote job submit script")
    ]
    
    for file_path, desc in script_files:
        check_file_exists(file_path, desc)
    
    print()
    
    # Check template content
    print("4. Template Content:")
    print("-" * 30)
    
    template_path = root / "config" / "af3_template.csv"
    if template_path.exists():
        try:
            with open(template_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) > 1:
                print(f"✓ Template has {len(lines) - 1} data rows")
                print("  Sample entries:")
                for i, line in enumerate(lines[:3]):
                    print(f"    {i}: {line.strip()}")
            else:
                print("✗ Template has no data rows")
        except Exception as e:
            print(f"✗ Error reading template: {e}")
    
    print()
    
    # Check script consistency
    print("5. Script Consistency Check:")
    print("-" * 30)
    
    consistency_results = check_script_consistency(root / "scripts")
    
    if "main_py_calls_correct_script" in consistency_results:
        if consistency_results["main_py_calls_correct_script"] is True:
            print("✓ main.py correctly calls remote_job_submit.py")
        elif consistency_results["main_py_calls_correct_script"] is False:
            print("✗ main.py incorrectly calls remote_job_submit.sh (should be .py)")
        else:
            print("? Could not determine if main.py calls correct script")
    
    if "submit_py_calls_main" in consistency_results:
        if consistency_results["submit_py_calls_main"] is True:
            print("✓ submit_benchmark.py calls main.py (NEW WORKFLOW)")
        elif consistency_results["submit_py_calls_main"] is False:
            print("✗ submit_benchmark.py does not call main.py (OLD WORKFLOW)")
        else:
            print("? Could not determine if submit_benchmark.py calls main.py")
    
    print()
    
    # Check JSON files
    print("6. JSON Files Status:")
    print("-" * 30)
    
    json_dir = root / "af3_json"
    if json_dir.exists():
        total_files, valid_files = check_json_validity(json_dir)
        
        print(f"Total JSON files: {total_files}")
        print(f"Valid input files: {valid_files}")
        
        if valid_files > 0:
            print("✓ Valid JSON input files found")
            json_files = list(json_dir.glob("*.json"))
            valid_json_files = [f for f in json_files if "_chain_mapping.json" not in f.name]
            print("Sample files:")
            for f in valid_json_files[:3]:
                print(f"  - {f.name}")
        else:
            print("✗ No valid JSON input files found")
    else:
        print("✗ JSON directory does not exist")
    
    print()
    
    # Check for remote_job_submit.sh vs .py
    print("7. Remote Job Submit Script Check:")
    print("-" * 30)
    
    remote_py = root / "remote_job_submit.py"
    remote_sh = root / "remote_job_submit.sh"
    
    if remote_py.exists() and not remote_sh.exists():
        print("✓ Using correct remote_job_submit.py script")
    elif remote_sh.exists() and not remote_py.exists():
        print("✗ Using incorrect remote_job_submit.sh script (should be .py)")
    elif remote_py.exists() and remote_sh.exists():
        print("! Both remote_job_submit.py and .sh exist (potential conflict)")
    else:
        print("✗ No remote job submit script found")
    
    print()
    
    # Recommendations
    print("8. Recommendations:")
    print("-" * 30)
    
    issues_found = False
    
    # Check for remote_job_submit.sh issue
    if remote_sh.exists() and not remote_py.exists():
        issues_found = True
        print("! CRITICAL ISSUE: remote_job_submit.sh should be renamed to remote_job_submit.py")
        print("  This is causing the workflow to fail.")
        print("  Solution: Rename the file and update references in main.py")
    
    # Check if main.py calls remote_job_submit.sh
    if "main_py_calls_correct_script" in consistency_results and consistency_results["main_py_calls_correct_script"] is False:
        issues_found = True
        print("! CRITICAL ISSUE: main.py is calling remote_job_submit.sh instead of .py")
        print("  Solution: Update main.py to call remote_job_submit.py instead")
    
    # Check if template is empty
    template_path = root / "config" / "af3_template.csv"
    if template_path.exists():
        try:
            with open(template_path, 'r') as f:
                lines = f.readlines()
            if len(lines) <= 1:
                issues_found = True
                print("! Template file is empty - add protein entries")
        except:
            pass
    
    if not issues_found:
        print("✓ No obvious issues detected")
        print("  Try running: python3 scripts/main.py")

if __name__ == "__main__":
    check_workflow_setup()

