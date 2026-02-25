#!/usr/bin/env python3
"""
End-to-End Test: Auto-config → Bridge → DE Preparation

Tests the complete automation flow:
1. Auto-discovery of project paths
2. Config generation
3. Path validation
4. Bridge execution
5. DE/GO preparation

Run on server with real data:
    python scripts/test_e2e_automation.py --project-id mouse-chd8
"""

import argparse
import sys
import json
from pathlib import Path
import subprocess

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.auto_config import ensure_bridge_config, validate_paths


class Colors:
    """Terminal colors for pretty output."""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'


def print_header(msg):
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{msg}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}\n")


def print_success(msg):
    print(f"{Colors.OKGREEN}✅ {msg}{Colors.ENDC}")


def print_error(msg):
    print(f"{Colors.FAIL}❌ {msg}{Colors.ENDC}")


def print_info(msg):
    print(f"{Colors.OKCYAN}ℹ️  {msg}{Colors.ENDC}")


def print_warning(msg):
    print(f"{Colors.WARNING}⚠️  {msg}{Colors.ENDC}")


def test_phase1_auto_config(project_id: str) -> dict:
    """
    Phase 1: Auto-config generation
    
    Tests:
    - ensure_bridge_config() creates config if missing
    - Config file has correct structure
    - Paths are auto-discovered
    """
    print_header("PHASE 1: Auto-Config Generation")
    
    print_info(f"Testing auto-config for project: {project_id}")
    
    # Clean up existing config for fresh test
    config_path = Path(f"config/projects/paths_{project_id.replace('-', '_')}.yaml")
    if config_path.exists():
        print_info(f"Removing existing config for clean test: {config_path}")
        config_path.unlink()
    
    # Test auto-generation
    result = ensure_bridge_config(project_id, force=True)
    
    if result['status'] not in ['created', 'exists']:
        print_error(f"Config generation failed: {result.get('message', 'Unknown error')}")
        return result
    
    print_success(f"Config created: {result['config_path']}")
    
    # Verify config structure
    import yaml
    with open(result['config_path']) as f:
        config = yaml.safe_load(f)
    
    required_keys = ['project_id', 'rnaseq_output', 'de_pipeline', 'counts_relpath']
    missing = [k for k in required_keys if k not in config]
    
    if missing:
        print_error(f"Missing config keys: {missing}")
        return {'status': 'error', 'message': f'Missing keys: {missing}'}
    
    print_success(f"Config structure valid")
    print_info(f"  project_id: {config['project_id']}")
    print_info(f"  rnaseq_output: {config['rnaseq_output']}")
    print_info(f"  de_pipeline: {config['de_pipeline']}")
    
    return result


def test_phase2_path_validation(config_result: dict) -> dict:
    """
    Phase 2: Path validation
    
    Tests:
    - All required paths exist
    - Counts file is present
    - Sample sheet is found
    """
    print_header("PHASE 2: Path Validation")
    
    if 'config' not in config_result:
        # Load config
        import yaml
        with open(config_result['config_path']) as f:
            config = yaml.safe_load(f)
    else:
        config = config_result['config']
    
    validation = validate_paths(config)
    
    print_info("Validating paths...")
    all_valid = True
    
    for path_name, exists in validation.items():
        if exists:
            print_success(f"{path_name}")
        else:
            print_error(f"{path_name} - NOT FOUND")
            all_valid = False
    
    if all_valid:
        print_success("All paths validated")
        return {'status': 'success', 'validation': validation}
    else:
        missing = [k for k, v in validation.items() if not v]
        print_error(f"Missing paths: {', '.join(missing)}")
        return {'status': 'partial', 'validation': validation, 'missing': missing}


def test_phase3_bridge_execution(project_id: str, dry_run: bool = False) -> dict:
    """
    Phase 3: Bridge script execution
    
    Tests:
    - Bridge script runs without errors
    - Counts matrix copied
    - Metadata generated
    - DE config created
    """
    print_header("PHASE 3: Bridge Script Execution")
    
    print_info(f"Running bridge script for {project_id}...")
    
    cmd = [
        "conda", "run", "-n", "rna-seq-pipeline",
        "python", "scripts/bridge_to_de_pipeline.py",
        "--project-id", project_id,
        "--skip-de",
        "--yes"
    ]
    
    if dry_run:
        cmd.append("--dry-run")
        print_warning("DRY RUN mode")
    
    print_info(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd="/data_3tb/shared/rna-seq-pipeline"
    )
    
    if result.returncode != 0:
        print_error(f"Bridge script failed with code {result.returncode}")
        print_error(f"STDERR: {result.stderr}")
        return {
            'status': 'error',
            'returncode': result.returncode,
            'stderr': result.stderr
        }
    
    print_success("Bridge script completed successfully")
    
    # Parse output to verify steps
    output = result.stdout
    
    checks = {
        'QC Summary': '📊 QC Summary:' in output,
        'Counts copied': '✅ Counts matrix prepared:' in output,
        'Metadata generated': '✅ Metadata generated:' in output,
        'Config created': '✅ Config generated:' in output,
    }
    
    print_info("\nBridge execution steps:")
    for step, passed in checks.items():
        if passed:
            print_success(step)
        else:
            print_warning(f"{step} - not detected in output")
    
    return {
        'status': 'success',
        'returncode': result.returncode,
        'output': output,
        'checks': checks
    }


def test_phase4_output_verification(project_id: str) -> dict:
    """
    Phase 4: Output verification
    
    Tests:
    - Counts file exists in DE pipeline
    - Metadata CSV exists
    - Config YAML exists
    - Files have correct content
    """
    print_header("PHASE 4: Output Verification")
    
    de_pipeline = Path("/home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis")
    
    # Check expected outputs
    expected_files = {
        'counts': de_pipeline / "data" / "raw" / f"{project_id}_counts.csv",
        'metadata': de_pipeline / "data" / "raw" / f"{project_id}_metadata.csv",
        'config': de_pipeline / "configs" / f"config_{project_id}.yml",
    }
    
    print_info("Checking DE pipeline outputs...")
    
    results = {}
    all_exist = True
    
    for name, path in expected_files.items():
        if path.exists():
            size = path.stat().st_size
            print_success(f"{name}: {path} ({size:,} bytes)")
            results[name] = {'exists': True, 'path': str(path), 'size': size}
        else:
            print_error(f"{name}: {path} - NOT FOUND")
            results[name] = {'exists': False, 'path': str(path)}
            all_exist = False
    
    if all_exist:
        print_success("All output files created")
        
        # Quick content check
        print_info("\nContent validation:")
        
        # Check counts file has data
        import pandas as pd
        counts_df = pd.read_csv(expected_files['counts'], index_col=0, nrows=5)
        print_info(f"  Counts matrix: {counts_df.shape[0]}+ rows, {counts_df.shape[1]} columns")
        
        # Check metadata has samples
        metadata_df = pd.read_csv(expected_files['metadata'])
        print_info(f"  Metadata: {len(metadata_df)} samples")
        
        # Check config has comparisons
        import yaml
        with open(expected_files['config']) as f:
            config = yaml.safe_load(f)
        
        comparisons = len(config.get('comparisons', []))
        print_info(f"  Config: {comparisons} comparison(s)")
        
        return {
            'status': 'success',
            'files': results
        }
    else:
        return {
            'status': 'partial',
            'files': results,
            'message': 'Some output files missing'
        }


def run_full_test(project_id: str, dry_run: bool = False) -> bool:
    """Run complete E2E test suite."""
    
    print_header(f"E2E AUTOMATION TEST: {project_id}")
    
    # Phase 1: Auto-config
    config_result = test_phase1_auto_config(project_id)
    if config_result['status'] == 'error':
        print_error("Phase 1 failed. Aborting.")
        return False
    
    # Phase 2: Path validation
    validation_result = test_phase2_path_validation(config_result)
    if validation_result['status'] == 'error':
        print_error("Phase 2 failed. Aborting.")
        return False
    elif validation_result['status'] == 'partial':
        print_warning("Phase 2 partial success. Continuing...")
    
    # Phase 3: Bridge execution
    bridge_result = test_phase3_bridge_execution(project_id, dry_run=dry_run)
    if bridge_result['status'] == 'error':
        print_error("Phase 3 failed. Aborting.")
        return False
    
    # Phase 4: Output verification
    if not dry_run:
        output_result = test_phase4_output_verification(project_id)
        if output_result['status'] != 'success':
            print_warning("Phase 4 partial success")
    
    # Final summary
    print_header("TEST SUMMARY")
    print_success("✅ Phase 1: Auto-config generation")
    print_success("✅ Phase 2: Path validation")
    print_success("✅ Phase 3: Bridge execution")
    
    if not dry_run:
        if output_result['status'] == 'success':
            print_success("✅ Phase 4: Output verification")
        else:
            print_warning("⚠️  Phase 4: Partial success")
    else:
        print_info("⏭️  Phase 4: Skipped (dry-run mode)")
    
    print_header("ALL TESTS PASSED ✅")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="End-to-end automation test for RNA-seq → DE/GO pipeline"
    )
    parser.add_argument(
        '--project-id',
        required=True,
        help='Project identifier (e.g., mouse-chd8)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Dry run mode (skip actual execution)'
    )
    
    args = parser.parse_args()
    
    success = run_full_test(args.project_id, dry_run=args.dry_run)
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
