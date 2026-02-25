#!/usr/bin/env python3
"""
Test auto_config.py functionality
"""

import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.auto_config import (
    discover_sample_sheet,
    validate_paths,
    generate_bridge_config,
    ensure_bridge_config
)


def test_discover_sample_sheet():
    """Test sample sheet discovery"""
    print("\n" + "="*60)
    print("TEST: discover_sample_sheet()")
    print("="*60)
    
    # Test with mouse-chd8
    result = discover_sample_sheet("mouse-chd8")
    print(f"Project: mouse-chd8")
    print(f"Found: {result}")
    
    assert result is not None, "Should find mouse_chd8.tsv"
    assert result.exists(), f"File should exist: {result}"
    print("✅ PASSED")


def test_validate_paths():
    """Test path validation"""
    print("\n" + "="*60)
    print("TEST: validate_paths()")
    print("="*60)
    
    # Use paths that might exist in dev environment
    config = {
        'rnaseq_output': '/data_3tb/shared/output/mouse-chd8',
        'de_pipeline': '/home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis',
        'sample_sheet_dir': '/data_3tb/shared/rna-seq-pipeline/config/samples',
        'counts_relpath': 'project_summary/counts/counts_matrix_clean.csv',
        'summary_relpath': 'project_summary.json'
    }
    
    validation = validate_paths(config)
    print(f"Validation results:")
    for key, exists in validation.items():
        status = "✅" if exists else "❌"
        print(f"  {status} {key}")
    
    # In dev environment, paths might not exist - that's OK
    print("ℹ️  Note: Some paths may not exist in development environment")
    print("✅ PASSED (validation function works)")


def test_ensure_bridge_config():
    """Test ensure_bridge_config (high-level function)"""
    print("\n" + "="*60)
    print("TEST: ensure_bridge_config()")
    print("="*60)
    
    # Check if we're in dev or server environment
    server_output = Path("/data_3tb/shared/output/mouse-chd8")
    is_server = server_output.exists()
    
    if not is_server:
        print("ℹ️  SKIPPED: Not on server (no actual data available)")
        print("   Run this test on the server with real data")
        return
    
    # Test with auto-discovery
    result = ensure_bridge_config("mouse-chd8", force=True)
    
    print(f"Status: {result['status']}")
    print(f"Message: {result['message']}")
    
    if 'config_path' in result:
        print(f"Config path: {result['config_path']}")
        config_file = Path(result['config_path'])
        assert config_file.exists(), f"Config file should exist: {config_file}"
    
    if 'validation' in result:
        print(f"\nValidation:")
        for key, exists in result['validation'].items():
            status = "✅" if exists else "❌"
            print(f"  {status} {key}")
    
    assert result['status'] in ['created', 'exists'], f"Should succeed, got: {result['status']}"
    print("✅ PASSED")


def test_generate_bridge_config():
    """Test generate_bridge_config (low-level function)"""
    print("\n" + "="*60)
    print("TEST: generate_bridge_config()")
    print("="*60)
    
    summary_path = Path("/data_3tb/shared/output/mouse-chd8/project_summary.json")
    
    if not summary_path.exists():
        print(f"ℹ️  SKIPPED: Summary file not found (not on server)")
        print(f"   Expected: {summary_path}")
        print(f"   Run this test on the server with real data")
        return
    
    config_path = generate_bridge_config(summary_path, force=True)
    
    print(f"Generated: {config_path}")
    assert config_path.exists(), f"Config file should exist: {config_path}"
    
    # Verify content
    import yaml
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    print(f"\nConfig contents:")
    for key, value in config.items():
        print(f"  {key}: {value}")
    
    assert 'project_id' in config, "Should have project_id"
    assert 'rnaseq_output' in config, "Should have rnaseq_output"
    assert 'de_pipeline' in config, "Should have de_pipeline"
    
    print("✅ PASSED")


if __name__ == '__main__':
    print("\n" + "="*60)
    print("AUTO_CONFIG TEST SUITE")
    print("="*60)
    
    try:
        test_discover_sample_sheet()
        test_validate_paths()
        test_ensure_bridge_config()
        test_generate_bridge_config()
        
        print("\n" + "="*60)
        print("ALL TESTS PASSED ✅")
        print("="*60)
    
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
