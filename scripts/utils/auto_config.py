#!/usr/bin/env python3
"""
Auto-generate bridge configuration from project summary.

This module automatically discovers paths and generates bridge config files
for connecting RNA-seq results to DE/GO analysis pipeline.
"""

import json
import yaml
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List


def discover_sample_sheet(project_id: str, search_dirs: List[Path] = None) -> Optional[Path]:
    """
    Discover sample sheet location by trying common locations.
    
    Args:
        project_id: Project identifier (e.g., "mouse-chd8" or "mouse_chd8")
        search_dirs: Optional list of directories to search
    
    Returns:
        Path to sample sheet if found, None otherwise
    """
    if search_dirs is None:
        search_dirs = [
            Path("/data_3tb/shared/rna-seq-pipeline/config/samples"),
            Path("/home/ygkim/ngs-pipeline/rna-seq-pipeline/config/samples"),
            Path.cwd() / "config" / "samples",
        ]
    
    # Try both hyphen and underscore versions
    project_id_variants = [
        project_id,
        project_id.replace('-', '_'),
        project_id.replace('_', '-')
    ]
    
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
            
        for variant in project_id_variants:
            candidate = search_dir / f"{variant}.tsv"
            if candidate.exists():
                return candidate
    
    return None


def validate_paths(config: Dict[str, str]) -> Dict[str, bool]:
    """
    Validate that all required paths exist.
    
    Args:
        config: Configuration dictionary with paths
    
    Returns:
        Dictionary mapping path names to existence status
    """
    validation = {}
    
    # Check rnaseq_output
    if 'rnaseq_output' in config:
        rnaseq_path = Path(config['rnaseq_output'])
        validation['rnaseq_output'] = rnaseq_path.exists()
        
        # Check counts file
        if 'counts_relpath' in config:
            counts_path = rnaseq_path / config['counts_relpath']
            validation['counts_file'] = counts_path.exists()
        
        # Check summary file
        if 'summary_relpath' in config:
            summary_path = rnaseq_path / config['summary_relpath']
            validation['summary_file'] = summary_path.exists()
    
    # Check de_pipeline
    if 'de_pipeline' in config:
        de_path = Path(config['de_pipeline'])
        validation['de_pipeline'] = de_path.exists()
    
    # Check sample_sheet_dir
    if 'sample_sheet_dir' in config:
        sample_dir = Path(config['sample_sheet_dir'])
        validation['sample_sheet_dir'] = sample_dir.exists()
    
    return validation


def generate_bridge_config(
    project_summary_path: Path,
    de_pipeline_dir: Path = None,
    output_dir: Path = None,
    force: bool = False
) -> Path:
    """
    Automatically generate bridge configuration from project summary.
    
    Args:
        project_summary_path: Path to project_summary.json
        de_pipeline_dir: DE/GO pipeline directory (auto-detect if None)
        output_dir: Where to save config file (auto-detect if None)
        force: Overwrite existing config file
    
    Returns:
        Path to generated config file
    
    Raises:
        FileNotFoundError: If project summary doesn't exist
        ValueError: If required information is missing
    """
    if not project_summary_path.exists():
        raise FileNotFoundError(f"Project summary not found: {project_summary_path}")
    
    # Read project summary
    with open(project_summary_path) as f:
        summary = json.load(f)
    
    project_id = summary.get('project_id')
    if not project_id:
        raise ValueError("Project ID not found in summary")
    
    # Auto-detect paths
    rnaseq_output = project_summary_path.parent.resolve()
    
    # Find DE pipeline
    if de_pipeline_dir is None:
        de_candidates = [
            Path("/home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis"),
            Path("/data_3tb/shared/RNA-Seq_DE_GO_analysis"),
            Path.home() / "ngs-pipeline" / "RNA-Seq_DE_GO_analysis",
        ]
        for candidate in de_candidates:
            if candidate.exists():
                de_pipeline_dir = candidate
                break
        
        if de_pipeline_dir is None:
            raise ValueError("Could not auto-detect DE/GO pipeline directory")
    
    de_pipeline_dir = Path(de_pipeline_dir).resolve()
    
    # Discover sample sheet
    sample_sheet_path = discover_sample_sheet(project_id)
    sample_sheet_dir = str(sample_sheet_path.parent) if sample_sheet_path else None
    
    # Determine counts file location
    counts_relpath = 'project_summary/counts/counts_matrix_clean.csv'
    counts_file = rnaseq_output / counts_relpath
    
    # Fallback to other locations if not found
    if not counts_file.exists():
        alternative_paths = [
            'counts/counts_matrix_clean.csv',
            'counts/counts_matrix.txt',
            'final_outputs/counts/counts_matrix.txt',
        ]
        for alt_path in alternative_paths:
            if (rnaseq_output / alt_path).exists():
                counts_relpath = alt_path
                break
    
    # Build configuration
    config = {
        'project_id': project_id,
        'rnaseq_output': str(rnaseq_output),
        'de_pipeline': str(de_pipeline_dir),
        'counts_relpath': counts_relpath,
        'summary_relpath': 'project_summary.json'
    }
    
    if sample_sheet_dir:
        config['sample_sheet_dir'] = sample_sheet_dir
    
    # Validate configuration
    validation = validate_paths(config)
    missing_paths = [k for k, v in validation.items() if not v]
    
    if missing_paths:
        print(f"⚠️  Warning: Some paths not found: {', '.join(missing_paths)}")
    
    # Determine output location
    if output_dir is None:
        output_dir = Path("/data_3tb/shared/rna-seq-pipeline/config/projects")
        if not output_dir.exists():
            output_dir = Path.cwd() / "config" / "projects"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate config file path
    config_filename = f"paths_{project_id.replace('-', '_')}.yaml"
    config_path = output_dir / config_filename
    
    # Check if file exists
    if config_path.exists() and not force:
        print(f"ℹ️  Config already exists: {config_path}")
        print(f"   Use force=True to overwrite")
        return config_path
    
    # Write config file
    with open(config_path, 'w') as f:
        # Write header
        f.write("# Bridge Script Path Configuration\n")
        f.write("# Auto-generated configuration\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Source: {project_summary_path}\n")
        f.write(f"# Used by: scripts/bridge_to_de_pipeline.py\n")
        f.write(f"#\n")
        f.write(f"# Project: {project_id}\n")
        
        # Add validation status
        if missing_paths:
            f.write(f"# WARNING: Missing paths: {', '.join(missing_paths)}\n")
        
        f.write("\n")
        
        # Write configuration
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print(f"✅ Generated config: {config_path}")
    print(f"   RNA-seq output: {config['rnaseq_output']}")
    print(f"   DE pipeline: {config['de_pipeline']}")
    if sample_sheet_dir:
        print(f"   Sample sheet: {sample_sheet_dir}")
    else:
        print(f"   ⚠️  Sample sheet not found - will need manual metadata creation")
    
    return config_path


def ensure_bridge_config(
    project_id: str,
    rnaseq_output: Optional[str] = None,
    de_pipeline: Optional[str] = None,
    force: bool = False
) -> Dict[str, any]:
    """
    Ensure bridge config exists, create if missing.
    
    This is the high-level function for agent automation.
    
    Args:
        project_id: Project identifier
        rnaseq_output: Optional RNA-seq output directory (auto-detect if None)
        de_pipeline: Optional DE pipeline directory (auto-detect if None)
        force: Overwrite existing config
    
    Returns:
        Dictionary with:
            - status: 'exists' | 'created' | 'error'
            - config_path: Path to config file
            - message: Description
            - validation: Path validation results (if created)
    """
    # Check if config already exists
    config_candidates = [
        Path("/data_3tb/shared/rna-seq-pipeline/config/projects") / f"paths_{project_id.replace('-', '_')}.yaml",
        Path.cwd() / "config" / "projects" / f"paths_{project_id.replace('-', '_')}.yaml",
    ]
    
    existing_config = None
    for candidate in config_candidates:
        if candidate.exists():
            existing_config = candidate
            break
    
    if existing_config and not force:
        return {
            'status': 'exists',
            'config_path': str(existing_config),
            'message': f'Config already exists: {existing_config}'
        }
    
    # Auto-discover rnaseq_output if not provided
    if not rnaseq_output:
        output_base = Path("/data_3tb/shared/output")
        candidates = [
            output_base / project_id,
            output_base / project_id.replace('-', '_'),
            output_base / project_id.replace('_', '-')
        ]
        
        for candidate in candidates:
            summary_file = candidate / "project_summary.json"
            if summary_file.exists():
                rnaseq_output = str(candidate)
                break
    
    if not rnaseq_output:
        return {
            'status': 'error',
            'message': f'Could not find RNA-seq output for project {project_id}. Tried: {output_base}'
        }
    
    # Generate config
    try:
        project_summary = Path(rnaseq_output) / "project_summary.json"
        
        config_path = generate_bridge_config(
            project_summary,
            de_pipeline_dir=Path(de_pipeline) if de_pipeline else None,
            force=force
        )
        
        # Load and validate
        with open(config_path) as f:
            config = yaml.safe_load(f)
        
        validation = validate_paths(config)
        
        return {
            'status': 'created',
            'config_path': str(config_path),
            'message': f'Config created: {config_path}',
            'validation': validation,
            'config': config
        }
    
    except Exception as e:
        return {
            'status': 'error',
            'message': f'Failed to generate config: {str(e)}'
        }


if __name__ == '__main__':
    """
    Command-line usage:
    
    # Auto-generate config from project summary
    python auto_config.py /data_3tb/shared/output/mouse-chd8/project_summary.json
    
    # Or use project ID for auto-discovery
    python auto_config.py --project-id mouse-chd8
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Auto-generate bridge configuration from project summary"
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        'summary_file',
        nargs='?',
        type=Path,
        help='Path to project_summary.json'
    )
    group.add_argument(
        '--project-id',
        help='Project ID (will auto-discover paths)'
    )
    
    parser.add_argument(
        '--de-pipeline',
        type=Path,
        help='DE/GO pipeline directory (auto-detect if not provided)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory for config file'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing config file'
    )
    
    args = parser.parse_args()
    
    try:
        if args.project_id:
            # Use high-level function
            result = ensure_bridge_config(
                args.project_id,
                de_pipeline=str(args.de_pipeline) if args.de_pipeline else None,
                force=args.force
            )
            
            print(f"\n{'='*60}")
            print(f"Status: {result['status'].upper()}")
            print(f"{'='*60}")
            print(result['message'])
            
            if 'validation' in result:
                print(f"\nPath validation:")
                for path_name, exists in result['validation'].items():
                    status = "✅" if exists else "❌"
                    print(f"  {status} {path_name}")
        else:
            # Use low-level function
            config_path = generate_bridge_config(
                args.summary_file,
                de_pipeline_dir=args.de_pipeline,
                output_dir=args.output_dir,
                force=args.force
            )
            
            print(f"\n{'='*60}")
            print(f"SUCCESS")
            print(f"{'='*60}")
            print(f"Config file: {config_path}")
    
    except Exception as e:
        print(f"\n{'='*60}")
        print(f"ERROR")
        print(f"{'='*60}")
        print(f"{str(e)}")
        exit(1)
