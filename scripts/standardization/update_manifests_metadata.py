#!/usr/bin/env python3
"""
Update existing manifest.json files with sample metadata from sample sheet.
Useful when samples were processed before sample sheet integration.
"""

import json
import argparse
from pathlib import Path
import pandas as pd
import sys

def load_sample_sheet(sample_sheet_path: Path) -> dict:
    """Load sample sheet and create sample_id -> metadata mapping."""
    df = pd.read_csv(sample_sheet_path, sep='\t')
    
    metadata_map = {}
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        
        # Extract relevant metadata fields
        metadata = {}
        metadata_fields = [
            'condition', 'replicate', 'tissue', 'age', 'sex',
            'sequencing_platform', 'library_type', 'read_type',
            'species', 'genome_build', 'batch', 'treatment', 'time_point', 'notes'
        ]
        
        for field in metadata_fields:
            if field in row and pd.notna(row[field]):
                # Convert to native Python type
                value = row[field]
                if isinstance(value, (int, float)):
                    if pd.isna(value):
                        continue
                    if isinstance(value, float) and value.is_integer():
                        value = int(value)
                metadata[field] = value
        
        metadata_map[sample_id] = metadata
    
    return metadata_map

def update_manifest(manifest_path: Path, sample_metadata: dict) -> bool:
    """Update a single manifest file with sample metadata."""
    try:
        # Load existing manifest
        with open(manifest_path) as f:
            manifest = json.load(f)
        
        sample_id = manifest.get('sample_id')
        if not sample_id:
            print(f"  ⚠️  No sample_id in manifest: {manifest_path}")
            return False
        
        # Add or update sample_metadata
        if sample_id in sample_metadata:
            manifest['sample_metadata'] = sample_metadata[sample_id]
            
            # Write updated manifest
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)
            
            print(f"  ✅ Updated: {sample_id}")
            return True
        else:
            print(f"  ⚠️  Sample not in sheet: {sample_id}")
            return False
            
    except Exception as e:
        print(f"  ❌ Error updating {manifest_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Retroactively add sample metadata to existing manifest files"
    )
    parser.add_argument(
        '--sample-sheet',
        type=Path,
        required=True,
        help='Path to sample sheet (TSV format)'
    )
    parser.add_argument(
        '--project-dir',
        type=Path,
        required=True,
        help='Project directory containing sample subdirectories'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be updated without modifying files'
    )
    
    args = parser.parse_args()
    
    if not args.sample_sheet.exists():
        print(f"ERROR: Sample sheet not found: {args.sample_sheet}", file=sys.stderr)
        sys.exit(1)
    
    if not args.project_dir.exists():
        print(f"ERROR: Project directory not found: {args.project_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading sample sheet: {args.sample_sheet}")
    sample_metadata = load_sample_sheet(args.sample_sheet)
    print(f"Found {len(sample_metadata)} samples in sheet\n")
    
    # Find all manifest files
    manifest_files = list(args.project_dir.rglob("*/rna-seq/final_outputs/manifest.json"))
    print(f"Found {len(manifest_files)} manifest files\n")
    
    if args.dry_run:
        print("DRY RUN MODE - No files will be modified\n")
    
    print("Updating manifests:")
    updated = 0
    for manifest_path in sorted(manifest_files):
        if args.dry_run:
            # Just check if update would happen
            with open(manifest_path) as f:
                manifest = json.load(f)
            sample_id = manifest.get('sample_id')
            if sample_id in sample_metadata:
                print(f"  Would update: {sample_id}")
                updated += 1
        else:
            if update_manifest(manifest_path, sample_metadata):
                updated += 1
    
    print(f"\n{'Would update' if args.dry_run else 'Updated'} {updated}/{len(manifest_files)} manifests")
    
    if args.dry_run:
        print("\nRun without --dry-run to apply changes")
    else:
        print("\n✅ Done! Re-run generate_project_summary.py to see updated conditions")

if __name__ == '__main__':
    main()
