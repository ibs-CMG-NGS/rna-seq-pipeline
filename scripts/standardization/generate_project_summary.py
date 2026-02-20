#!/usr/bin/env python3
"""
Generate project-level summary from all sample manifests.
This enables agent to query project status at a glance.
"""

import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any
import sys

def collect_sample_manifests(project_dir: Path) -> List[Dict]:
    """Collect all manifest.json files in project directory."""
    manifests = []
    manifest_files = list(project_dir.rglob("*/rna-seq/final_outputs/manifest.json"))
    
    for mf in sorted(manifest_files):
        try:
            with open(mf) as f:
                manifest = json.load(f)
                manifests.append(manifest)
        except Exception as e:
            print(f"Warning: Failed to load {mf}: {e}", file=sys.stderr)
    
    return manifests

def summarize_qc_status(manifests: List[Dict]) -> Dict[str, Any]:
    """Summarize QC status across all samples."""
    total = len(manifests)
    passed = sum(1 for m in manifests if m.get('qc_metrics', {}).get('overall_status') == 'PASS')
    failed = total - passed
    
    statuses = {}
    for m in manifests:
        status = m.get('qc_metrics', {}).get('overall_status', 'UNKNOWN')
        sample_id = m.get('sample_id', 'unknown')
        statuses[sample_id] = status
    
    return {
        "total_samples": total,
        "passed": passed,
        "failed": failed,
        "pass_rate": round(passed / total * 100, 2) if total > 0 else 0,
        "sample_statuses": statuses
    }

def group_by_condition(manifests: List[Dict]) -> Dict[str, List[str]]:
    """Group samples by experimental condition."""
    groups = {}
    for m in manifests:
        condition = m.get('sample_metadata', {}).get('condition', 'Unknown')
        sample_id = m.get('sample_id', 'unknown')
        
        if condition not in groups:
            groups[condition] = []
        groups[condition].append(sample_id)
    
    return groups

def extract_qc_metrics(manifests: List[Dict]) -> Dict[str, Any]:
    """Extract key QC metrics for all samples."""
    metrics = {}
    
    for m in manifests:
        sample_id = m.get('sample_id', 'unknown')
        qc = m.get('qc_metrics', {})
        
        metrics[sample_id] = {
            "uniquely_mapped_pct": qc.get('alignment', {}).get('uniquely_mapped_pct', 0),
            "assignment_rate": qc.get('quantification', {}).get('assignment_rate', 0),
            "total_reads": qc.get('sequencing', {}).get('total_reads', 0),
            "status": qc.get('overall_status', 'UNKNOWN')
        }
    
    return metrics

def calculate_aggregate_stats(manifests: List[Dict]) -> Dict[str, float]:
    """Calculate aggregate statistics across all samples."""
    if not manifests:
        return {}
    
    unique_map_rates = []
    assign_rates = []
    
    for m in manifests:
        qc = m.get('qc_metrics', {})
        unique_map_rates.append(qc.get('alignment', {}).get('uniquely_mapped_pct', 0))
        assign_rates.append(qc.get('quantification', {}).get('assignment_rate', 0))
    
    return {
        "mean_uniquely_mapped_pct": round(sum(unique_map_rates) / len(unique_map_rates), 2),
        "mean_assignment_rate": round(sum(assign_rates) / len(assign_rates), 2),
        "min_uniquely_mapped_pct": round(min(unique_map_rates), 2),
        "max_uniquely_mapped_pct": round(max(unique_map_rates), 2),
        "min_assignment_rate": round(min(assign_rates), 2),
        "max_assignment_rate": round(max(assign_rates), 2)
    }

def identify_issues(manifests: List[Dict]) -> List[Dict[str, str]]:
    """Identify samples with issues."""
    issues = []
    
    for m in manifests:
        sample_id = m.get('sample_id', 'unknown')
        qc = m.get('qc_metrics', {})
        
        if qc.get('overall_status') == 'FAIL':
            sample_issues = qc.get('issues', [])
            for issue in sample_issues:
                issues.append({
                    "sample_id": sample_id,
                    "severity": issue.get('severity', 'UNKNOWN'),
                    "message": issue.get('message', '')
                })
    
    return issues

def generate_project_summary(project_dir: Path, project_id: str, output_file: Path) -> None:
    """Generate comprehensive project summary."""
    
    print(f"Collecting manifests from: {project_dir}")
    manifests = collect_sample_manifests(project_dir)
    
    if not manifests:
        print("ERROR: No manifest files found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(manifests)} sample manifests")
    
    # Build summary
    summary = {
        "project_id": project_id,
        "pipeline_type": manifests[0].get('pipeline_type', 'rna-seq') if manifests else 'rna-seq',
        "generated_date": datetime.now().isoformat(),
        
        "qc_summary": summarize_qc_status(manifests),
        "condition_groups": group_by_condition(manifests),
        "aggregate_stats": calculate_aggregate_stats(manifests),
        "sample_metrics": extract_qc_metrics(manifests),
        "issues": identify_issues(manifests),
        
        "files": {
            "project_dir": str(project_dir),
            "sample_manifests": [m.get('sample_id') for m in manifests]
        }
    }
    
    # Write summary
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✅ Project summary written to: {output_file}")
    print(f"\nSummary:")
    print(f"  Total samples: {summary['qc_summary']['total_samples']}")
    print(f"  Passed QC: {summary['qc_summary']['passed']}")
    print(f"  Failed QC: {summary['qc_summary']['failed']}")
    print(f"  Pass rate: {summary['qc_summary']['pass_rate']}%")
    
    if summary['condition_groups']:
        print(f"\nConditions:")
        for cond, samples in summary['condition_groups'].items():
            print(f"  {cond}: {len(samples)} samples")

def main():
    parser = argparse.ArgumentParser(
        description="Generate project-level summary from sample manifests"
    )
    parser.add_argument(
        '--project-dir',
        type=Path,
        required=True,
        help='Project directory containing sample subdirectories'
    )
    parser.add_argument(
        '--project-id',
        required=True,
        help='Project identifier'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output path for project_summary.json'
    )
    
    args = parser.parse_args()
    
    if not args.project_dir.exists():
        print(f"ERROR: Project directory not found: {args.project_dir}", file=sys.stderr)
        sys.exit(1)
    
    generate_project_summary(args.project_dir, args.project_id, args.output)

if __name__ == '__main__':
    main()
