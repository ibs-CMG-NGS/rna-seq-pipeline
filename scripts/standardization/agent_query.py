#!/usr/bin/env python3
"""
Agent Query Interface for Pipeline Results

Enables natural language queries about pipeline status:
- "What's the QC status?"
- "Show failed samples"
- "Compare conditions"
- "Which samples need reprocessing?"
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Any
import sys

class PipelineQueryAgent:
    """Agent interface for querying pipeline results."""
    
    def __init__(self, project_summary_path: Path):
        """Initialize with project summary."""
        with open(project_summary_path) as f:
            self.summary = json.load(f)
        
        self.project_dir = Path(self.summary['files']['project_dir'])
        self.project_id = self.summary['project_id']
    
    def get_overall_status(self) -> Dict[str, Any]:
        """Get overall project status."""
        qc = self.summary['qc_summary']
        return {
            "project_id": self.project_id,
            "total_samples": qc['total_samples'],
            "passed": qc['passed'],
            "failed": qc['failed'],
            "pass_rate": f"{qc['pass_rate']}%"
        }
    
    def get_failed_samples(self) -> List[str]:
        """Get list of failed samples."""
        failed = []
        for sample, status in self.summary['qc_summary']['sample_statuses'].items():
            if status != 'PASS':
                failed.append(sample)
        return failed
    
    def get_samples_by_condition(self, condition: str) -> List[str]:
        """Get samples for specific condition."""
        groups = self.summary['condition_groups']
        return groups.get(condition, [])
    
    def get_all_conditions(self) -> List[str]:
        """Get list of all conditions."""
        return list(self.summary['condition_groups'].keys())

    def get_sample_axes(self) -> Dict[str, Any]:
        """
        Get all experimental axes and their values.
        Returns genotype, tissue, sex groupings.
        """
        axes = self.summary.get('sample_axes', {})
        result = {}
        for axis, groups in axes.items():
            result[axis] = {
                "values": list(groups.keys()),
                "counts": {v: len(s) for v, s in groups.items()}
            }
        return result

    def filter_samples(self, filters: Dict[str, str]) -> List[str]:
        """
        Filter samples by multiple axes simultaneously.
        
        Args:
            filters: e.g. {"tissue": "HPC", "sex": "Male", "genotype": "wildtype"}
        
        Returns list of sample IDs matching ALL filters.
        """
        axes = self.summary.get('sample_axes', {})
        
        # Start with all samples, then intersect per filter
        matching = set(self.summary['qc_summary']['sample_statuses'].keys())
        
        for axis, value in filters.items():
            axis_groups = axes.get(axis, {})
            # Case-insensitive match
            matched_value = next(
                (v for v in axis_groups if v.lower() == value.lower()), None
            )
            if matched_value:
                matching &= set(axis_groups[matched_value])
            else:
                matching = set()  # No match → empty
                break
        
        return sorted(matching)

    def compare_by_axis(self, axis: str, filters: Dict[str, str] = None) -> Dict[str, Any]:
        """
        Compare QC metrics grouped by a specific axis, with optional pre-filtering.
        
        Args:
            axis:    "genotype" | "tissue" | "sex"
            filters: optional pre-filter e.g. {"tissue": "HPC"}
        
        Examples:
            compare_by_axis("genotype", {"tissue": "HPC"})
                → wildtype vs heterozygous, HPC samples only
            compare_by_axis("tissue", {"genotype": "wildtype"})
                → HPC vs PFC, wildtype only
        """
        axes = self.summary.get('sample_axes', {})
        axis_groups = axes.get(axis, {})
        
        if not axis_groups:
            return {"error": f"Unknown axis: {axis}. Available: {list(axes.keys())}"}
        
        comparison = {}
        for value, samples in axis_groups.items():
            # Apply optional pre-filter
            if filters:
                samples = [s for s in samples if s in self.filter_samples(filters)]
            
            if not samples:
                continue
            
            metrics = [
                self.summary['sample_metrics'][s]
                for s in samples
                if s in self.summary['sample_metrics']
            ]
            
            if not metrics:
                continue
            
            comparison[value] = {
                "n_samples": len(samples),
                "mean_uniquely_mapped": round(
                    sum(m['uniquely_mapped_pct'] for m in metrics) / len(metrics), 2
                ),
                "mean_assignment_rate": round(
                    sum(m['assignment_rate'] for m in metrics) / len(metrics), 2
                ),
                "passed": sum(1 for m in metrics if m['status'] == 'PASS'),
                "failed": sum(1 for m in metrics if m['status'] != 'PASS'),
                "samples": samples
            }
        
        return {
            "axis": axis,
            "filters_applied": filters or {},
            "comparison": comparison
        }
    
    def get_sample_details(self, sample_id: str) -> Dict[str, Any]:
        """Get detailed metrics for a sample."""
        metrics = self.summary['sample_metrics'].get(sample_id, {})
        
        # Load full manifest for more details
        manifest_path = self.project_dir / sample_id / "rna-seq" / "final_outputs" / "manifest.json"
        if manifest_path.exists():
            with open(manifest_path) as f:
                manifest = json.load(f)
            
            return {
                "sample_id": sample_id,
                "status": metrics.get('status', 'UNKNOWN'),
                "qc_metrics": metrics,
                "metadata": manifest.get('sample_metadata', {}),
                "issues": manifest.get('qc_metrics', {}).get('issues', []),
                "next_steps": manifest.get('next_steps', [])
            }
        
        return {"sample_id": sample_id, "metrics": metrics}
    
    def compare_conditions(self) -> Dict[str, Dict[str, float]]:
        """Compare QC metrics across conditions."""
        comparison = {}
        
        for condition, samples in self.summary['condition_groups'].items():
            metrics = [self.summary['sample_metrics'][s] for s in samples]
            
            comparison[condition] = {
                "n_samples": len(samples),
                "mean_uniquely_mapped": round(sum(m['uniquely_mapped_pct'] for m in metrics) / len(metrics), 2),
                "mean_assignment_rate": round(sum(m['assignment_rate'] for m in metrics) / len(metrics), 2),
                "passed": sum(1 for m in metrics if m['status'] == 'PASS'),
                "failed": sum(1 for m in metrics if m['status'] != 'PASS')
            }
        
        return comparison
    
    def get_all_issues(self) -> List[Dict[str, str]]:
        """Get all issues across project."""
        return self.summary['issues']
    
    def get_aggregate_stats(self) -> Dict[str, float]:
        """Get aggregate statistics."""
        return self.summary['aggregate_stats']
    
    def query(self, query_type: str, **kwargs) -> Any:
        """General query interface."""
        
        query_map = {
            'status': self.get_overall_status,
            'failed': self.get_failed_samples,
            'conditions': self.get_all_conditions,
            'compare': self.compare_conditions,
            'issues': self.get_all_issues,
            'stats': self.get_aggregate_stats,
            'axes': self.get_sample_axes,
        }
        
        if query_type == 'sample':
            return self.get_sample_details(kwargs.get('sample_id'))
        elif query_type == 'condition':
            return self.get_samples_by_condition(kwargs.get('condition'))
        elif query_type == 'filter':
            return self.filter_samples(kwargs.get('filters', {}))
        elif query_type == 'compare_axis':
            return self.compare_by_axis(
                axis=kwargs.get('axis', 'genotype'),
                filters=kwargs.get('filters')
            )
        elif query_type in query_map:
            return query_map[query_type]()
        else:
            raise ValueError(f"Unknown query type: {query_type}")

def main():
    parser = argparse.ArgumentParser(
        description="Query pipeline results for agent-based analysis"
    )
    parser.add_argument(
        '--project-summary',
        type=Path,
        required=True,
        help='Path to project_summary.json'
    )
    parser.add_argument(
        '--query',
        required=True,
        choices=['status', 'failed', 'conditions', 'compare', 'issues', 'stats',
                 'sample', 'condition', 'axes', 'filter', 'compare_axis'],
        help='Type of query to perform'
    )
    parser.add_argument(
        '--sample-id',
        help='Sample ID for sample query'
    )
    parser.add_argument(
        '--condition',
        help='Condition name for condition query'
    )
    parser.add_argument(
        '--axis',
        default='condition',
        help='Axis for compare_axis query (e.g. condition, tissue, sex). '
             'Any column from the samplesheet is valid.'
    )
    parser.add_argument(
        '--filter',
        action='append',
        metavar='AXIS=VALUE',
        help='Filter samples by axis=value (repeatable). e.g. --filter tissue=HPC --filter sex=Male'
    )
    parser.add_argument(
        '--format',
        choices=['json', 'text'],
        default='json',
        help='Output format'
    )
    
    args = parser.parse_args()
    
    if not args.project_summary.exists():
        print(f"ERROR: Project summary not found: {args.project_summary}", file=sys.stderr)
        sys.exit(1)
    
    # Initialize agent
    agent = PipelineQueryAgent(args.project_summary)
    
    # Execute query
    try:
        if args.query == 'sample':
            if not args.sample_id:
                print("ERROR: --sample-id required for sample query", file=sys.stderr)
                sys.exit(1)
            result = agent.query('sample', sample_id=args.sample_id)
        elif args.query == 'condition':
            if not args.condition:
                print("ERROR: --condition required for condition query", file=sys.stderr)
                sys.exit(1)
            result = agent.query('condition', condition=args.condition)
        elif args.query == 'filter':
            filters = {}
            for f in (args.filter or []):
                if '=' in f:
                    k, v = f.split('=', 1)
                    filters[k.strip()] = v.strip()
            result = agent.query('filter', filters=filters)
        elif args.query == 'compare_axis':
            filters = {}
            for f in (args.filter or []):
                if '=' in f:
                    k, v = f.split('=', 1)
                    filters[k.strip()] = v.strip()
            result = agent.query('compare_axis', axis=args.axis, filters=filters or None)
        else:
            result = agent.query(args.query)
        
        # Output
        if args.format == 'json':
            print(json.dumps(result, indent=2))
        else:
            # Text format
            if args.query == 'status':
                print(f"Project: {result['project_id']}")
                print(f"Total: {result['total_samples']} samples")
                print(f"Passed: {result['passed']} | Failed: {result['failed']}")
                print(f"Pass rate: {result['pass_rate']}")
            elif args.query == 'failed':
                if result:
                    print("Failed samples:")
                    for s in result:
                        print(f"  - {s}")
                else:
                    print("No failed samples")
            elif args.query == 'conditions':
                print("Available conditions:")
                for c in result:
                    print(f"  - {c}")
            else:
                print(json.dumps(result, indent=2))
    
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
