#!/usr/bin/env python3
"""
Pipeline Execution Tools for LLM Agent

Tools to enable full pipeline execution from raw FASTQ files:
- create_project_config: Generate config.yaml from parameters
- detect_fastq_files: Scan directory and identify FASTQ files
- validate_input_data: Pre-flight checks before execution
- run_pipeline: Execute Snakemake workflow
"""

import json
import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional
import subprocess
import re
import os


def create_project_config(
    project_id: str,
    data_dir: str,
    results_dir: str,
    species: str = "human",
    read_type: str = "paired-end",
    use_sample_sheet: bool = False,
    template: str = "config/config.yaml"
) -> Dict[str, Any]:
    """
    Generate project-specific config.yaml from parameters.
    
    Args:
        project_id: Project identifier
        data_dir: Directory containing FASTQ files
        results_dir: Output directory for results
        species: Species (human, mouse, etc.)
        read_type: paired-end or single-end
        use_sample_sheet: Whether to use sample sheet
        template: Template config file to use as base
    
    Returns:
        {
            "status": "success"|"error",
            "config_path": str,
            "message": str
        }
    """
    try:
        # Load template config
        if Path(template).exists():
            with open(template) as f:
                config = yaml.safe_load(f)
        else:
            config = {}
        
        # Update with project-specific values
        config.update({
            'project_id': project_id,
            'use_standard_structure': True,
            'data_dir': str(Path(data_dir).resolve()),
            'base_results_dir': str(Path(results_dir).resolve()),
            'species': species,
            'use_sample_sheet': use_sample_sheet,
            'pipeline_type': 'rna-seq'
        })
        
        # Set species-specific reference paths
        if species == "human":
            config['genome_fasta'] = "genome/human/GRCh38.fa"
            config['genes_gtf'] = "genome/human/gencode.v44.annotation.gtf"
        elif species == "mouse":
            config['genome_fasta'] = "genome/mouse/GRCm39.fa"
            config['genes_gtf'] = "genome/mouse/gencode.vM33.annotation.gtf"
        
        # Save config
        config_dir = Path("config/projects")
        config_dir.mkdir(parents=True, exist_ok=True)
        
        config_path = config_dir / f"{project_id}.yaml"
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        
        return {
            "status": "success",
            "config_path": str(config_path),
            "message": f"Config created: {config_path}",
            "config": config
        }
    
    except Exception as e:
        return {
            "status": "error",
            "message": f"Failed to create config: {str(e)}"
        }


def detect_fastq_files(
    data_dir: str,
    pattern: str = "*{,_,-}{R,read,}{1,2}{,_001}.{fastq,fq}{,.gz}"
) -> Dict[str, Any]:
    """
    Scan directory and detect FASTQ files.
    
    Args:
        data_dir: Directory to scan
        pattern: Glob pattern for FASTQ files
    
    Returns:
        {
            "status": "success"|"error",
            "total_files": int,
            "samples": [{"sample_id": str, "R1": str, "R2": str, "size_gb": float}],
            "read_type": "paired-end"|"single-end",
            "total_size_gb": float
        }
    """
    try:
        data_path = Path(data_dir).resolve()
        
        if not data_path.exists():
            return {
                "status": "error",
                "message": f"Directory not found: {data_dir}"
            }
        
        # Find all FASTQ files
        fastq_files = []
        for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
            fastq_files.extend(data_path.glob(ext))
        
        if not fastq_files:
            return {
                "status": "error",
                "message": f"No FASTQ files found in {data_dir}"
            }
        
        # Group by sample
        samples = {}
        read_type = "single-end"
        
        for fq in fastq_files:
            # Try to extract sample name and read number
            name = fq.stem.replace('.fastq', '').replace('.fq', '').replace('.gz', '')
            
            # Common patterns for R1/R2
            if re.search(r'[_-]R?[12](?:_001)?$', name):
                # Paired-end
                read_type = "paired-end"
                sample_id = re.sub(r'[_-]R?[12](?:_001)?$', '', name)
                read_num = 'R1' if re.search(r'[_-]R?1(?:_001)?$', name) else 'R2'
            else:
                # Single-end
                sample_id = name
                read_num = 'R1'
            
            if sample_id not in samples:
                samples[sample_id] = {}
            
            samples[sample_id][read_num] = str(fq)
        
        # Calculate sizes and build result
        sample_list = []
        total_size = 0
        
        for sample_id, reads in sorted(samples.items()):
            sample_info = {"sample_id": sample_id}
            
            for read_key in ['R1', 'R2']:
                if read_key in reads:
                    fq_path = Path(reads[read_key])
                    size_gb = fq_path.stat().st_size / (1024**3)
                    sample_info[read_key] = reads[read_key]
                    total_size += size_gb
            
            sample_info['size_gb'] = sum(
                Path(r).stat().st_size / (1024**3) 
                for r in reads.values()
            )
            
            sample_list.append(sample_info)
        
        return {
            "status": "success",
            "total_files": len(fastq_files),
            "samples": sample_list,
            "total_samples": len(sample_list),
            "read_type": read_type,
            "total_size_gb": round(total_size, 2),
            "data_dir": str(data_path)
        }
    
    except Exception as e:
        return {
            "status": "error",
            "message": f"Failed to detect FASTQ files: {str(e)}"
        }


def validate_input_data(
    config_file: str,
    required_disk_gb: Optional[float] = None
) -> Dict[str, Any]:
    """
    Pre-flight validation before pipeline execution.
    
    Args:
        config_file: Path to config.yaml
        required_disk_gb: Required disk space (auto-estimated if None)
    
    Returns:
        {
            "status": "valid"|"invalid",
            "checks": {...},
            "warnings": [...],
            "errors": [...],
            "estimated_runtime": str
        }
    """
    try:
        config_path = Path(config_file)
        
        if not config_path.exists():
            return {
                "status": "invalid",
                "errors": [f"Config file not found: {config_file}"]
            }
        
        # Load config
        with open(config_path) as f:
            config = yaml.safe_load(f)
        
        checks = {}
        warnings = []
        errors = []
        
        # Check 1: Data directory exists
        data_dir = Path(config.get('data_dir', ''))
        if data_dir.exists():
            checks['data_dir'] = True
            
            # Count FASTQ files
            fastq_result = detect_fastq_files(str(data_dir))
            if fastq_result['status'] == 'success':
                checks['fastq_files'] = True
                num_samples = fastq_result['total_samples']
                total_size_gb = fastq_result['total_size_gb']
            else:
                checks['fastq_files'] = False
                errors.append(fastq_result['message'])
        else:
            checks['data_dir'] = False
            errors.append(f"Data directory not found: {data_dir}")
        
        # Check 2: Output directory writable
        results_dir = Path(config.get('base_results_dir', 'results'))
        try:
            results_dir.mkdir(parents=True, exist_ok=True)
            checks['output_dir'] = True
        except Exception as e:
            checks['output_dir'] = False
            errors.append(f"Cannot create output directory: {e}")
        
        # Check 3: Disk space
        if checks.get('fastq_files'):
            # Estimate required space (input + ~3x for intermediate + output)
            required_space = total_size_gb * 4
            
            stat = os.statvfs(results_dir)
            available_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)
            
            checks['disk_space'] = available_gb >= required_space
            
            if available_gb < required_space:
                errors.append(
                    f"Insufficient disk space: {required_space:.1f}GB required, "
                    f"{available_gb:.1f}GB available"
                )
            elif available_gb < required_space * 1.5:
                warnings.append(
                    f"Low disk space: {required_space:.1f}GB required, "
                    f"{available_gb:.1f}GB available (recommend 1.5x buffer)"
                )
        
        # Check 4: Reference genome exists (if specified)
        genome_fasta = config.get('genome_fasta', '')
        if genome_fasta:
            genome_path = Path(genome_fasta)
            if genome_path.exists():
                checks['reference_genome'] = True
            else:
                checks['reference_genome'] = False
                errors.append(f"Reference genome not found: {genome_fasta}")
        
        # Check 5: GTF annotation exists
        genes_gtf = config.get('genes_gtf', '')
        if genes_gtf:
            gtf_path = Path(genes_gtf)
            if gtf_path.exists():
                checks['annotation'] = True
            else:
                checks['annotation'] = False
                errors.append(f"GTF annotation not found: {genes_gtf}")
        
        # Estimate runtime
        if checks.get('fastq_files'):
            # Rough estimate: ~30 min per sample on 8 cores
            estimated_hours = (num_samples * 30) / 60 / 8 * config.get('cores', 8)
            estimated_runtime = f"{estimated_hours:.1f}-{estimated_hours*1.5:.1f} hours ({num_samples} samples, {config.get('cores', 8)} cores)"
        else:
            estimated_runtime = "Unknown"
        
        status = "valid" if not errors else "invalid"
        
        return {
            "status": status,
            "checks": checks,
            "warnings": warnings,
            "errors": errors,
            "estimated_runtime": estimated_runtime,
            "summary": {
                "samples": num_samples if checks.get('fastq_files') else 0,
                "total_size_gb": total_size_gb if checks.get('fastq_files') else 0,
                "available_disk_gb": available_gb if 'available_gb' in locals() else 0
            }
        }
    
    except Exception as e:
        return {
            "status": "error",
            "message": f"Validation failed: {str(e)}"
        }


def run_pipeline(
    config_file: str,
    cores: int = 8,
    dry_run: bool = True,
    until_rule: Optional[str] = None,
    force: bool = False
) -> Dict[str, Any]:
    """
    Execute Snakemake workflow.
    
    Args:
        config_file: Path to config.yaml
        cores: Number of cores to use
        dry_run: If True, show what would be run without executing
        until_rule: Stop at specific rule (e.g., 'fastqc', 'multiqc')
        force: Force re-run even if outputs exist
    
    Returns:
        {
            "status": "success"|"error"|"dry_run",
            "returncode": int,
            "jobs": int (for dry-run),
            "output": str,
            "command": str
        }
    """
    try:
        # Validate config exists
        config_path = Path(config_file)
        if not config_path.exists():
            return {
                "status": "error",
                "message": f"Config file not found: {config_file}"
            }
        
        # Build snakemake command
        cmd = [
            "snakemake",
            "--configfile", str(config_path),
            "--cores", str(cores)
        ]
        
        if dry_run:
            cmd.append("--dry-run")
        
        if until_rule:
            cmd.extend(["--until", until_rule])
        
        if force:
            cmd.append("--forceall")
        
        # Execute
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path.cwd()
        )
        
        # Parse output
        output = result.stdout + result.stderr
        
        # Extract job count from dry-run
        jobs = 0
        if dry_run:
            match = re.search(r'Job counts:.*?(\d+)', output, re.DOTALL)
            if match:
                # Sum all job counts
                job_counts = re.findall(r'\s+(\d+)', output[match.start():])
                jobs = sum(int(j) for j in job_counts)
        
        status = "dry_run" if dry_run else ("success" if result.returncode == 0 else "error")
        
        return {
            "status": status,
            "returncode": result.returncode,
            "jobs": jobs if dry_run else None,
            "output": output,
            "command": " ".join(cmd)
        }
    
    except Exception as e:
        return {
            "status": "error",
            "message": f"Pipeline execution failed: {str(e)}"
        }


# CLI for testing
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Pipeline execution tools")
    subparsers = parser.add_subparsers(dest='command')
    
    # create-config subcommand
    create_parser = subparsers.add_parser('create-config')
    create_parser.add_argument('--project-id', required=True)
    create_parser.add_argument('--data-dir', required=True)
    create_parser.add_argument('--results-dir', required=True)
    create_parser.add_argument('--species', default='human')
    
    # detect-fastq subcommand
    detect_parser = subparsers.add_parser('detect-fastq')
    detect_parser.add_argument('--data-dir', required=True)
    
    # validate subcommand
    validate_parser = subparsers.add_parser('validate')
    validate_parser.add_argument('--config', required=True)
    
    # run subcommand
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('--config', required=True)
    run_parser.add_argument('--cores', type=int, default=8)
    run_parser.add_argument('--dry-run', action='store_true')
    
    args = parser.parse_args()
    
    if args.command == 'create-config':
        result = create_project_config(
            args.project_id, args.data_dir, args.results_dir, args.species
        )
    elif args.command == 'detect-fastq':
        result = detect_fastq_files(args.data_dir)
    elif args.command == 'validate':
        result = validate_input_data(args.config)
    elif args.command == 'run':
        result = run_pipeline(args.config, args.cores, args.dry_run)
    else:
        parser.print_help()
        exit(1)
    
    print(json.dumps(result, indent=2, ensure_ascii=False))
