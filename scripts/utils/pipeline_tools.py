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
        stdout = result.stdout
        stderr = result.stderr
        combined = stdout + stderr

        if dry_run:
            # Parse dry-run summary instead of dumping raw output
            jobs_by_rule = {}

            # Snakemake ≥7: "Job stats:" table
            # Format:
            #   Job stats:
            #   job                      count
            #   ---------------------  -------
            #   all                          1
            #   cutadapt                    38
            #   ...
            stats_match = re.search(r'Job stats:\s*\n', combined)
            if stats_match:
                block = combined[stats_match.end():]
                for line in block.splitlines():
                    # Stop at blank line or next section header
                    if not line.strip() or line.startswith('=') or line.startswith('Building'):
                        break
                    # Skip header / separator lines
                    if re.match(r'^[-\s]+$', line) or line.lstrip().startswith('job') or line.lstrip().startswith('rule'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            jobs_by_rule[parts[0]] = int(parts[-1])
                        except ValueError:
                            pass

            # Fallback: "Job counts:" table (Snakemake ≤6)
            if not jobs_by_rule:
                counts_match = re.search(r'Job counts:\s*\ncount\s+jobs\s*\n([\s\S]+?)(?=\n\S|\Z)', combined)
                if counts_match:
                    for line in counts_match.group(1).strip().splitlines():
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                jobs_by_rule[parts[1]] = int(parts[0])
                            except ValueError:
                                pass

            total_jobs = sum(jobs_by_rule.values())

            # In a dry-run, "missing output" lines are NORMAL (jobs haven't run yet).
            # Only collect genuine errors: exceptions, tracebacks, config errors.
            issues = [
                ln.strip() for ln in combined.splitlines()
                if any(kw in ln.lower() for kw in ('error:', 'exception', 'traceback', 'syntaxerror'))
                and 'missing' not in ln.lower()  # exclude normal "Missing files" dry-run output
            ][:20]

            status = "dry_run_ok" if result.returncode == 0 else "dry_run_error"

            result_dict = {
                "status": status,
                "returncode": result.returncode,
                "dry_run": True,
                "total_jobs": total_jobs,
                "jobs_by_rule": jobs_by_rule,
                "issues": issues,
                "note": (
                    "dry_run_ok: All jobs listed above will be executed when the pipeline runs. "
                    "Missing output files are EXPECTED at this stage (pipeline has not run yet)."
                    if status == "dry_run_ok" else
                    "dry_run_error: Snakemake encountered a configuration error. Check issues[]."
                ),
                "command": " ".join(cmd),
            }

            # If parsing found nothing, include a snippet for debugging
            if not jobs_by_rule:
                result_dict["output_snippet"] = combined[:2000]

            return result_dict
        else:
            # Real run: return compact summary + last 50 lines of output
            output_lines = combined.splitlines()
            status = "success" if result.returncode == 0 else "error"
            return {
                "status": status,
                "returncode": result.returncode,
                "dry_run": False,
                "output_tail": "\n".join(output_lines[-50:]),
                "command": " ".join(cmd),
            }
    
    except Exception as e:
        return {
            "status": "error",
            "message": f"Pipeline execution failed: {str(e)}"
        }



# ─────────────────────────────────────────────────────────────────────────────
# Phase 8B tools
# ─────────────────────────────────────────────────────────────────────────────

def monitor_pipeline(
    project_id: str,
    base_results_dir: str,
    config_file: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Check pipeline execution status by inspecting completed output files and logs.

    Strategy (no live Snakemake process needed):
    1. Load config to determine expected rules / sample count.
    2. Count output files that exist vs expected per rule.
    3. Classify status: not_started / running / completed / error.

    Returns:
        {
            "status": "not_started"|"running"|"completed"|"error",
            "progress_pct": float,
            "completed_rules": {rule: n_done},
            "total_rules": {rule: n_expected},
            "current_stage": str,
            "errors": [...]
        }
    """
    try:
        base_path = Path(base_results_dir)
        project_dir = base_path / project_id

        # Path overlap fallback: if user passed the full project path as base_results_dir
        # (e.g. base=/data_3tb/output/mouse-chd8/ and project_id=mouse-chd8),
        # project_dir becomes .../mouse-chd8/mouse-chd8 which doesn't exist.
        # In that case, check if base_path itself contains sample-like subdirs and use it directly.
        if not project_dir.exists() and base_path.exists():
            candidate_dirs = [
                p for p in base_path.iterdir()
                if p.is_dir()
                and not p.name.startswith('.')
                and p.name not in ('logs', 'project_summary')
            ]
            if candidate_dirs:
                project_dir = base_path

        # Try to get sample count from config
        n_samples = 0
        if config_file and Path(config_file).exists():
            with open(config_file) as f:
                cfg = yaml.safe_load(f)
            data_dir = cfg.get('data_dir', '')
            if data_dir:
                fq = detect_fastq_files(data_dir)
                if fq['status'] == 'success':
                    n_samples = fq['total_samples']

        if n_samples == 0:
            # Fallback: count sample dirs
            if project_dir.exists():
                n_samples = sum(
                    1 for p in project_dir.iterdir()
                    if p.is_dir() and not p.name.startswith('.')
                    and p.name not in ('logs', 'project_summary')
                )

        if n_samples == 0:
            return {"status": "not_started", "message": f"Project directory not found: {project_dir}"}

        # Rules and their expected output counts
        expected = {
            'cutadapt':          n_samples,   # trimmed fastq logs
            'fastqc_raw':        n_samples * 2,
            'star_align':        n_samples,
            'featurecounts_quant': 1,
            'multiqc':           1,
            'generate_manifest': n_samples,
        }

        completed: Dict[str, int] = {}

        # Count per-rule outputs via manifest files (most reliable indicator)
        manifests = list(project_dir.glob('*/rna-seq/final_outputs/manifest.json'))
        completed['generate_manifest'] = len(manifests)

        # Count BAM files as proxy for star_align
        bams = list(project_dir.glob('*/rna-seq/final_outputs/*.bam'))
        completed['star_align'] = len(bams)

        # Count trimmed log files as proxy for cutadapt
        logs_dir = project_dir / 'logs'
        cutadapt_logs = list(logs_dir.glob('cutadapt/*.log')) if logs_dir.exists() else []
        completed['cutadapt'] = len(cutadapt_logs)

        # Count fastqc zips
        fastqc_zips = list(project_dir.glob('*/rna-seq/final_outputs/qc/*_fastqc.zip'))
        # also check project_summary qc dir
        fastqc_zips += list((project_dir / 'project_summary' / 'qc').glob('*_fastqc.zip'))
        completed['fastqc_raw'] = len(set(str(f) for f in fastqc_zips))

        # MultiQC report
        mqc = list(project_dir.glob('project_summary/qc/multiqc_report.html'))
        completed['multiqc'] = len(mqc)

        # featureCounts
        counts = list(project_dir.glob('project_summary/counts/*.txt'))
        completed['featurecounts_quant'] = min(len(counts), 1)

        # Compute overall progress
        total_exp = sum(expected.values())
        total_done = sum(min(completed.get(r, 0), expected[r]) for r in expected)
        progress = round(total_done / total_exp * 100, 1) if total_exp > 0 else 0.0

        # Determine current stage (last rule with partial completion)
        pipeline_order = ['cutadapt', 'fastqc_raw', 'star_align',
                          'featurecounts_quant', 'generate_manifest', 'multiqc']
        current_stage = 'not_started'
        for rule in pipeline_order:
            if completed.get(rule, 0) > 0:
                current_stage = rule
        if completed.get('multiqc', 0) >= 1:
            current_stage = 'completed'

        # Check for error logs
        error_logs = []
        if logs_dir.exists():
            for log_file in list(logs_dir.rglob('*.log'))[:50]:
                try:
                    text = log_file.read_text(errors='ignore')
                    if 'error' in text.lower() or 'exception' in text.lower():
                        error_logs.append(str(log_file.relative_to(project_dir)))
                except Exception:
                    pass

        status = (
            'not_started' if progress == 0
            else 'completed'  if progress >= 99.0
            else 'error'      if error_logs
            else 'running'
        )

        return {
            "status": status,
            "progress_pct": progress,
            "current_stage": current_stage,
            "completed_rules": completed,
            "total_rules": expected,
            "n_samples": n_samples,
            "project_dir": str(project_dir),
            "error_logs": error_logs[:5],
        }

    except Exception as e:
        return {"status": "error", "message": f"monitor_pipeline failed: {e}"}


def create_sample_sheet(
    project_id: str,
    data_dir: str,
    conditions: Optional[Dict[str, List[str]]] = None,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate a sample metadata TSV from detected FASTQ files.

    If `conditions` is provided, each key is a condition name and each value
    is a list of sample_id substrings to match (case-insensitive prefix/substr).
    If not provided, tries to infer condition from sample name patterns.

    Returns:
        {
            "status": "success"|"error",
            "tsv_path": str,
            "n_samples": int,
            "conditions_assigned": {condition: count},
            "unassigned": [sample_ids...]
        }
    """
    try:
        # Detect FASTQ files
        fq_result = detect_fastq_files(data_dir)
        if fq_result['status'] != 'success':
            return {"status": "error", "message": fq_result['message']}

        samples = fq_result['samples']

        # Build condition lookup: sample_id -> condition
        cond_map: Dict[str, str] = {}

        if conditions:
            for cond_name, patterns in conditions.items():
                for sample in samples:
                    sid = sample['sample_id']
                    for pat in patterns:
                        if pat.lower() in sid.lower():
                            cond_map[sid] = cond_name
                            break
        else:
            # Auto-infer: look for common suffixes like _WT _KO _Ctrl _Treat etc.
            auto_patterns = {
                'wildtype':      ['_W', '_WT', '_wt', '_wildtype', '_Ctrl', '_ctrl', '_control'],
                'heterozygous':  ['_H', '_Het', '_het', '_heterozygous', '_KO', '_ko'],
                'treatment':     ['_T', '_Treat', '_treat', '_treatment'],
            }
            for sample in samples:
                sid = sample['sample_id']
                for cond_name, pats in auto_patterns.items():
                    matched = False
                    for p in pats:
                        if sid.endswith(p) or f'{p}_' in sid:
                            cond_map[sid] = cond_name
                            matched = True
                            break
                    if matched:
                        break

        # Build TSV rows
        tsv_path = output_path or f"config/samples/{project_id}.tsv"
        Path(tsv_path).parent.mkdir(parents=True, exist_ok=True)

        rows = []
        unassigned = []
        for sample in sorted(samples, key=lambda s: s['sample_id']):
            sid = sample['sample_id']
            cond = cond_map.get(sid, 'unknown')
            if cond == 'unknown':
                unassigned.append(sid)
            rows.append({
                'sample_id': sid,
                'condition': cond,
                'replicate': '',
                'fastq_r1': sample.get('R1', ''),
                'fastq_r2': sample.get('R2', ''),
            })

        with open(tsv_path, 'w') as f:
            f.write('sample_id\tcondition\treplicate\tfastq_r1\tfastq_r2\n')
            for row in rows:
                f.write('\t'.join([
                    row['sample_id'], row['condition'], row['replicate'],
                    row['fastq_r1'], row['fastq_r2'],
                ]) + '\n')

        conditions_assigned: Dict[str, int] = {}
        for cond in cond_map.values():
            conditions_assigned[cond] = conditions_assigned.get(cond, 0) + 1

        return {
            "status": "success",
            "tsv_path": tsv_path,
            "n_samples": len(rows),
            "conditions_assigned": conditions_assigned,
            "unassigned": unassigned,
            "message": (
                f"Sample sheet created: {tsv_path}. "
                + (f"{len(unassigned)} samples unassigned — specify conditions manually." if unassigned else "All samples assigned.")
            ),
        }

    except Exception as e:
        return {"status": "error", "message": f"create_sample_sheet failed: {e}"}


def estimate_resources(
    config_file: str,
    cores: int = 8,
) -> Dict[str, Any]:
    """
    Estimate runtime, disk, and memory requirements before pipeline execution.

    Returns:
        {
            "status": "success"|"error",
            "n_samples": int,
            "total_input_gb": float,
            "estimated_output_gb": float,
            "estimated_hours_min": float,
            "estimated_hours_max": float,
            "recommended_cores": int,
            "available_disk_gb": float,
            "available_ram_gb": float,
            "disk_ok": bool,
            "ram_ok": bool,
        }
    """
    try:
        config_path = Path(config_file)
        if not config_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        data_dir = cfg.get('data_dir', '')
        results_dir = Path(cfg.get('base_results_dir', 'results'))

        # FASTQ scan
        fq = detect_fastq_files(data_dir)
        if fq['status'] != 'success':
            return {"status": "error", "message": fq['message']}

        n_samples = fq['total_samples']
        total_input_gb = fq['total_size_gb']

        # Output estimate: ~5x input (trimmed + aligned BAM + counts + QC)
        estimated_output_gb = round(total_input_gb * 5, 1)

        # Runtime estimate per sample on 1 core: ~45 min
        # With `cores` available, samples run in parallel (up to cores/2 at once)
        parallel = max(1, cores // 2)
        batches = (n_samples + parallel - 1) // parallel
        hours_min = round(batches * 0.5, 1)   # optimistic
        hours_max = round(batches * 1.0, 1)   # conservative

        # Disk availability
        try:
            results_dir.mkdir(parents=True, exist_ok=True)
            st = os.statvfs(str(results_dir))
            available_disk_gb = round(st.f_bavail * st.f_frsize / 1024**3, 1)
        except Exception:
            available_disk_gb = -1.0

        # RAM availability
        try:
            available_ram_gb = -1.0
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemAvailable'):
                        available_ram_gb = round(int(line.split()[1]) / 1024**2, 1)
                        break
        except Exception:
            available_ram_gb = -1.0

        # STAR needs ~30 GB RAM for mouse/human genome index
        star_ram_gb = 30.0
        ram_ok = available_ram_gb < 0 or available_ram_gb >= star_ram_gb
        disk_ok = available_disk_gb < 0 or available_disk_gb >= estimated_output_gb

        recommended_cores = min(cores, n_samples * 2)

        return {
            "status": "success",
            "n_samples": n_samples,
            "total_input_gb": total_input_gb,
            "estimated_output_gb": estimated_output_gb,
            "estimated_hours_min": hours_min,
            "estimated_hours_max": hours_max,
            "recommended_cores": recommended_cores,
            "available_disk_gb": available_disk_gb,
            "available_ram_gb": available_ram_gb,
            "disk_ok": disk_ok,
            "ram_ok": ram_ok,
            "warnings": (
                ([f"Low RAM: STAR needs ~{star_ram_gb}GB, only {available_ram_gb}GB available"]
                 if not ram_ok else [])
                + ([f"Low disk: need ~{estimated_output_gb}GB, only {available_disk_gb}GB available"]
                   if not disk_ok else [])
            ),
        }

    except Exception as e:
        return {"status": "error", "message": f"estimate_resources failed: {e}"}


# ─────────────────────────────────────────────────────────────────────────────
# CLI for testing
# ─────────────────────────────────────────────────────────────────────────────
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
