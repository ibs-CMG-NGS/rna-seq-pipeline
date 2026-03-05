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
    sample_sheet: str = None,
    genome_dir: str = None,
    genome_fasta: str = None,
    annotation_gtf: str = None,
    star_index: str = None,
    genome_build: str = None,
    threads: int = 12,
    memory_gb: int = 48,
    strandedness: int = 0,
    template: str = "config/config.yaml"
) -> Dict[str, Any]:
    """
    Generate a complete project-specific config.yaml.

    Produces a fully self-contained config that overrides ALL keys from
    config.yaml so no '/path/to/...' placeholder leaks through.

    Args:
        project_id: Project identifier (used as output subdirectory name)
        data_dir: Directory containing FASTQ files
        results_dir: Base output directory (project subdir created inside)
        species: 'human' | 'mouse' | 'rat' | other
        read_type: 'paired-end' | 'single-end'
        use_sample_sheet: Whether to use a sample sheet TSV
        sample_sheet: Path to sample sheet TSV (required when use_sample_sheet=True)
        genome_dir: Path to genome reference directory
        genome_fasta: Path to genome FASTA file
        annotation_gtf: Path to gene annotation GTF file
        star_index: Path to STAR genome index directory
        genome_build: Genome build string (e.g. 'GRCh38', 'GRCm38')
        threads: CPU threads for alignment/quantification
        memory_gb: RAM limit in GB
        strandedness: 0=unstranded, 1=forward, 2=reverse
        template: Unused (kept for backwards compat)

    Returns:
        {"status": "success"|"error", "config_path": str, "message": str, "config": dict}
    """
    try:
        # ── Species defaults ──────────────────────────────────────────────
        _species_defaults = {
            "human": {
                "genome_build": "GRCh38",
                "gc_min": 35, "gc_max": 65,
                "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            },
            "mouse": {
                "genome_build": "GRCm38",
                "gc_min": 40, "gc_max": 55,
                "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            },
            "rat": {
                "genome_build": "Rnor_6.0",
                "gc_min": 40, "gc_max": 55,
                "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            },
        }
        sp = _species_defaults.get(species.lower(), _species_defaults["human"])

        base_results = str(Path(results_dir).resolve())
        project_results = str(Path(results_dir).resolve() / project_id)

        config = {
            # ── Project ───────────────────────────────────────────────────
            'project_id': project_id,
            'project_name': project_id.replace('-', ' ').replace('_', ' ').title(),
            'pipeline_type': 'rna-seq',

            # ── Input ─────────────────────────────────────────────────────
            'data_dir': str(Path(data_dir).resolve()),
            'raw_data_subdir': '',

            # ── Sample sheet ──────────────────────────────────────────────
            'use_sample_sheet': use_sample_sheet,
            'sample_sheet': sample_sheet or f'config/samples/{project_id}.tsv',

            # ── Output structure ──────────────────────────────────────────
            'use_standard_structure': True,
            'base_results_dir': base_results,
            'results_dir': project_results,
            'logs_dir': f'{project_results}/logs',

            # ── Reference genome ──────────────────────────────────────────
            'genome_dir':     genome_dir     or '',
            'genome_fasta':   genome_fasta   or '',
            'annotation_gtf': annotation_gtf or '',
            'star_index':     star_index     or '',

            # ── Species ───────────────────────────────────────────────────
            'species': species,
            'genome_build': genome_build or sp['genome_build'],

            # ── Adapters ──────────────────────────────────────────────────
            'adapter_r1': sp['adapter_r1'],
            'adapter_r2': sp['adapter_r2'],

            # ── QC parameters ─────────────────────────────────────────────
            'quality_cutoff': 20,
            'min_read_length': 20,
            'fastqc_threads': 2,

            # ── Trimming ──────────────────────────────────────────────────
            'cutadapt_threads': 4,

            # ── Alignment ─────────────────────────────────────────────────
            'star_threads': threads,
            'star_memory_gb': min(35, memory_gb - 13),
            'star_sort_memory_bytes': 30000000000,
            'gzipped_fastq': True,
            'star_params': (
                '--outFilterMultimapNmax 20 --alignSJoverhangMin 8 '
                '--alignSJDBoverhangMin 1'
            ),

            # ── Quantification ────────────────────────────────────────────
            'featurecounts_threads': min(8, threads),
            'feature_type': 'exon',
            'attribute_type': 'gene_id',
            'strandedness': strandedness,
            'featurecounts_params': '-p -B -C',

            # ── QC report ─────────────────────────────────────────────────
            'generate_qc_report': True,
            'generate_multiqc': True,
            'qc_report_filename': 'qc_report.html',
            'qc_top_genes': 10,

            # ── FastQC auto-evaluation ────────────────────────────────────
            'fastqc_evaluation': {
                'enabled': True,
                'min_total_sequences': 5000000,
                'min_gc_content':  sp['gc_min'],
                'max_gc_content':  sp['gc_max'],
                'min_median_quality': 28,
                'min_lower_quartile': 20,
                'min_q30_percentage': 75,
                'max_n_content': 5.0,
                'critical_n_content': 10.0,
                'max_adapter_trimmed': 1.0,
                'warn_adapter_raw': 10.0,
                'evaluation_report': 'fastqc_evaluation.txt',
                'evaluation_json':   'fastqc_evaluation.json',
            },

            # ── QC thresholds (STAR alignment) ────────────────────────────
            'qc_thresholds': {
                'min_uniquely_mapped_pct': 70.0,
                'min_assignment_rate':     60.0,
                'max_mismatch_rate':        2.0,
            },

            # ── Resources ─────────────────────────────────────────────────
            'threads': threads,
            'memory_gb': memory_gb,
        }

        # ── Write config ──────────────────────────────────────────────────
        config_dir = Path("config/projects")
        config_dir.mkdir(parents=True, exist_ok=True)
        config_path = config_dir / f"{project_id}.yaml"

        with open(config_path, 'w') as f:
            f.write(f"# {config['project_name']} — generated by create_project_config\n")
            f.write(f"# Species: {species}  Genome: {config['genome_build']}\n")
            f.write(f"# ⚠️  Review genome paths before running!\n\n")
            yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

        return {
            "status": "success",
            "config_path": str(config_path),
            "message": (
                f"Config created: {config_path}\n"
                f"  → Review genome paths (genome_fasta, annotation_gtf, star_index)\n"
                f"  → strandedness={strandedness} — verify with RSeQC if unsure"
            ),
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
        
        # Check 1: Data directory / sample sheet FASTQ files
        data_dir_str = config.get('data_dir', '')
        sample_sheet_path = config.get('sample_sheet', '')

        # If data_dir not specified but sample_sheet is, derive data_dir from it
        if not data_dir_str and sample_sheet_path:
            ss_path = Path(sample_sheet_path)
            if not ss_path.is_absolute():
                # resolve relative to config file location
                ss_path = config_path.parent.parent / ss_path  # config/projects/ -> repo root
                if not ss_path.exists():
                    ss_path = config_path.parent / ss_path.name
            if ss_path.exists():
                import csv
                with open(ss_path) as sf:
                    reader = csv.DictReader(sf, delimiter='\t')
                    rows = list(reader)
                if rows and 'fastq_r1' in rows[0]:
                    first_fastq = Path(rows[0]['fastq_r1'])
                    data_dir_str = str(first_fastq.parent)

        data_dir = Path(data_dir_str) if data_dir_str else Path('')

        if data_dir_str and data_dir.exists():
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
        elif data_dir_str:
            checks['data_dir'] = False
            errors.append(f"Data directory not found: {data_dir}")
        else:
            checks['data_dir'] = False
            errors.append("No data_dir or sample_sheet specified in config")
        
        # Check 2: Output directory writable
        # Support both 'results_dir' and legacy 'base_results_dir' keys
        results_dir = Path(config.get('results_dir', config.get('base_results_dir', 'results')))
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
        
        # Check 5: GTF annotation exists (support both 'genes_gtf' and 'annotation_gtf' keys)
        genes_gtf = config.get('genes_gtf', '') or config.get('annotation_gtf', '')
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
        
        # Resolve pipeline root: walk up from config_file until Snakefile is found
        # This ensures snakemake is launched from the repo root regardless of CWD
        pipeline_root = config_path.resolve().parent
        for _ in range(5):
            if (pipeline_root / "Snakefile").exists():
                break
            pipeline_root = pipeline_root.parent
        else:
            # Fallback: use current working directory
            pipeline_root = Path.cwd()

        # Use absolute path for config so it works from any cwd
        cmd = [
            "snakemake",
            "--configfile", str(config_path.resolve()),
            "--cores", str(cores)
        ]
        
        if dry_run:
            cmd.append("--dry-run")
        
        if until_rule:
            cmd.extend(["--until", until_rule])
        
        if force:
            cmd.append("--forceall")
        
        # Execute from pipeline root (where Snakefile lives)
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(pipeline_root)
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
                "pipeline_root": str(pipeline_root),
                "note": (
                    "dry_run_ok: All jobs listed above will be executed when the pipeline runs. "
                    "Missing output files are EXPECTED at this stage (pipeline has not run yet)."
                    if status == "dry_run_ok" else
                    "dry_run_error: Snakemake encountered a configuration error. Check issues[] and stderr_snippet."
                ),
                "command": " ".join(cmd),
            }

            # Always include stderr snippet for debugging (first 3000 chars)
            if result.returncode != 0 or not jobs_by_rule:
                result_dict["stderr_snippet"] = (stderr or stdout)[:3000]

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


# ─────────────────────────────────────────────────────────────────────────────
# Phase 8C: Project/Pipeline information reading tools
# ─────────────────────────────────────────────────────────────────────────────

def read_project_config(config_file: str) -> Dict[str, Any]:
    """
    Read and summarise a project config YAML for the agent.

    Returns a human-readable summary of the project setup so the LLM can
    answer questions like "이 프로젝트 설정이 어떻게 돼 있어?" without the user
    having to provide every detail.
    """
    try:
        path = Path(config_file)
        if not path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(path) as f:
            cfg = yaml.safe_load(f)

        # Count samples via sample sheet if available
        n_samples = "unknown"
        conditions: List[str] = []
        sample_sheet_path = cfg.get("sample_sheet", "")
        if sample_sheet_path:
            ss = Path(sample_sheet_path)
            if not ss.is_absolute():
                ss = path.parent.parent / sample_sheet_path
            if ss.exists():
                import csv
                with open(ss) as f2:
                    reader = csv.DictReader(f2, delimiter="\t")
                    rows = list(reader)
                n_samples = len(rows)
                if rows and "condition" in rows[0]:
                    conditions = sorted(set(r["condition"] for r in rows))

        return {
            "status": "success",
            "project_id":      cfg.get("project_id"),
            "species":         cfg.get("species"),
            "genome_build":    cfg.get("genome_build"),
            "n_samples":       n_samples,
            "conditions":      conditions,
            "data_dir":        cfg.get("data_dir"),
            "results_dir":     cfg.get("base_results_dir") or cfg.get("results_dir"),
            "star_index":      cfg.get("star_index"),
            "annotation_gtf":  cfg.get("annotation_gtf") or cfg.get("genes_gtf"),
            "strandedness":    cfg.get("strandedness"),
            "threads":         cfg.get("threads") or cfg.get("star_threads"),
            "sample_sheet":    str(ss) if sample_sheet_path else None,
            "use_sample_sheet": cfg.get("use_sample_sheet"),
            "raw_config":      cfg,
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}


def read_sample_sheet(config_file: str) -> Dict[str, Any]:
    """
    Read the sample sheet referenced in a config file and return a structured
    summary: sample list, conditions, replicates, FASTQ paths.

    Answers: "샘플 목록 보여줘", "어떤 조건이 있어?", "wildtype 샘플이 몇 개야?"
    """
    try:
        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(cfg_path) as f:
            cfg = yaml.safe_load(f)

        sheet_path = cfg.get("sample_sheet", "")
        if not sheet_path:
            return {"status": "error", "message": "No sample_sheet key in config"}

        ss = Path(sheet_path)
        if not ss.is_absolute():
            ss = cfg_path.parent.parent / sheet_path
        if not ss.exists():
            return {"status": "error", "message": f"Sample sheet not found: {ss}"}

        import csv
        with open(ss) as f2:
            reader = csv.DictReader(f2, delimiter="\t")
            rows = list(reader)

        if not rows:
            return {"status": "error", "message": "Sample sheet is empty"}

        # Build summary
        conditions: Dict[str, List[str]] = {}
        tissues: Dict[str, List[str]] = {}
        for row in rows:
            sid = row.get("sample_id", "")
            cond = row.get("condition", "unknown")
            tissue = row.get("tissue", "")
            conditions.setdefault(cond, []).append(sid)
            if tissue:
                tissues.setdefault(tissue, []).append(sid)

        return {
            "status": "success",
            "n_samples": len(rows),
            "columns": list(rows[0].keys()),
            "conditions": {k: {"count": len(v), "samples": v} for k, v in conditions.items()},
            "tissues": {k: {"count": len(v), "samples": v} for k, v in tissues.items()} if tissues else {},
            "samples": rows,  # full list for detailed queries
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}


def read_qc_results(config_file: str) -> Dict[str, Any]:
    """
    Read QC results from the output directory:
    - FastQC evaluation JSON (per-sample pass/fail + metrics)
    - MultiQC general_stats (mapping rates, duplication, GC)
    - STAR log summary (uniquely mapped %)

    Answers: "QC 결과 어때?", "매핑률 낮은 샘플 있어?", "어떤 샘플이 실패했어?"
    """
    try:
        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(cfg_path) as f:
            cfg = yaml.safe_load(f)

        base = Path(cfg.get("base_results_dir") or cfg.get("results_dir", "results"))
        project_id = cfg.get("project_id", "")
        project_dir = base / project_id if (base / project_id).exists() else base
        qc_dir = project_dir / "project_summary" / "qc"

        result: Dict[str, Any] = {"status": "success", "project_dir": str(project_dir)}

        # 1. FastQC evaluation JSON
        eval_json = qc_dir / "fastqc_evaluation.json"
        if eval_json.exists():
            with open(eval_json) as f:
                eval_data = json.load(f)
            result["fastqc_evaluation"] = eval_data
        else:
            result["fastqc_evaluation"] = None

        # 2. MultiQC general stats TSV
        mqc_stats = qc_dir / "multiqc_data" / "multiqc_general_stats.txt"
        if not mqc_stats.exists():
            mqc_stats = next(qc_dir.glob("**/multiqc_general_stats.txt"), None)  # type: ignore

        if mqc_stats and Path(mqc_stats).exists():
            import csv
            with open(mqc_stats) as f:
                reader = csv.DictReader(f, delimiter="\t")
                mqc_rows = list(reader)
            result["multiqc_stats"] = mqc_rows
            result["multiqc_n_samples"] = len(mqc_rows)
        else:
            result["multiqc_stats"] = None

        # 3. STAR logs — extract uniquely mapped %
        star_logs = list(project_dir.glob("*/rna-seq/final_outputs/*Log.final.out"))
        star_logs += list(project_dir.glob("logs/star/*.log"))
        star_summary: List[Dict] = []
        for log in star_logs[:38]:
            try:
                text = Path(log).read_text()
                uniq_m = re.search(r'Uniquely mapped reads %\s+\|\s+([\d.]+)%', text)
                multi_m = re.search(r'% of reads mapped to multiple loci\s+\|\s+([\d.]+)%', text)
                unmapped_m = re.search(r'% of reads unmapped: too short\s+\|\s+([\d.]+)%', text)
                sample_name = log.parent.parent.parent.name if 'rna-seq' in str(log) else log.stem
                star_summary.append({
                    "sample": sample_name,
                    "uniquely_mapped_pct": float(uniq_m.group(1)) if uniq_m else None,
                    "multi_mapped_pct": float(multi_m.group(1)) if multi_m else None,
                    "unmapped_too_short_pct": float(unmapped_m.group(1)) if unmapped_m else None,
                })
            except Exception:
                pass

        if star_summary:
            avg_uniq = sum(s["uniquely_mapped_pct"] or 0 for s in star_summary) / len(star_summary)
            low_mapping = [s for s in star_summary if (s["uniquely_mapped_pct"] or 100) < 70]
            result["star_alignment"] = {
                "n_samples": len(star_summary),
                "avg_uniquely_mapped_pct": round(avg_uniq, 2),
                "low_mapping_samples": low_mapping,
                "per_sample": star_summary,
            }
        else:
            result["star_alignment"] = None

        # Overall summary
        n_pass = n_fail = 0
        if result["fastqc_evaluation"] and isinstance(result["fastqc_evaluation"], dict):
            for v in result["fastqc_evaluation"].values():
                if isinstance(v, dict):
                    if v.get("overall_pass"):
                        n_pass += 1
                    else:
                        n_fail += 1

        result["summary"] = {
            "fastqc_pass": n_pass,
            "fastqc_fail": n_fail,
            "star_avg_mapping_pct": result["star_alignment"]["avg_uniquely_mapped_pct"] if result["star_alignment"] else None,
            "multiqc_available": result["multiqc_stats"] is not None,
        }

        return result
    except Exception as e:
        return {"status": "error", "message": str(e)}


def read_counts(
    config_file: str,
    genes: Optional[List[str]] = None,
    top_n: int = 20,
) -> Dict[str, Any]:
    """
    Read the featureCounts expression matrix.

    Args:
        config_file: Path to project config
        genes: Specific gene IDs or names to look up (optional)
        top_n: Return top N most variable genes if genes not specified

    Answers: "CHD8 발현량 보여줘", "가장 많이 발현된 유전자 20개", "counts matrix 어디에 있어?"
    """
    try:
        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(cfg_path) as f:
            cfg = yaml.safe_load(f)

        base = Path(cfg.get("base_results_dir") or cfg.get("results_dir", "results"))
        project_id = cfg.get("project_id", "")
        project_dir = base / project_id if (base / project_id).exists() else base
        counts_dir = project_dir / "project_summary" / "counts"

        # Find counts matrix
        candidates = (
            list(counts_dir.glob("*.txt"))
            + list(counts_dir.glob("*.tsv"))
            + list(counts_dir.glob("*.csv"))
        )
        if not candidates:
            return {
                "status": "not_found",
                "message": f"No counts files found in {counts_dir}",
                "counts_dir": str(counts_dir),
            }

        counts_file = candidates[0]
        import csv

        # Detect delimiter
        sep = "\t" if counts_file.suffix in (".txt", ".tsv") else ","

        with open(counts_file) as f:
            reader = csv.DictReader(f, delimiter=sep)
            rows = list(reader)

        if not rows:
            return {"status": "error", "message": "Counts file is empty"}

        cols = list(rows[0].keys())
        # First column is typically gene_id; rest are sample counts
        gene_col = cols[0]
        sample_cols = cols[1:]

        # Filter to requested genes
        if genes:
            matched = [r for r in rows if any(
                g.lower() in r[gene_col].lower() for g in genes
            )]
            return {
                "status": "success",
                "counts_file": str(counts_file),
                "n_genes_total": len(rows),
                "genes_queried": genes,
                "results": matched,
                "samples": sample_cols,
            }

        # Return top N most variable genes
        def variance(row: Dict) -> float:
            vals = []
            for s in sample_cols:
                try:
                    vals.append(float(row[s]))
                except (ValueError, KeyError):
                    pass
            if len(vals) < 2:
                return 0.0
            mean = sum(vals) / len(vals)
            return sum((v - mean) ** 2 for v in vals) / len(vals)

        top_rows = sorted(rows, key=variance, reverse=True)[:top_n]

        return {
            "status": "success",
            "counts_file": str(counts_file),
            "n_genes_total": len(rows),
            "n_samples": len(sample_cols),
            "samples": sample_cols,
            "top_variable_genes": top_rows,
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}


def read_pipeline_logs(
    config_file: str,
    rule: Optional[str] = None,
    sample_id: Optional[str] = None,
    tail_lines: int = 50,
) -> Dict[str, Any]:
    """
    Read pipeline execution logs.

    Args:
        config_file: Path to project config
        rule: Specific rule to read logs for (cutadapt/star/fastqc/featurecounts)
        sample_id: Specific sample to read logs for
        tail_lines: Number of tail lines to return per log

    Answers: "로그 보여줘", "Chd8_HPC_1F_W 샘플 STAR 로그 어때?", "에러 있어?"
    """
    try:
        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(cfg_path) as f:
            cfg = yaml.safe_load(f)

        base = Path(cfg.get("base_results_dir") or cfg.get("results_dir", "results"))
        project_id = cfg.get("project_id", "")
        project_dir = base / project_id if (base / project_id).exists() else base
        logs_dir = project_dir / "logs"

        if not logs_dir.exists():
            # Also try config-specified logs_dir
            logs_dir = Path(cfg.get("logs_dir", str(project_dir / "logs")))

        if not logs_dir.exists():
            return {"status": "not_found", "message": f"Logs directory not found: {logs_dir}"}

        # Find matching log files
        if rule and sample_id:
            pattern = f"{rule}/{sample_id}*.log"
        elif rule:
            pattern = f"{rule}/*.log"
        elif sample_id:
            pattern = f"**/{sample_id}*.log"
        else:
            pattern = "**/*.log"

        log_files = list(logs_dir.glob(pattern))[:20]

        if not log_files:
            return {
                "status": "not_found",
                "message": f"No logs found matching rule={rule}, sample={sample_id}",
                "logs_dir": str(logs_dir),
                "available_rules": [p.name for p in logs_dir.iterdir() if p.is_dir()],
            }

        logs_out = []
        errors_found = []
        for lf in sorted(log_files):
            try:
                text = lf.read_text(errors="ignore")
                lines = text.splitlines()
                tail = "\n".join(lines[-tail_lines:])
                has_error = any(
                    kw in text.lower()
                    for kw in ("error", "exception", "failed", "killed", "oom")
                )
                if has_error:
                    errors_found.append(str(lf.relative_to(logs_dir)))
                logs_out.append({
                    "file": str(lf.relative_to(logs_dir)),
                    "lines_total": len(lines),
                    "has_error": has_error,
                    "tail": tail,
                })
            except Exception:
                pass

        return {
            "status": "success",
            "n_logs": len(logs_out),
            "errors_found": errors_found,
            "logs": logs_out,
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}
