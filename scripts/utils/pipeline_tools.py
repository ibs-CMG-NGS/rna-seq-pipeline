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
import sys
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

        # Avoid double-appending project_id when the user already included it
        # e.g. results_dir="/results/my-project" + project_id="my-project"
        _rp = Path(results_dir).resolve()
        if _rp.name == project_id:
            project_results = str(_rp)
            base_results = str(_rp.parent)
        else:
            project_results = str(_rp / project_id)
            base_results = str(_rp)

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
        
        # Find all FASTQ files — search recursively into subdirectories.
        # is_file() is REQUIRED: directories can be named *.fastq in per-sample dir layouts
        # (e.g. SampleName_R1_001.fastq/ containing SampleName_R1_001.fastq).
        fastq_files = []
        for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
            fastq_files.extend(f for f in data_path.rglob(ext) if f.is_file())

        if not fastq_files:
            return {
                "status": "error",
                "message": f"No FASTQ files found in {data_dir} (searched recursively)"
            }
        
        # Check if files span multiple subdirectories
        subdirs = set(fq.parent for fq in fastq_files)
        has_subdirs = len(subdirs) > 1 or (len(subdirs) == 1 and next(iter(subdirs)) != data_path)

        # Group by sample
        samples = {}
        read_type = "single-end"
        name_conflicts = []  # track sample_ids where files from different dirs collide

        for fq in fastq_files:
            # Strip known double extensions without corrupting sample names
            # (using Path.stem removes only one suffix, e.g. .gz from .fastq.gz)
            fname = fq.name
            for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
                if fname.lower().endswith(ext):
                    name = fname[:-len(ext)]
                    break
            else:
                name = fq.stem  # fallback: remove only last suffix

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
            elif read_num in samples[sample_id]:
                # Duplicate sample_id from a different subdirectory — flag it
                name_conflicts.append(f"{sample_id} ({fq} vs {samples[sample_id][read_num]})")

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
        
        warnings = []
        if has_subdirs:
            warnings.append(
                f"Files found across {len(subdirs)} subdirectories — "
                "sample sheet mode (use_sample_sheet: true) is recommended."
            )
        if name_conflicts:
            warnings.append(
                f"Duplicate sample names from different subdirectories: {name_conflicts[:3]}"
                + (" ..." if len(name_conflicts) > 3 else "")
            )

        return {
            "status": "success",
            "total_files": len(fastq_files),
            "samples": sample_list,
            "total_samples": len(sample_list),
            "read_type": read_type,
            "total_size_gb": round(total_size, 2),
            "data_dir": str(data_path),
            "has_subdirs": has_subdirs,
            "warnings": warnings,
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
    force: bool = False,
    background: bool = False,
) -> Dict[str, Any]:
    """
    Execute Snakemake workflow.

    Args:
        config_file: Path to config.yaml
        cores: Number of cores to use
        dry_run: If True, show what would be run without executing (always blocking)
        until_rule: Stop at specific rule (e.g., 'fastqc', 'multiqc')
        force: Force re-run even if outputs exist
        background: If True (and dry_run=False), launch in background and return
                    immediately with pid + log_file path. Agent should use True;
                    batch_runner uses False to stay blocking.

    Returns:
        Blocking:   {"status": "success"|"error", "returncode": int, ...}
        Background: {"status": "running", "pid": int, "log_file": str, ...}
        Dry-run:    {"status": "dry_run_ok"|"dry_run_error", "total_jobs": int, ...}
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

        # Resolve snakemake binary. Try multiple locations in priority order:
        # 1. Same conda env as running Python (agent started with correct env)
        # 2. CONDA_PREFIX env var (conda env is activated)
        # 3. Known rna-seq-pipeline conda env (server-specific fallback)
        # 4. System PATH
        import shutil
        _snakemake_candidates = [
            Path(sys.executable).parent / "snakemake",
            Path(os.environ.get("CONDA_PREFIX", "/nonexistent")) / "bin" / "snakemake",
            Path("/home/ygkim/program/anaconda3/envs/rna-seq-pipeline/bin/snakemake"),
        ]
        snakemake_cmd = next(
            (str(p) for p in _snakemake_candidates if p.exists()),
            shutil.which("snakemake") or "snakemake",
        )

        # Use absolute path for config so it works from any cwd
        cmd = [
            snakemake_cmd,
            "--configfile", str(config_path.resolve()),
            "--cores", str(cores)
        ]
        
        if dry_run:
            cmd.append("--dry-run")
        
        if until_rule:
            cmd.extend(["--until", until_rule])
        
        if force:
            cmd.append("--forceall")
        
        # ── Background launch (real run only) ──────────────────────────────────
        if background and not dry_run:
            log_dir = pipeline_root / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / f"snakemake_{config_path.stem}.log"
            pid_file = log_dir / f"snakemake_{config_path.stem}.pid"

            with open(log_file, "w") as lf:
                proc = subprocess.Popen(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    cwd=str(pipeline_root),
                )
            pid_file.write_text(str(proc.pid))

            return {
                "status": "running",
                "pid": proc.pid,
                "log_file": str(log_file),
                "pid_file": str(pid_file),
                "command": " ".join(cmd),
                "note": (
                    "Pipeline launched in background. "
                    "Use monitor_pipeline to check progress, "
                    f"or: tail -f {log_file}"
                ),
            }

        # ── Blocking execution (dry-run or background=False) ────────────────
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
            # Blocking real run: return compact summary
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

def _parse_snakemake_log(log_path: Path) -> Dict[str, Any]:
    """
    Parse a Snakemake log file (from background run or .snakemake/log/).

    Key patterns:
      - "X of Y steps (Z%) done"  → progress
      - "rule NAME:" / "localrule NAME:"  → active rules
      - "Finished job N."  → completed jobs counter
      - "Error in rule NAME:"  → rule-level errors
      - "Exiting because a job execution failed"  → fatal error
      - "Complete log:"  → pipeline finished successfully
    """
    try:
        text = log_path.read_text(errors='ignore')
    except OSError:
        return {}

    # Latest progress: "X of Y steps (Z%) done"
    progress_matches = re.findall(r'(\d+) of (\d+) steps \([\d.]+%\) done', text)
    jobs_done = jobs_total = progress_pct = None
    if progress_matches:
        done_s, total_s = progress_matches[-1]
        jobs_done = int(done_s)
        jobs_total = int(total_s)
        progress_pct = round(jobs_done / jobs_total * 100, 1) if jobs_total else 0.0

    # Rules seen (last one = currently running or most recently finished)
    rule_names = re.findall(r'(?:^|\n)(?:local)?rule (\w+):\n', text)
    current_rule = rule_names[-1] if rule_names else None

    # Rules with at least one completion (rule name precedes "Finished job N.")
    # Build mapping: segment between consecutive rule mentions
    finished_count = len(re.findall(r'Finished job \d+\.', text))

    # Error detection
    error_rules = re.findall(r'Error in rule (\w+):', text)
    fatal_error = 'Exiting because a job execution failed' in text

    # Completion signal
    is_complete = (
        'Complete log:' in text
        or bool(re.search(r'\d+ of \d+ steps \(100%\) done', text))
    )

    log_size_kb = round(log_path.stat().st_size / 1024, 1)

    return {
        "progress_pct": progress_pct,
        "jobs_done": jobs_done,
        "jobs_total": jobs_total,
        "finished_count": finished_count,
        "current_rule": current_rule,
        "error_rules": error_rules,
        "fatal_error": fatal_error,
        "is_complete": is_complete,
        "log_file": str(log_path),
        "log_size_kb": log_size_kb,
    }


def _parse_star_final_log(log_path: Path) -> Dict[str, Any]:
    """
    Parse STAR Log.final.out and return key mapping metrics.
    Format: "   Label | <TAB> value"
    """
    try:
        text = log_path.read_text(errors='ignore')
    except OSError:
        return {}

    def _extract(label: str) -> Optional[str]:
        m = re.search(rf'{re.escape(label)}\s*\|\s*(\S+)', text)
        return m.group(1) if m else None

    uniquely_pct = _extract('Uniquely mapped reads %')
    multi_pct    = _extract('% of reads mapped to multiple loci')
    unmapped_pct = _extract('% of reads unmapped: too short')
    total_reads  = _extract('Number of input reads')

    return {
        "uniquely_mapped_pct": uniquely_pct,
        "multi_mapped_pct": multi_pct,
        "unmapped_too_short_pct": unmapped_pct,
        "total_reads": int(total_reads) if total_reads and total_reads.isdigit() else total_reads,
    }


def _find_snakemake_log(config_file: Optional[str], pipeline_root: Optional[Path]) -> Optional[Path]:
    """
    Find the most recent/relevant Snakemake log:
    1. Background log written by run_pipeline(background=True):
       <pipeline_root>/logs/snakemake_<config_stem>.log
    2. Latest file in <pipeline_root>/.snakemake/log/
    """
    if config_file and pipeline_root:
        bg_log = pipeline_root / "logs" / f"snakemake_{Path(config_file).stem}.log"
        if bg_log.exists() and bg_log.stat().st_size > 0:
            return bg_log

    if pipeline_root:
        smk_log_dir = pipeline_root / ".snakemake" / "log"
        if smk_log_dir.exists():
            logs = sorted(smk_log_dir.glob("*.snakemake.log"), key=lambda p: p.stat().st_mtime)
            if logs:
                return logs[-1]  # most recent

    return None


def monitor_pipeline(
    project_id: str,
    base_results_dir: str,
    config_file: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Check pipeline execution status by parsing Snakemake/STAR logs and
    inspecting output files.

    Priority:
    1. Snakemake log (background log or .snakemake/log/) → exact progress %
    2. STAR Log.final.out per sample → per-sample mapping rates
    3. Output file counting → fallback when no log available

    Returns:
        {
            "status": "not_started"|"running"|"completed"|"error",
            "progress_pct": float,
            "current_stage": str,
            "completed_rules": {rule: n_done},
            "total_rules": {rule: n_expected},
            "n_samples": int,
            "snakemake_log": {...},   # present when log found
            "star_per_sample": {...}, # per-sample mapping stats
            "project_dir": str,
        }
    """
    try:
        base_path = Path(base_results_dir)
        project_dir = base_path / project_id

        # Path overlap fallback: user passed full project path as base_results_dir
        if not project_dir.exists() and base_path.exists():
            candidate_dirs = [
                p for p in base_path.iterdir()
                if p.is_dir()
                and not p.name.startswith('.')
                and p.name not in ('logs', 'project_summary')
            ]
            if candidate_dirs:
                project_dir = base_path

        # Locate pipeline root (walk up from config_file until Snakefile found)
        pipeline_root: Optional[Path] = None
        if config_file and Path(config_file).exists():
            root = Path(config_file).resolve().parent
            for _ in range(6):
                if (root / "Snakefile").exists():
                    pipeline_root = root
                    break
                root = root.parent

        # ── 1. Snakemake log parsing ──────────────────────────────────────────
        smk_log_path = _find_snakemake_log(config_file, pipeline_root)
        smk = _parse_snakemake_log(smk_log_path) if smk_log_path else {}

        # ── 2. STAR Log.final.out parsing ─────────────────────────────────────
        star_logs = sorted(project_dir.glob("intermediate/aligned/*/Log.final.out")) \
            if project_dir.exists() else []
        star_per_sample: Dict[str, Any] = {}
        for sl in star_logs:
            sample_id = sl.parent.name
            star_per_sample[sample_id] = _parse_star_final_log(sl)

        # ── 3. File-based counts (fallback + supplement) ───────────────────────
        # Sample count: from STAR logs first, then config, then dir listing
        n_samples = len(star_logs)
        if n_samples == 0 and config_file and Path(config_file).exists():
            with open(config_file) as f:
                cfg = yaml.safe_load(f)
            data_dir = cfg.get('data_dir', '')
            if data_dir:
                fq = detect_fastq_files(data_dir)
                if fq.get('status') == 'success':
                    n_samples = fq['total_samples']
        if n_samples == 0 and project_dir.exists():
            n_samples = sum(
                1 for p in project_dir.iterdir()
                if p.is_dir() and not p.name.startswith('.')
                and p.name not in ('logs', 'project_summary')
            )

        if n_samples == 0 and not smk:
            return {"status": "not_started",
                    "message": f"No pipeline output found under: {project_dir}"}

        expected = {
            'cutadapt':            n_samples,
            'fastqc_raw':          n_samples * 2,
            'star_align':          n_samples,
            'featurecounts_quant': 1,
            'multiqc':             1,
            'generate_manifest':   n_samples,
        }

        # Count output files
        logs_dir = project_dir / 'logs'
        completed: Dict[str, int] = {}
        completed['star_align']          = len(star_logs)  # most reliable
        completed['generate_manifest']   = len(list(project_dir.glob('*/rna-seq/final_outputs/manifest.json')))
        completed['cutadapt']            = len(list(logs_dir.glob('cutadapt/*.log'))) if logs_dir.exists() else 0
        fqc = list(project_dir.glob('*/rna-seq/final_outputs/qc/*_fastqc.zip'))
        fqc += list((project_dir / 'project_summary' / 'qc').glob('*_fastqc.zip'))
        completed['fastqc_raw']          = len(set(str(f) for f in fqc))
        completed['multiqc']             = len(list(project_dir.glob('project_summary/qc/multiqc_report.html')))
        completed['featurecounts_quant'] = min(len(list(project_dir.glob('project_summary/counts/*.txt'))), 1)

        # ── 4. Determine progress and status ──────────────────────────────────
        # Prefer Snakemake log progress (exact), fallback to file counting
        if smk.get('progress_pct') is not None:
            progress = smk['progress_pct']
        else:
            total_exp  = sum(expected.values())
            total_done = sum(min(completed.get(r, 0), expected[r]) for r in expected)
            progress   = round(total_done / total_exp * 100, 1) if total_exp > 0 else 0.0

        # Current stage: most advanced rule with at least one output
        pipeline_order = ['cutadapt', 'fastqc_raw', 'star_align',
                          'featurecounts_quant', 'generate_manifest', 'multiqc']
        current_stage = smk.get('current_rule') or 'not_started'
        for rule in pipeline_order:
            if completed.get(rule, 0) > 0:
                current_stage = rule
        if completed.get('multiqc', 0) >= 1 or smk.get('is_complete'):
            current_stage = 'completed'

        # Status classification
        has_fatal = smk.get('fatal_error', False)
        is_complete = smk.get('is_complete', False) or progress >= 99.5 or current_stage == 'completed'
        status = (
            'error'       if has_fatal
            else 'completed' if is_complete
            else 'not_started' if progress == 0 and not smk
            else 'running'
        )

        # ── 5. STAR summary stats ─────────────────────────────────────────────
        star_mapping_pcts = []
        for s in star_per_sample.values():
            pct_str = s.get('uniquely_mapped_pct', '')
            if pct_str and pct_str.endswith('%'):
                try:
                    star_mapping_pcts.append(float(pct_str.rstrip('%')))
                except ValueError:
                    pass
        star_summary: Dict[str, Any] = {
            "n_completed": len(star_logs),
            "n_expected": n_samples,
        }
        if star_mapping_pcts:
            star_summary["avg_uniquely_mapped_pct"] = round(sum(star_mapping_pcts) / len(star_mapping_pcts), 2)
            star_summary["min_uniquely_mapped_pct"] = round(min(star_mapping_pcts), 2)

        result: Dict[str, Any] = {
            "status": status,
            "progress_pct": progress,
            "current_stage": current_stage,
            "completed_rules": completed,
            "total_rules": expected,
            "n_samples": n_samples,
            "star_alignment": star_summary,
            "project_dir": str(project_dir),
        }
        if star_per_sample:
            result["star_per_sample"] = star_per_sample
        if smk:
            result["snakemake_log"] = smk
        if smk.get('error_rules'):
            result["error_rules"] = smk['error_rules']

        return result

    except Exception as e:
        return {"status": "error", "message": f"monitor_pipeline failed: {e}"}


def create_sample_sheet(
    project_id: str,
    data_dir: str,
    conditions: Optional[Dict[str, List[str]]] = None,
    output_path: Optional[str] = None,
    overwrite: bool = True,
) -> Dict[str, Any]:
    """
    Generate a sample metadata TSV from detected FASTQ files.

    If `overwrite=False` and the TSV already exists, returns its current content
    without regenerating (preserves manually revised conditions).

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
        tsv_path = output_path or f"config/samples/{project_id}.tsv"

        # If overwrite=False and file already exists, return existing content
        if not overwrite and Path(tsv_path).exists():
            import pandas as _pd
            existing = _pd.read_csv(tsv_path, sep='\t')
            conds = existing['condition'].unique().tolist() if 'condition' in existing.columns else []
            n = len(existing)
            return {
                "status": "success",
                "tsv_path": str(tsv_path),
                "n_samples": n,
                "conditions_assigned": {c: int((existing['condition'] == c).sum()) for c in conds},
                "unassigned": [],
                "skipped": True,
                "message": f"Sample sheet already exists ({n} samples, conditions: {conds}). Skipped regeneration to preserve manual edits.",
            }

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
            # Auto-infer from sample name patterns.
            # Priority: longer/more-specific prefixes first to avoid partial matches
            # (e.g. "nMonTg" must match before "MonTg").
            # suffix_patterns: exact suffix match (e.g. _WT, _Ctrl)
            # prefix_patterns: prefix match on the basename before first digit/_S
            suffix_patterns = {
                'wildtype':      ['_W', '_WT', '_wt', '_wildtype', '_Ctrl', '_ctrl', '_control'],
                'heterozygous':  ['_H', '_Het', '_het', '_heterozygous', '_KO', '_ko'],
                'treatment':     ['_T', '_Treat', '_treat', '_treatment'],
            }
            # Prefix patterns — ordered longest-first to prevent shorter prefixes
            # from shadowing more-specific ones (nMonTg must beat MonTg).
            prefix_patterns = [
                ('non_stimulated', ['nMonTg', 'nTg', 'NonTg', 'non_tg']),
                ('wildtype',       ['WT', 'Wt', 'Wild', 'Ctrl', 'Control', 'MonWT']),
                ('transgenic',     ['MonTg', 'Tg', 'TG']),
                ('heterozygous',   ['Het', 'HET', 'KO', 'ko']),
                ('treatment',      ['Treat', 'treat', 'STIM', 'Stim']),
            ]
            for sample in samples:
                sid = sample['sample_id']
                matched = False
                # 1. Try suffix patterns (exact)
                for cond_name, pats in suffix_patterns.items():
                    for p in pats:
                        if sid.endswith(p):
                            cond_map[sid] = cond_name
                            matched = True
                            break
                    if matched:
                        break
                # 2. Try prefix patterns if suffix didn't match
                if not matched:
                    for cond_name, pats in prefix_patterns:
                        for p in pats:
                            if sid.startswith(p):
                                cond_map[sid] = cond_name
                                matched = True
                                break
                        if matched:
                            break

        # Build TSV rows
        Path(tsv_path).parent.mkdir(parents=True, exist_ok=True)

        # Auto-assign replicate numbers: within each condition, rank by sample_id
        cond_counter: Dict[str, int] = {}
        rep_map: Dict[str, str] = {}
        for sample in sorted(samples, key=lambda s: s['sample_id']):
            sid = sample['sample_id']
            cond = cond_map.get(sid, 'unknown')
            cond_counter[cond] = cond_counter.get(cond, 0) + 1
            rep_map[sid] = str(cond_counter[cond])

        # Check for missing R2 in paired-end data
        missing_r2 = [
            s['sample_id'] for s in samples
            if 'R2' not in s and fq_result.get('read_type') == 'paired-end'
        ]

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
                'replicate': rep_map.get(sid, ''),
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

        warnings = []
        if unassigned:
            warnings.append(
                f"{len(unassigned)} samples unassigned — pass conditions={{...}} to assign them."
            )
        if missing_r2:
            warnings.append(
                f"Missing R2 for {len(missing_r2)} sample(s): {missing_r2[:5]}"
                + (" ..." if len(missing_r2) > 5 else "")
                + " — check FASTQ files or set read_type: single-end."
            )

        return {
            "status": "success",
            "tsv_path": tsv_path,
            "n_samples": len(rows),
            "conditions_assigned": conditions_assigned,
            "unassigned": unassigned,
            "missing_r2": missing_r2,
            "warnings": warnings,
            "message": (
                f"Sample sheet created: {tsv_path}. "
                + (" ".join(warnings) if warnings else "All samples assigned.")
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
            # sample_sheet is relative to the pipeline root.
            # Config files live in config/projects/, so go up 3 levels;
            # fall back to 2 levels for configs placed directly in config/.
            for parent in (cfg_path.parent.parent.parent, cfg_path.parent.parent,
                           cfg_path.parent, Path.cwd()):
                candidate = parent / sheet_path
                if candidate.exists():
                    ss = candidate
                    break
            else:
                ss = cfg_path.parent.parent.parent / sheet_path  # show meaningful error path
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
        summary_dir = project_dir / "project_summary"

        result: Dict[str, Any] = {"status": "success", "project_dir": str(project_dir)}

        # 0. Pipeline QC Evaluation — comprehensive PASS/WARN/FAIL per sample
        pipeline_qc_json = summary_dir / "pipeline_qc_evaluation.json"
        if pipeline_qc_json.exists():
            with open(pipeline_qc_json) as f:
                pqc = json.load(f)
            samples = pqc.get("samples", [])
            n_pass = sum(1 for s in samples if s["status"] == "PASS")
            n_warn = sum(1 for s in samples if s["status"] == "WARN")
            n_fail = sum(1 for s in samples if s["status"] == "FAIL")
            result["pipeline_qc"] = {
                "total_samples": len(samples),
                "passed": n_pass,
                "warned": n_warn,
                "failed": n_fail,
                "pass_rate": round(n_pass / len(samples) * 100, 1) if samples else 0,
                "thresholds": pqc.get("thresholds", {}),
                "per_sample": [
                    {
                        "sample": s["sample"],
                        "status": s["status"],
                        "recommendation": s.get("recommendation", ""),
                        "metrics": s.get("metrics", {}),
                        "issues": s.get("issues", []),
                        "warnings": s.get("warnings", []),
                    }
                    for s in samples
                ],
                "failed_samples": [s["sample"] for s in samples if s["status"] == "FAIL"],
                "warned_samples": [s["sample"] for s in samples if s["status"] == "WARN"],
            }
        else:
            result["pipeline_qc"] = None

        # 1. FastQC evaluation — parse multiqc_fastqc.txt (fastqc_evaluation.json is not generated)
        mqc_fastqc = summary_dir / "multiqc_report_data" / "multiqc_fastqc.txt"
        if not mqc_fastqc.exists():
            mqc_fastqc = next(summary_dir.glob("**/multiqc_fastqc.txt"), None)  # type: ignore
        if mqc_fastqc and Path(mqc_fastqc).exists():
            import csv
            with open(mqc_fastqc) as f:
                reader = csv.DictReader(f, delimiter="\t")
                fastqc_rows = list(reader)
            skip_cols = {"Sample", "Filename", "File type", "Encoding",
                         "Total Sequences", "Total Bases", "Sequences flagged as poor quality",
                         "Sequence length", "%GC", "total_deduplicated_percentage",
                         "avg_sequence_length", "median_sequence_length"}
            # These checks commonly fail for RNA-seq raw data and are NOT quality problems
            expected_fail_checks = {
                "per_base_sequence_content",   # random hexamer priming bias
                "sequence_duplication_levels", # high expression of a few genes
                "adapter_content",             # expected before trimming
                "per_sequence_gc_content",     # RNA-seq GC distribution differs from genome
            }
            fastqc_eval = {}
            for row in fastqc_rows:
                sample = row.get("Sample", "")
                checks = {k: v for k, v in row.items() if k not in skip_cols}
                critical_checks = {k: v for k, v in checks.items()
                                   if k not in expected_fail_checks}
                overall_pass = all(v in ("pass", "warn", "") for v in critical_checks.values())
                expected_fails = [k for k in expected_fail_checks
                                  if checks.get(k) == "fail"]
                fastqc_eval[sample] = {
                    "overall_pass": overall_pass,
                    "checks": checks,
                    "expected_fails": expected_fails,
                }
            result["fastqc_evaluation"] = fastqc_eval
        else:
            result["fastqc_evaluation"] = None

        # 2. MultiQC general stats TSV
        mqc_stats = summary_dir / "multiqc_report_data" / "multiqc_general_stats.txt"
        if not mqc_stats.exists():
            mqc_stats = next(summary_dir.glob("**/multiqc_general_stats.txt"), None)  # type: ignore

        if mqc_stats and Path(mqc_stats).exists():
            import csv
            with open(mqc_stats) as f:
                reader = csv.DictReader(f, delimiter="\t")
                mqc_rows = list(reader)
            result["multiqc_stats"] = mqc_rows
            # multiqc_general_stats has one row per tool×sample
            # (BAM paths for STAR/featureCounts + sample names for FastQC/Cutadapt).
            # Count only unique biological samples: non-BAM rows, strip _1/_2 read suffix.
            import re as _re
            sample_names = set()
            for row in mqc_rows:
                s = row.get("Sample", "")
                if ".bam" in s:
                    continue  # BAM path rows — skip
                s = _re.sub(r'[_-][12]$', '', s)  # strip _1/_2 read suffix
                if s:
                    sample_names.add(s)
            result["multiqc_n_samples"] = len(sample_names) if sample_names else len(mqc_rows)
        else:
            result["multiqc_stats"] = None

        # 3. STAR logs — extract uniquely mapped %
        star_logs = list(project_dir.glob("intermediate/aligned/*/Log.final.out"))
        if not star_logs:
            star_logs = list(project_dir.glob("*/rna-seq/final_outputs/*Log.final.out"))
        if not star_logs:
            star_logs = list(project_dir.glob("logs/star/*.log"))
        star_summary: List[Dict] = []
        for log in star_logs[:38]:
            try:
                text = Path(log).read_text()
                uniq_m = re.search(r'Uniquely mapped reads %\s+\|\s+([\d.]+)%', text)
                multi_m = re.search(r'% of reads mapped to multiple loci\s+\|\s+([\d.]+)%', text)
                unmapped_m = re.search(r'% of reads unmapped: too short\s+\|\s+([\d.]+)%', text)
                if "intermediate/aligned" in str(log):
                    sample_name = log.parent.name
                elif 'rna-seq' in str(log):
                    sample_name = log.parent.parent.parent.name
                else:
                    sample_name = log.stem
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

        # Overall summary — count at file level (R1+R2) and sample level
        n_pass_files = n_fail_files = 0
        if result["fastqc_evaluation"] and isinstance(result["fastqc_evaluation"], dict):
            for v in result["fastqc_evaluation"].values():
                if isinstance(v, dict):
                    if v.get("overall_pass"):
                        n_pass_files += 1
                    else:
                        n_fail_files += 1
        n_fastqc_files = n_pass_files + n_fail_files

        pqc = result.get("pipeline_qc") or {}
        result["summary"] = {
            "pipeline_qc_total": pqc.get("total_samples"),
            "pipeline_qc_passed": pqc.get("passed"),
            "pipeline_qc_warned": pqc.get("warned"),
            "pipeline_qc_failed": pqc.get("failed"),
            "pipeline_qc_pass_rate": f"{pqc.get('pass_rate', 'N/A')}%",
            "pipeline_qc_failed_samples": pqc.get("failed_samples", []),
            "fastqc_files_evaluated": n_fastqc_files,
            "fastqc_critical_pass": n_pass_files,
            "fastqc_critical_fail": n_fail_files,
            "fastqc_note": (
                "per_base_sequence_content, sequence_duplication_levels, adapter_content "
                "failures are expected for RNA-seq raw data and excluded from pass/fail."
            ),
            "multiqc_n_samples": result.get("multiqc_n_samples"),
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


# ─────────────────────────────────────────────────────────────────────────────
# DE analysis bridge
# ─────────────────────────────────────────────────────────────────────────────

def run_bridge(
    config_file: str,
    de_pipeline_dir: str,
    skip_de: bool = True,
    dry_run: bool = False,
    exclude_samples: List[str] = None,
) -> Dict[str, Any]:
    """
    Transfer RNA-seq counts matrix to DE/GO analysis pipeline.

    Reads rnaseq output dir from config_file (results_dir or base_results_dir/project_id),
    then runs PipelineBridge to:
      1. Verify counts matrix exists
      2. Copy counts to DE pipeline input
      3. Generate metadata TSV from sample sheet
      4. Generate DE pipeline config.yml
      5. Optionally trigger DE/GO Snakemake (skip_de=False)

    Args:
        config_file:     Path to RNA-seq project config YAML
        de_pipeline_dir: Path to DE/GO analysis pipeline directory
        skip_de:         If True (default), only prepare files — do NOT run DE pipeline
        dry_run:         If True, show what would be done without writing files

    Returns:
        {
            "status": "success" | "error",
            "rnaseq_output": str,
            "de_pipeline": str,
            "counts_file": str | None,
            "metadata_file": str | None,
            "de_config_file": str | None,
            "message": str
        }
    """
    try:
        import yaml as _yaml

        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        with open(cfg_path) as f:
            cfg = _yaml.safe_load(f)

        project_id = cfg.get("project_id", "project")

        # Resolve RNA-seq output directory (same logic as Snakefile)
        if "results_dir" in cfg:
            rnaseq_output = Path(cfg["results_dir"])
        else:
            base = cfg.get("base_results_dir", "results")
            rnaseq_output = Path(base) / project_id

        de_pipeline = Path(de_pipeline_dir)

        if not rnaseq_output.exists():
            return {
                "status": "error",
                "message": f"RNA-seq output not found: {rnaseq_output}. "
                           f"Has the pipeline completed?",
                "rnaseq_output": str(rnaseq_output),
            }

        if not de_pipeline.exists():
            return {
                "status": "error",
                "message": f"DE/GO pipeline directory not found: {de_pipeline}",
                "de_pipeline": str(de_pipeline),
            }

        # Import PipelineBridge from scripts/bridge_to_de_pipeline.py
        bridge_script = Path(__file__).parent.parent / "bridge_to_de_pipeline.py"
        import importlib.util
        spec = importlib.util.spec_from_file_location("bridge_to_de_pipeline", bridge_script)
        bridge_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(bridge_mod)
        PipelineBridge = bridge_mod.PipelineBridge

        bridge = PipelineBridge(
            rnaseq_output_dir=rnaseq_output,
            de_pipeline_dir=de_pipeline,
            project_id=project_id,
            config={
                "sample_sheet_dir": str(cfg_path.parent.parent / "samples"),
                "counts_relpath": cfg.get(
                    "counts_relpath",
                    "project_summary/counts/counts_matrix_clean.csv"
                ),
                "species": cfg.get("species"),
            },
            assume_yes=True,
            exclude_samples=exclude_samples or [],
        )

        # Step 1: check completion
        if not bridge.check_rnaseq_completion():
            return {
                "status": "error",
                "message": "RNA-seq pipeline not complete. counts_matrix_clean.csv not found.",
                "rnaseq_output": str(rnaseq_output),
            }

        result: Dict[str, Any] = {
            "status": "success",
            "rnaseq_output": str(rnaseq_output),
            "de_pipeline": str(de_pipeline),
            "project_id": project_id,
            "dry_run": dry_run,
        }

        if dry_run:
            counts_relpath = cfg.get(
                "counts_relpath",
                "project_summary/counts/counts_matrix_clean.csv"
            )
            counts_file = rnaseq_output / counts_relpath
            excl_note = f", excluding {exclude_samples}" if exclude_samples else ""
            result["message"] = (
                f"DRY RUN — would copy {counts_file} to {de_pipeline}/data/raw/{excl_note}, "
                f"generate metadata and DE config for project '{project_id}'."
            )
            result["counts_exists"] = counts_file.exists()
            result["excluded_samples"] = exclude_samples or []
            return result

        # Step 2: copy counts
        counts_dest = bridge.prepare_de_input()
        result["counts_file"] = str(counts_dest)

        # Step 3: metadata
        metadata_dest = bridge.generate_metadata_template()
        result["metadata_file"] = str(metadata_dest) if metadata_dest else None

        # Step 4: DE config
        if metadata_dest:
            de_config = bridge.generate_de_config(counts_dest, metadata_dest)
            result["de_config_file"] = str(de_config)
        else:
            result["de_config_file"] = None

        # Step 5: optionally run DE
        if not skip_de and metadata_dest and de_config:
            bridge.trigger_de_analysis(dry_run=False, config_file=de_config)
            result["de_triggered"] = True
        else:
            result["de_triggered"] = False

        result["message"] = (
            f"Bridge complete. counts → {result['counts_file']}, "
            f"metadata → {result['metadata_file']}, "
            f"DE config → {result['de_config_file']}. "
            + ("DE analysis triggered." if result["de_triggered"]
               else "DE analysis NOT triggered (skip_de=True). Run manually when ready.")
        )
        return result

    except Exception as e:
        import traceback
        return {"status": "error", "message": str(e), "traceback": traceback.format_exc()}


def get_pipeline_status(
    config_file: Optional[str] = None,
    data_dir: Optional[str] = None,
    results_dir: Optional[str] = None,
    de_pipeline_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    파일시스템을 읽어 RNA-seq → DE/GO → seqviewer 파이프라인의 현재 진행 상태를 추론한다.

    별도 done flag를 새로 만들지 않고 기존 결과물(config, sample_sheet, counts matrix,
    final_de_results.csv, .seqviewer_done.flag 등)의 존재 여부로 각 단계 완료를 판단한다.

    Args:
        config_file:     RNA-seq 프로젝트 config YAML 경로 (optional)
        data_dir:        FASTQ 원본 디렉터리 경로 (optional)
        results_dir:     RNA-seq 결과 출력 디렉터리 (optional; config에서 자동 파싱)
        de_pipeline_dir: DE/GO 분석 파이프라인 루트 경로 (optional)

    Returns:
        {
          "project_id": str | None,
          "steps": {step_id: {"label": str, "done": bool}},
          "completed": [step_id, ...],
          "next_step": step_id | None,
          "next_suggestion": str,
          "summary": str   # 예: "5/8 단계 완료 ■■■■■□□□"
        }
    """
    # ── 1. config 파싱으로 경로 보완 ──────────────────────────────
    project_id = None
    sample_sheet_path = None
    config_path_obj = Path(config_file) if config_file else None

    # pipeline root = /data_3tb/shared/rna-seq-pipeline  (scripts/utils/ → 2 levels up)
    pipeline_root = Path(__file__).parent.parent.parent

    if config_path_obj and config_path_obj.exists():
        try:
            import yaml as _yaml
            cfg = _yaml.safe_load(config_path_obj.read_text())
            project_id = cfg.get("project_id")
            if not results_dir:
                results_dir = cfg.get("results_dir") or cfg.get("base_results_dir")
            # data_dir: config의 data_dir 필드에서 파싱 (파라미터 우선)
            if not data_dir:
                data_dir = cfg.get("data_dir")
            # sample_sheet 경로: 상대 경로는 pipeline root 기준
            ss_raw = cfg.get("sample_sheet")
            if ss_raw:
                ss_p = Path(ss_raw)
                if not ss_p.is_absolute():
                    ss_p = pipeline_root / ss_p
                sample_sheet_path = ss_p
        except Exception:
            pass

    results_dir_p = Path(results_dir) if results_dir else None
    de_dir_p = Path(de_pipeline_dir) if de_pipeline_dir else None

    # ── 2. 각 단계 판단 ───────────────────────────────────────────
    def _has_fastq(d: Optional[Path]) -> bool:
        if not d or not d.exists():
            return False
        return any(d.rglob("*.fastq.gz")) or any(d.rglob("*.fq.gz")) or \
               any(d.rglob("*.fastq")) or any(d.rglob("*.fq"))

    def _pipeline_started(pid: Optional[str]) -> bool:
        """pipeline root의 logs/ 디렉터리에서 snakemake log 존재 여부 확인"""
        logs_dir = pipeline_root / "logs"
        if not logs_dir.exists():
            return False
        if pid:
            return (logs_dir / f"snakemake_{pid}.log").exists()
        return any(logs_dir.glob("snakemake_*.log"))

    def _pipeline_done(res: Optional[Path]) -> bool:
        """counts_matrix_clean.csv 존재 → RNA-seq 완료"""
        if not res or not res.exists():
            return False
        return (res / "project_summary" / "counts" / "counts_matrix_clean.csv").exists()

    def _de_prepared(de: Optional[Path], pid: Optional[str]) -> bool:
        """DE pipeline data/raw/에 counts CSV 존재 → bridge 완료"""
        if not de or not de.exists():
            return False
        raw = de / "data" / "raw"
        if pid:
            return (raw / f"{pid}_counts.csv").exists()
        return any(raw.glob("*_counts.csv"))

    def _de_done(de: Optional[Path], pid: Optional[str]) -> bool:
        """pairwise/*/final_de_results.csv 존재 → DE 분석 완료"""
        if not de or not de.exists():
            return False
        out = de / "output"
        if pid:
            out = out / pid
        if not out.exists():
            return False
        return any(out.rglob("final_de_results.csv"))

    def _seqviewer_done(de: Optional[Path], pid: Optional[str]) -> bool:
        if not de or not de.exists():
            return False
        out = de / "output"
        if pid:
            out = out / pid
        return (out / "seqviewer" / ".seqviewer_done.flag").exists()

    data_dir_p = Path(data_dir) if data_dir else None

    STEPS = [
        ("fastq_detected",      "FASTQ 파일 탐지",         _has_fastq(data_dir_p)),
        ("config_created",      "프로젝트 설정(config) 생성", bool(config_path_obj and config_path_obj.exists())),
        ("sample_sheet_created","샘플시트 생성",             bool(sample_sheet_path and sample_sheet_path.exists())),
        ("pipeline_started",    "RNA-seq 파이프라인 시작",   _pipeline_started(project_id)),
        ("pipeline_done",       "RNA-seq 파이프라인 완료",   _pipeline_done(results_dir_p)),
        ("de_prepared",         "DE 분석 준비 (bridge)",    _de_prepared(de_dir_p, project_id)),
        ("de_done",             "DE/GO 분석 완료",           _de_done(de_dir_p, project_id)),
        ("seqviewer_done",      "seqviewer 내보내기 완료",   _seqviewer_done(de_dir_p, project_id)),
    ]

    SUGGESTIONS = {
        "fastq_detected":       "FASTQ 파일을 탐지하세요. (detect_fastq_files)",
        "config_created":       "프로젝트 config를 생성하세요. (create_project_config)",
        "sample_sheet_created": "샘플시트를 생성하세요. (create_sample_sheet)",
        "pipeline_started":     "RNA-seq 파이프라인을 실행하세요. (run_pipeline)",
        "pipeline_done":        "파이프라인 실행 중입니다. 완료를 기다리세요. (monitor_pipeline)",
        "de_prepared":          "DE 분석 파이프라인으로 결과를 전달하세요. (run_bridge)",
        "de_done":              "DE/GO 분석을 실행하세요. (DE pipeline Snakemake)",
        "seqviewer_done":       "seqviewer export를 실행하세요. (DE pipeline seqviewer export)",
    }

    steps_dict = {sid: {"label": lbl, "done": done} for sid, lbl, done in STEPS}
    completed = [sid for sid, _, done in STEPS if done]

    # 첫 번째 미완료 단계 찾기
    next_step = None
    next_suggestion = "모든 단계가 완료되었습니다."
    for sid, _, done in STEPS:
        if not done:
            next_step = sid
            next_suggestion = SUGGESTIONS.get(sid, f"{sid} 단계를 진행하세요.")
            break

    n_done = len(completed)
    n_total = len(STEPS)
    bar = "■" * n_done + "□" * (n_total - n_done)
    summary = f"{n_done}/{n_total} 단계 완료 [{bar}]"

    return {
        "project_id": project_id,
        "steps": steps_dict,
        "completed": completed,
        "next_step": next_step,
        "next_suggestion": next_suggestion,
        "summary": summary,
    }


# ── RSeQC Library Validation ──────────────────────────────────────────────────

def _gtf_to_bed12(gtf_path: str, bed_path: str) -> None:
    """
    Convert GTF annotation to BED12 format required by RSeQC.
    Groups exons by transcript_id, outputs one BED12 line per transcript.
    """
    import re

    def _attr(attrs: str, key: str) -> str:
        m = re.search(rf'{key}\s+"([^"]+)"', attrs)
        return m.group(1) if m else ""

    # Collect exons grouped by transcript_id
    transcripts: Dict[str, Dict] = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("exon", "transcript"):
                continue
            chrom, start, end, strand, attrs = (
                parts[0], int(parts[3]) - 1, int(parts[4]), parts[6], parts[8]
            )
            tid = _attr(attrs, "transcript_id")
            if not tid:
                continue
            if feature == "transcript":
                transcripts.setdefault(tid, {"chrom": chrom, "start": start,
                                             "end": end, "strand": strand,
                                             "name": _attr(attrs, "gene_name") or tid,
                                             "exons": []})
                transcripts[tid].update({"chrom": chrom, "start": start,
                                         "end": end, "strand": strand})
            else:  # exon
                transcripts.setdefault(tid, {"chrom": chrom, "start": start,
                                             "end": end, "strand": strand,
                                             "name": tid, "exons": []})
                transcripts[tid]["exons"].append((start, end))

    with open(bed_path, "w") as out:
        for tid, tx in transcripts.items():
            exons = sorted(tx["exons"])
            if not exons:
                continue
            tx_start = tx["start"]
            tx_end = tx["end"]
            block_count = len(exons)
            block_sizes = ",".join(str(e - s) for s, e in exons) + ","
            block_starts = ",".join(str(s - tx_start) for s, e in exons) + ","
            out.write(
                f"{tx['chrom']}\t{tx_start}\t{tx_end}\t{tx['name']}\t0\t"
                f"{tx['strand']}\t{tx_start}\t{tx_end}\t0\t"
                f"{block_count}\t{block_sizes}\t{block_starts}\n"
            )


def validate_library_type(
    config_file: str,
    n_reads: int = 500_000,
    cores: int = 8,
    sample_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Pre-flight mRNA-seq library validation using RSeQC.

    Steps:
    1. Convert GTF → BED12 (cached at genome_dir/annotation_rseqc.bed)
    2. Select first sample, subset first n_reads pairs from FASTQ
    3. Run STAR alignment on the subset
    4. Run infer_experiment.py → auto-detect strandedness
    5. Run read_distribution.py → verify mRNA-seq (exon enrichment)
    6. Update strandedness in config YAML if changed
    7. Clean up temporary files

    Returns:
        {
            "status": "success"|"warning"|"error",
            "strandedness_detected": 0|1|2,
            "strandedness_updated": bool,
            "strandedness_before": int,
            "library_type_verdict": "mRNA-seq"|"suspicious"|"likely_not_mRNA-seq",
            "exon_pct": float,
            "intergenic_pct": float,
            "infer_experiment_raw": str,
            "read_distribution_raw": str,
            "warnings": [...],
            "sample_used": str,
            "bed_file": str,
        }
    """
    import tempfile
    import shutil as _shutil

    warnings_list: List[str] = []

    try:
        # ── Load config ────────────────────────────────────────────────────
        with open(config_file) as f:
            cfg = yaml.safe_load(f)

        annotation_gtf = cfg.get("annotation_gtf") or cfg.get("genes_gtf", "")
        star_index     = cfg.get("star_index", "")
        genome_dir     = cfg.get("genome_dir", "") or str(Path(star_index).parent)
        strandedness_before = int(cfg.get("strandedness", 0))

        if not annotation_gtf or not Path(annotation_gtf).exists():
            return {"status": "error", "message": f"annotation_gtf not found: {annotation_gtf}"}
        if not star_index or not Path(star_index).exists():
            return {"status": "error", "message": f"star_index not found: {star_index}"}

        # ── Resolve tool paths (same conda env as this Python) ─────────────
        _conda_bin = Path(sys.executable).parent
        def _tool(name: str) -> str:
            p = _conda_bin / name
            return str(p) if p.exists() else name

        infer_exp_bin    = _tool("infer_experiment.py")
        read_dist_bin    = _tool("read_distribution.py")
        star_bin         = _tool("STAR")
        samtools_bin     = _tool("samtools")

        for bin_path, label in [
            (infer_exp_bin, "infer_experiment.py"),
            (read_dist_bin, "read_distribution.py"),
            (star_bin, "STAR"),
            (samtools_bin, "samtools"),
        ]:
            if not Path(bin_path).exists() and _shutil.which(bin_path) is None:
                return {"status": "error",
                        "message": f"{label} not found. Install RSeQC: conda install -c bioconda rseqc"}

        # ── Step 1: GTF → BED12 (cache) ────────────────────────────────────
        bed_file = str(Path(genome_dir) / "annotation_rseqc.bed")
        if not Path(bed_file).exists():
            print(f"[validate_library_type] Converting GTF → BED12: {bed_file}")
            _gtf_to_bed12(annotation_gtf, bed_file)
        else:
            print(f"[validate_library_type] Using cached BED12: {bed_file}")

        # ── Step 2: Select sample ──────────────────────────────────────────
        # If sample_id is specified, find that sample; otherwise use first available.
        sample_sheet = cfg.get("sample_sheet", "")
        r1_path = r2_path = ""
        selected_id = ""

        if sample_sheet and Path(sample_sheet).exists():
            import csv
            with open(sample_sheet) as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
            # Find target row: match sample_id if given, else use first valid row
            for row in rows:
                sid = row.get("sample_id", "")
                if sample_id and sid != sample_id:
                    continue
                r1 = row.get("fastq_r1", "")
                r2 = row.get("fastq_r2", "")
                if r1 and Path(r1).exists():
                    r1_path, r2_path, selected_id = r1, r2, sid
                    break
            if sample_id and not selected_id:
                return {"status": "error",
                        "message": f"sample_id '{sample_id}' not found in sample sheet or FASTQ missing"}
        else:
            fq_result = detect_fastq_files(cfg.get("data_dir", ""))
            samples = fq_result.get("samples", [])
            target = next((s for s in samples if not sample_id or s["sample_id"] == sample_id), None)
            if not target:
                return {"status": "error",
                        "message": f"sample_id '{sample_id}' not found in data_dir"}
            r1_path = target.get("R1", "")
            r2_path = target.get("R2", "")
            selected_id = target.get("sample_id", "sample1")

        if not r1_path or not Path(r1_path).exists():
            return {"status": "error", "message": f"Cannot find R1 FASTQ for validation (path: {r1_path})"}

        print(f"[validate_library_type] Using sample: {selected_id}")

        # ── Step 3: Subset FASTQ ────────────────────────────────────────────
        tmp_dir = Path(tempfile.mkdtemp(prefix="rseqc_validate_"))
        try:
            n_lines = n_reads * 4  # 4 lines per FASTQ record
            sub_r1 = str(tmp_dir / "subset_R1.fastq")
            sub_r2 = str(tmp_dir / "subset_R2.fastq") if r2_path else None

            def _subset_fastq(src: str, dst: str, n: int) -> None:
                opener = "zcat" if src.endswith(".gz") else "cat"
                subprocess.run(
                    f"{opener} {src} | head -{n} > {dst}",
                    shell=True, check=True
                )

            print(f"[validate_library_type] Subsetting {n_reads:,} reads...")
            _subset_fastq(r1_path, sub_r1, n_lines)
            if sub_r2:
                _subset_fastq(r2_path, sub_r2, n_lines)

            # ── Step 4: STAR alignment on subset ───────────────────────────
            # STAR requires trailing slash in prefix to write files inside directory
            star_out_dir = tmp_dir / "star_out"
            star_out_dir.mkdir(parents=True, exist_ok=True)
            star_out_prefix = str(star_out_dir) + "/"

            star_cmd = [
                star_bin,
                "--runThreadN", str(cores),
                "--genomeDir", star_index,
                "--readFilesIn", sub_r1,
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outFileNamePrefix", star_out_prefix,
                "--outSAMattributes", "NH", "HI", "AS", "NM",
                "--limitBAMsortRAM", "4000000000",
            ]
            if sub_r2:
                star_cmd[star_cmd.index(sub_r1) + 1:star_cmd.index(sub_r1) + 1] = [sub_r2]

            print(f"[validate_library_type] Running STAR on subset...")
            star_result = subprocess.run(
                star_cmd, capture_output=True, text=True, cwd=str(tmp_dir)
            )
            if star_result.returncode != 0:
                return {
                    "status": "error",
                    "message": "STAR alignment failed",
                    "stderr": star_result.stderr[-1000:],
                }

            bam_file = str(star_out_dir / "Aligned.sortedByCoord.out.bam")
            subprocess.run([samtools_bin, "index", bam_file], check=True)

            # ── Step 5: infer_experiment.py ────────────────────────────────
            print(f"[validate_library_type] Running infer_experiment.py...")
            infer_result = subprocess.run(
                [infer_exp_bin, "-r", bed_file, "-i", bam_file],
                capture_output=True, text=True
            )
            infer_raw = infer_result.stdout + infer_result.stderr

            # Parse strandedness from output
            # RSeQC outputs fractions for 1++,1--,2+-,2-+ patterns
            frac_failed = frac_fwd = frac_rev = 0.0
            for line in infer_raw.splitlines():
                if "failed to determine" in line.lower():
                    m = re.search(r"([\d.]+)$", line.strip())
                    if m:
                        frac_failed = float(m.group(1))
                elif "1++,1--,2+-,2-+" in line:
                    m = re.search(r"([\d.]+)$", line.strip())
                    if m:
                        frac_fwd = float(m.group(1))
                elif "1+-,1-+,2++,2--" in line:
                    m = re.search(r"([\d.]+)$", line.strip())
                    if m:
                        frac_rev = float(m.group(1))

            if frac_failed > 0.75 or (frac_fwd < 0.6 and frac_rev < 0.6):
                strandedness_detected = 0  # unstranded
            elif frac_fwd >= 0.6:
                strandedness_detected = 1  # forward
            else:
                strandedness_detected = 2  # reverse (e.g. TruSeq)

            # ── Step 6: read_distribution.py ───────────────────────────────
            print(f"[validate_library_type] Running read_distribution.py...")
            dist_result = subprocess.run(
                [read_dist_bin, "-r", bed_file, "-i", bam_file],
                capture_output=True, text=True
            )
            dist_raw = dist_result.stdout + dist_result.stderr

            # Parse exon/intergenic percentages
            # Output format: Group  Total_bases  Tag_count  Tags/Kb
            # Tag_count is column index 2 (0-based) after splitting by whitespace
            total_tags = total_assigned = exon_tags = 0
            for line in dist_raw.splitlines():
                parts = line.split()
                if "Total Tags" in line and len(parts) >= 3:
                    total_tags = int(parts[2])
                elif "Total Assigned Tags" in line and len(parts) >= 4:
                    total_assigned = int(parts[3])
                elif parts and parts[0] in ("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"):
                    if len(parts) >= 3:
                        exon_tags += int(parts[2])

            # Intergenic = unassigned reads (not in any annotated region)
            intergenic_tags = total_tags - total_assigned if total_tags > total_assigned else 0
            exon_pct       = round(exon_tags / total_tags * 100, 1) if total_tags else 0.0
            intergenic_pct = round(intergenic_tags / total_tags * 100, 1) if total_tags else 0.0

            if exon_pct >= 60:
                library_verdict = "mRNA-seq"
            elif exon_pct >= 30:
                library_verdict = "suspicious"
                warnings_list.append(
                    f"Exon coverage {exon_pct}% is lower than expected for mRNA-seq (≥60%). "
                    "Could be pre-mRNA, low-quality, or ATAC-seq/WGS."
                )
            else:
                library_verdict = "likely_not_mRNA-seq"
                warnings_list.append(
                    f"Exon coverage {exon_pct}% is very low — likely NOT mRNA-seq. "
                    f"Intergenic: {intergenic_pct}%. Consider checking if data is ATAC-seq or WGS."
                )

            if intergenic_pct > 40:
                warnings_list.append(
                    f"Intergenic reads {intergenic_pct}% is high — possible ATAC-seq or contamination."
                )

            # ── Step 7: Update config if strandedness changed ──────────────
            strandedness_updated = strandedness_detected != strandedness_before
            if strandedness_updated:
                with open(config_file) as f:
                    config_text = f.read()
                config_text = re.sub(
                    r"^(strandedness:\s*)(\d+)",
                    rf"\g<1>{strandedness_detected}",
                    config_text,
                    flags=re.MULTILINE,
                )
                with open(config_file, "w") as f:
                    f.write(config_text)
                warnings_list.append(
                    f"strandedness updated: {strandedness_before} → {strandedness_detected} in {config_file}"
                )

        finally:
            _shutil.rmtree(tmp_dir, ignore_errors=True)

        status = "warning" if warnings_list else "success"
        return {
            "status": status,
            "strandedness_detected": strandedness_detected,
            "strandedness_updated": strandedness_updated,
            "strandedness_before": strandedness_before,
            "library_type_verdict": library_verdict,
            "exon_pct": exon_pct,
            "intergenic_pct": intergenic_pct,
            "infer_experiment_raw": infer_raw,
            "read_distribution_raw": dist_raw,
            "warnings": warnings_list,
            "sample_used": selected_id,
            "bed_file": bed_file,
        }

    except Exception as e:
        return {"status": "error", "message": f"validate_library_type failed: {e}"}


def validate_de_config_conditions(de_config_file: str) -> Dict[str, Any]:
    """
    Validate that pairwise_comparisons in a DE config YAML match the actual
    condition values present in the referenced metadata CSV.

    Detects two common mismatches:
    - Case difference:          'control' vs 'Control'
    - Alphanumeric order swap:  '1D' vs 'D1', '3D' vs 'D3'

    Returns:
        {
            "status": "ok" | "mismatch" | "error",
            "metadata_conditions": [str, ...],
            "config_conditions":   [str, ...],
            "mismatches":          [str, ...],   # in config but not in metadata
            "suggestions":         {wrong: correct, ...},
            "unresolved":          [str, ...],   # no automatic fix found
            "message":             str
        }
    """
    try:
        import re as _re
        import yaml as _yaml
        import pandas as _pd

        de_cfg_path = Path(de_config_file)
        if not de_cfg_path.exists():
            return {"status": "error", "message": f"DE config not found: {de_config_file}"}

        with open(de_cfg_path) as f:
            de_cfg = _yaml.safe_load(f)

        # Resolve metadata path (relative to DE pipeline root = parent of configs/)
        metadata_raw = de_cfg.get("metadata_path", "")
        if not Path(metadata_raw).is_absolute():
            metadata_path = de_cfg_path.parent.parent / metadata_raw
        else:
            metadata_path = Path(metadata_raw)

        if not metadata_path.exists():
            return {"status": "error", "message": f"Metadata not found: {metadata_path}"}

        meta_df = _pd.read_csv(metadata_path)
        cond_col = "condition" if "condition" in meta_df.columns else "group"
        metadata_conditions = set(meta_df[cond_col].unique().tolist())

        # Collect all condition names used in pairwise_comparisons
        pairwise = de_cfg.get("de_analysis", {}).get("pairwise_comparisons", [])
        config_conditions: set = set()
        for pair in pairwise:
            config_conditions.update(pair)

        mismatches = config_conditions - metadata_conditions
        if not mismatches:
            return {
                "status": "ok",
                "metadata_conditions": sorted(metadata_conditions),
                "config_conditions": sorted(config_conditions),
                "mismatches": [],
                "suggestions": {},
                "unresolved": [],
                "message": (
                    f"✅ All conditions match.\n"
                    f"  Metadata: {sorted(metadata_conditions)}\n"
                    f"  Config:   {sorted(config_conditions)}"
                ),
            }

        # Build suggestions for each mismatch
        suggestions: Dict[str, str] = {}
        for m in sorted(mismatches):
            # 1. Case-insensitive exact match  (e.g., 'control' → 'Control')
            ci = next((c for c in metadata_conditions if c.lower() == m.lower()), None)
            if ci:
                suggestions[m] = ci
                continue
            # 2. Alphanumeric order swap  (e.g., '1D' ↔ 'D1', '3D' ↔ 'D3')
            nums  = "".join(_re.findall(r"\d+", m))
            chars = "".join(_re.findall(r"[A-Za-z]+", m))
            if nums and chars:
                for cand in [chars + nums, nums + chars]:
                    hit = next((c for c in metadata_conditions
                                if c == cand or c.lower() == cand.lower()), None)
                    if hit:
                        suggestions[m] = hit
                        break

        unresolved = sorted(m for m in mismatches if m not in suggestions)

        lines = ["⚠️  Condition mismatch in DE config:"]
        lines.append(f"  Metadata conditions : {sorted(metadata_conditions)}")
        lines.append(f"  Config conditions   : {sorted(config_conditions)}")
        lines.append(f"  Mismatched (config) : {sorted(mismatches)}")
        if suggestions:
            lines.append("  Suggested corrections:")
            for wrong, correct in sorted(suggestions.items()):
                lines.append(f"    '{wrong}' → '{correct}'")
        if unresolved:
            lines.append(f"  Unresolved (manual fix needed): {unresolved}")

        return {
            "status": "mismatch",
            "metadata_conditions": sorted(metadata_conditions),
            "config_conditions": sorted(config_conditions),
            "mismatches": sorted(mismatches),
            "suggestions": suggestions,
            "unresolved": unresolved,
            "message": "\n".join(lines),
        }

    except Exception as e:
        return {"status": "error", "message": f"validate_de_config_conditions failed: {e}"}


def apply_de_config_corrections(
    de_config_file: str,
    corrections: Dict[str, str],
) -> Dict[str, Any]:
    """
    Apply condition name corrections to pairwise_comparisons in a DE config YAML.

    Args:
        de_config_file: Path to the DE config YAML (will be updated in-place).
        corrections:    {wrong_name: correct_name, ...}  — from validate_de_config_conditions.

    Returns:
        {"status": "success"|"error", "applied": {wrong: correct}, "message": str}
    """
    try:
        import yaml as _yaml

        de_cfg_path = Path(de_config_file)
        if not de_cfg_path.exists():
            return {"status": "error", "message": f"DE config not found: {de_config_file}"}

        with open(de_cfg_path) as f:
            content = f.read()
            de_cfg = _yaml.safe_load(content)

        pairwise = de_cfg.get("de_analysis", {}).get("pairwise_comparisons", [])
        if not pairwise:
            return {"status": "error", "message": "No pairwise_comparisons found in DE config."}

        applied: Dict[str, str] = {}
        new_pairwise = []
        for pair in pairwise:
            new_pair = []
            for cond in pair:
                if cond in corrections:
                    new_pair.append(corrections[cond])
                    applied[cond] = corrections[cond]
                else:
                    new_pair.append(cond)
            new_pairwise.append(new_pair)

        de_cfg["de_analysis"]["pairwise_comparisons"] = new_pairwise

        with open(de_cfg_path, "w") as f:
            _yaml.dump(de_cfg, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

        lines = [f"✅ DE config updated: {de_cfg_path.name}"]
        for wrong, correct in sorted(applied.items()):
            lines.append(f"  '{wrong}' → '{correct}'")
        lines.append("Re-run DE pipeline to apply changes.")

        return {
            "status": "success",
            "applied": applied,
            "de_config_file": str(de_cfg_path),
            "message": "\n".join(lines),
        }

    except Exception as e:
        return {"status": "error", "message": f"apply_de_config_corrections failed: {e}"}


def read_de_results_summary(
    de_config_file: str,
    pair: Optional[str] = None,
    top_n: int = 10,
) -> Dict[str, Any]:
    """
    DE 분석 결과 요약 — LLM이 해석할 수 있는 구조화된 데이터 반환.

    각 pairwise comparison에 대해:
    - 유의미한 유전자 수 (up / down / total)
    - 상위 유전자 목록 (|log2FC| 기준, padj 포함)
    - 상위 GO BP 경로 (up / down 각 최대 5개)
    - 경고 (유전자 0개, n이 적은 그룹 등)

    Args:
        de_config_file: DE 파이프라인 config YAML 경로
        pair:           특정 비교 pair만 조회 (예: "D1_vs_Control"). None이면 전체.
        top_n:          상위 유전자 반환 수 (기본 10)

    Returns:
        {
            "status": "success" | "partial" | "error",
            "pairs": {
                "D1_vs_Control": {
                    "sig_total": int,
                    "sig_up": int,
                    "sig_down": int,
                    "top_up_genes":   [{"gene": str, "log2FC": float, "padj": float}, ...],
                    "top_down_genes": [{"gene": str, "log2FC": float, "padj": float}, ...],
                    "top_go_up":   [{"term": str, "padj": float, "count": int}, ...],
                    "top_go_down": [{"term": str, "padj": float, "count": int}, ...],
                    "warnings": [str, ...]
                }, ...
            },
            "summary_text": str,   # LLM이 직접 활용할 수 있는 한국어 요약
            "not_found": [str, ...]  # 결과 파일이 없는 pair
        }
    """
    try:
        import yaml as _yaml
        import pandas as _pd

        de_cfg_path = Path(de_config_file)
        if not de_cfg_path.exists():
            return {"status": "error", "message": f"DE config not found: {de_config_file}"}

        with open(de_cfg_path) as f:
            de_cfg = _yaml.safe_load(f)

        de_pipeline_dir = de_cfg_path.parent.parent
        output_dir = de_pipeline_dir / de_cfg.get("output_dir", "output")
        padj_cut  = de_cfg.get("de_analysis", {}).get("padj_cutoff",  0.05)
        lfc_cut   = de_cfg.get("de_analysis", {}).get("log2fc_cutoff", 0.0)

        pairwise = de_cfg.get("de_analysis", {}).get("pairwise_comparisons", [])
        all_pairs = [f"{t}_vs_{b}" for t, b in pairwise]

        if pair:
            target_pairs = [pair]
        else:
            target_pairs = all_pairs

        results: Dict[str, Any] = {}
        not_found: List[str] = []

        for p in target_pairs:
            pair_dir  = output_dir / "pairwise" / p
            de_file   = pair_dir / "final_de_results.csv"

            if not de_file.exists():
                not_found.append(p)
                continue

            de_df = _pd.read_csv(de_file, index_col=0)

            # 유의미한 유전자 필터
            sig = de_df.dropna(subset=["padj", "log2FoldChange"])
            sig = sig[(sig["padj"] < padj_cut) & (sig["log2FoldChange"].abs() > lfc_cut)]
            up   = sig[sig["log2FoldChange"] > 0]
            down = sig[sig["log2FoldChange"] < 0]

            gene_col = "symbol" if "symbol" in de_df.columns else de_df.index.name or "gene"

            def _top_genes(df: "_pd.DataFrame", ascending: bool) -> List[Dict]:
                sub = df.sort_values("log2FoldChange", ascending=ascending).head(top_n)
                rows = []
                for _, row in sub.iterrows():
                    name = row.get("symbol", row.name) if "symbol" in sub.columns else row.name
                    rows.append({
                        "gene":    str(name),
                        "log2FC":  round(float(row["log2FoldChange"]), 3),
                        "padj":    float(f"{row['padj']:.2e}"),
                    })
                return rows

            top_up   = _top_genes(up,   ascending=False)
            top_down = _top_genes(down, ascending=True)

            # GO BP 상위 경로
            def _top_go(direction: str) -> List[Dict]:
                go_file = pair_dir / f"go_enrichment_{direction}_BP.csv"
                if not go_file.exists():
                    return []
                raw = Path(go_file).read_text().strip()
                if not raw or raw == '""':
                    return []
                try:
                    gdf = _pd.read_csv(go_file)
                    if gdf.empty:
                        return []
                    gdf = gdf.sort_values("p.adjust").head(5)
                    return [
                        {
                            "term":  row.get("Description", row.get("ID", "")),
                            "padj":  float(f"{row['p.adjust']:.2e}"),
                            "count": int(row.get("Count", 0)),
                        }
                        for _, row in gdf.iterrows()
                    ]
                except Exception:
                    return []

            top_go_up   = _top_go("up")
            top_go_down = _top_go("down")

            warnings: List[str] = []
            if len(sig) == 0:
                warnings.append("유의미한 DE 유전자 없음 — padj/log2FC 컷오프 재검토 권장")
            elif len(sig) < 10:
                warnings.append(f"유의미한 유전자가 {len(sig)}개로 매우 적음")

            results[p] = {
                "sig_total":      int(len(sig)),
                "sig_up":         int(len(up)),
                "sig_down":       int(len(down)),
                "top_up_genes":   top_up,
                "top_down_genes": top_down,
                "top_go_up":      top_go_up,
                "top_go_down":    top_go_down,
                "padj_cutoff":    padj_cut,
                "log2fc_cutoff":  lfc_cut,
                "warnings":       warnings,
            }

        # ── 한국어 요약 텍스트 생성 ─────────────────────────────────────────
        lines = [f"[DE 결과 요약]  (padj < {padj_cut}, |log2FC| > {lfc_cut})"]
        for p, r in results.items():
            lines.append(f"\n◆ {p}")
            lines.append(f"  유의미한 유전자: {r['sig_total']}개  (↑{r['sig_up']} / ↓{r['sig_down']})")
            if r["top_up_genes"]:
                top_names = ", ".join(g["gene"] for g in r["top_up_genes"][:5])
                lines.append(f"  상위 Up-regulated: {top_names}")
            if r["top_down_genes"]:
                top_names = ", ".join(g["gene"] for g in r["top_down_genes"][:5])
                lines.append(f"  상위 Down-regulated: {top_names}")
            if r["top_go_up"]:
                go_names = "; ".join(g["term"] for g in r["top_go_up"][:3])
                lines.append(f"  GO UP (BP): {go_names}")
            if r["top_go_down"]:
                go_names = "; ".join(g["term"] for g in r["top_go_down"][:3])
                lines.append(f"  GO DOWN (BP): {go_names}")
            for w in r["warnings"]:
                lines.append(f"  ⚠️  {w}")
        if not_found:
            lines.append(f"\n결과 없음 (아직 미실행): {', '.join(not_found)}")

        status = "error" if not results and not_found else (
            "partial" if not_found else "success"
        )
        return {
            "status":       status,
            "pairs":        results,
            "not_found":    not_found,
            "summary_text": "\n".join(lines),
        }

    except Exception as e:
        return {"status": "error", "message": f"read_de_results_summary failed: {e}"}


def check_analysis_readiness(de_config_file: str) -> Dict[str, Any]:
    """
    DE 분석 실행 전 종합 진단 — 실행 전에 발견할 수 있는 모든 문제를 한 번에 점검.

    점검 항목:
    1. 파일 존재 여부 (counts CSV, metadata CSV)
    2. counts ↔ metadata 샘플 ID 일치 여부
    3. pairwise_comparisons 조건명 ↔ metadata condition 일치 여부
    4. 조건별 샘플 수 (n < 3 경고)
    5. counts 기본 통계 (전체 유전자 수, 저발현 유전자 비율)

    Returns:
        {
            "status":          "ok" | "warning" | "error",
            "issues":          [{"level": "error"|"warning", "message": str}, ...],
            "sample_stats": {
                "n_samples_metadata": int,
                "n_samples_counts":   int,
                "samples_only_in_metadata": [str, ...],
                "samples_only_in_counts":   [str, ...],
            },
            "condition_stats": {condition: {"n": int, "samples": [str]}, ...},
            "counts_stats": {
                "n_genes": int,
                "n_zero_genes": int,
                "zero_pct": float,
            },
            "condition_check": {...},   # validate_de_config_conditions 결과
            "summary_text": str
        }
    """
    try:
        import yaml as _yaml
        import pandas as _pd

        de_cfg_path = Path(de_config_file)
        if not de_cfg_path.exists():
            return {"status": "error", "message": f"DE config not found: {de_config_file}"}

        with open(de_cfg_path) as f:
            de_cfg = _yaml.safe_load(f)

        de_pipeline_dir = de_cfg_path.parent.parent
        issues: List[Dict[str, str]] = []

        # ── 1. 파일 존재 여부 ────────────────────────────────────────────────
        def _resolve(raw: str) -> Path:
            p = Path(raw)
            return p if p.is_absolute() else de_pipeline_dir / p

        counts_path   = _resolve(de_cfg.get("count_data_path",  ""))
        metadata_path = _resolve(de_cfg.get("metadata_path",    ""))

        if not counts_path.exists():
            return {"status": "error",
                    "message": f"Counts file not found: {counts_path}",
                    "issues": [{"level": "error", "message": f"counts 파일 없음: {counts_path}"}]}
        if not metadata_path.exists():
            return {"status": "error",
                    "message": f"Metadata file not found: {metadata_path}",
                    "issues": [{"level": "error", "message": f"metadata 파일 없음: {metadata_path}"}]}

        # ── 2. 데이터 로드 ───────────────────────────────────────────────────
        counts_df = _pd.read_csv(counts_path, index_col=0)
        meta_df   = _pd.read_csv(metadata_path)

        counts_samples  = set(counts_df.columns.tolist())
        cond_col        = "condition" if "condition" in meta_df.columns else "group"
        meta_samples    = set(meta_df["sample_id"].tolist()) if "sample_id" in meta_df.columns else set()

        # ── 3. 샘플 ID 일치 ──────────────────────────────────────────────────
        only_meta   = sorted(meta_samples - counts_samples)
        only_counts = sorted(counts_samples - meta_samples)

        if only_meta:
            issues.append({"level": "error",
                           "message": f"metadata에만 있고 counts에 없는 샘플: {only_meta}"})
        if only_counts:
            issues.append({"level": "warning",
                           "message": f"counts에만 있고 metadata에 없는 샘플 (분석에서 제외됨): {only_counts}"})

        # ── 4. 조건명 일치 (validate_de_config_conditions 재사용) ────────────
        cond_check = validate_de_config_conditions(de_config_file)
        if cond_check["status"] == "mismatch":
            for m in cond_check["mismatches"]:
                suggestion = cond_check["suggestions"].get(m, "?")
                issues.append({"level": "error",
                               "message": f"조건명 불일치: config '{m}' → metadata에는 '{suggestion}'"})

        # ── 5. 조건별 샘플 수 ────────────────────────────────────────────────
        condition_stats: Dict[str, Any] = {}
        for cond, grp in meta_df.groupby(cond_col):
            sids = grp["sample_id"].tolist() if "sample_id" in grp.columns else []
            n    = len(grp)
            condition_stats[str(cond)] = {"n": n, "samples": sids}
            if n < 3:
                issues.append({"level": "warning",
                               "message": f"'{cond}' 그룹 샘플 수 n={n} (n≥3 권장, 통계적 신뢰도 낮음)"})
            if n < 2:
                issues.append({"level": "error",
                               "message": f"'{cond}' 그룹 샘플 수 n={n} — DESeq2 실행 불가"})

        # ── 6. counts 기본 통계 ──────────────────────────────────────────────
        n_genes      = len(counts_df)
        n_zero_genes = int((counts_df.sum(axis=1) == 0).sum())
        zero_pct     = round(n_zero_genes / n_genes * 100, 1) if n_genes else 0.0
        if zero_pct > 30:
            issues.append({"level": "warning",
                           "message": f"모든 샘플에서 발현량 0인 유전자 {zero_pct}% — 필터링 권장"})

        counts_stats = {
            "n_genes":      n_genes,
            "n_zero_genes": n_zero_genes,
            "zero_pct":     zero_pct,
        }
        sample_stats = {
            "n_samples_metadata":        len(meta_samples),
            "n_samples_counts":          len(counts_samples),
            "samples_only_in_metadata":  only_meta,
            "samples_only_in_counts":    only_counts,
        }

        # ── 요약 텍스트 ──────────────────────────────────────────────────────
        errors   = [i for i in issues if i["level"] == "error"]
        warnings = [i for i in issues if i["level"] == "warning"]
        status   = "error" if errors else ("warning" if warnings else "ok")

        lines = ["[DE 실행 전 진단]"]
        lines.append(f"  counts : {counts_path.name}  ({n_genes}개 유전자, {len(counts_samples)}개 샘플)")
        lines.append(f"  metadata: {metadata_path.name}  ({len(meta_samples)}개 샘플)")
        lines.append(f"\n조건별 샘플 수:")
        for cond, s in condition_stats.items():
            flag = " ⚠️" if s["n"] < 3 else ""
            lines.append(f"  {cond}: n={s['n']}{flag}")
        lines.append(f"\ncounts 통계: {n_genes}개 유전자, zero-count {zero_pct}%")

        if not issues:
            lines.append("\n✅ 이상 없음 — DE 분석 실행 가능")
        else:
            if errors:
                lines.append(f"\n🔴 오류 {len(errors)}개:")
                for i in errors:
                    lines.append(f"  - {i['message']}")
            if warnings:
                lines.append(f"\n🟡 경고 {len(warnings)}개:")
                for i in warnings:
                    lines.append(f"  - {i['message']}")

        return {
            "status":          status,
            "issues":          issues,
            "sample_stats":    sample_stats,
            "condition_stats": condition_stats,
            "counts_stats":    counts_stats,
            "condition_check": cond_check,
            "summary_text":    "\n".join(lines),
        }

    except Exception as e:
        return {"status": "error", "message": f"check_analysis_readiness failed: {e}"}


def setup_and_validate(
    config_file: str,
    sample_id: Optional[str] = None,
    n_reads: int = 500_000,
    cores: int = 8,
) -> Dict[str, Any]:
    """
    Combined pre-flight setup: detect FASTQs, create sample sheet, validate library type.

    Steps:
    1. Read config to get data_dir, project_id, sample_sheet path
    2. detect_fastq_files(data_dir)
    3. create_sample_sheet(project_id, data_dir, output_path=sample_sheet)
    4. validate_library_type(config_file, sample_id, n_reads, cores)

    Returns:
        {
            "status": "success"|"warning"|"error",
            "fastq_detection": {...},   # detect_fastq_files result
            "sample_sheet": {...},      # create_sample_sheet result
            "library_validation": {...},# validate_library_type result
            "summary": str,             # human-readable summary
        }
    """
    try:
        cfg_path = Path(config_file)
        if not cfg_path.exists():
            return {"status": "error", "message": f"Config not found: {config_file}"}

        import yaml as _yaml
        with open(cfg_path) as f:
            cfg = _yaml.safe_load(f)

        data_dir = cfg.get("data_dir", "")
        project_id = cfg.get("project_id", "")
        sample_sheet_path = cfg.get("sample_sheet", f"config/samples/{project_id}.tsv")

        # Resolve relative sample_sheet path relative to config file location
        if not Path(sample_sheet_path).is_absolute():
            pipeline_root = cfg_path.parent.parent.parent  # config/projects/ → root
            sample_sheet_path = str(pipeline_root / sample_sheet_path)

        # ── Step 1: Detect FASTQs ────────────────────────────────────────────
        fq_result = detect_fastq_files(data_dir)

        # ── Step 2: Create sample sheet ──────────────────────────────────────
        if fq_result.get("status") == "success":
            ss_result = create_sample_sheet(
                project_id=project_id,
                data_dir=data_dir,
                output_path=sample_sheet_path,
                overwrite=False,  # preserve manually revised conditions
            )
        else:
            ss_result = {"status": "error", "message": "Skipped — FASTQ detection failed"}

        # ── Step 3: Library validation ───────────────────────────────────────
        lib_result = validate_library_type(
            config_file=config_file,
            n_reads=n_reads,
            cores=cores,
            sample_id=sample_id,
        )

        # ── Build summary ────────────────────────────────────────────────────
        lines = []

        # FASTQ detection summary
        if fq_result.get("status") == "success":
            n_s = fq_result.get("n_samples", 0)
            n_f = fq_result.get("n_files", 0)
            lines.append(f"[FASTQ] {n_s}개 샘플, {n_f}개 파일 감지 완료")
        else:
            lines.append(f"[FASTQ] 오류: {fq_result.get('message', '알 수 없음')}")

        # Sample sheet summary
        if ss_result.get("status") == "success":
            n_ss = ss_result.get("n_samples", 0)
            conds = ss_result.get("conditions_assigned", {})
            cond_str = ", ".join(f"{k}: {v}개" for k, v in conds.items())
            skipped_note = " (기존 파일 유지 — 수동 수정 보존)" if ss_result.get("skipped") else ""
            lines.append(f"[샘플시트] {n_ss}개 샘플 등록 ({cond_str}){skipped_note}")
            unassigned = ss_result.get("unassigned", [])
            if unassigned:
                lines.append(f"  ⚠ 조건 미지정 샘플: {', '.join(unassigned)}")
        else:
            lines.append(f"[샘플시트] 오류: {ss_result.get('message', '알 수 없음')}")

        # Library validation summary
        if lib_result.get("status") in ("success", "warning"):
            strand_map = {0: "unstranded (0)", 1: "forward (1)", 2: "reverse (2)"}
            detected = lib_result.get("strandedness_detected", "?")
            verdict = lib_result.get("library_type_verdict", "?")
            exon_pct = lib_result.get("exon_pct", 0)
            updated = lib_result.get("strandedness_updated", False)
            before = lib_result.get("strandedness_before", "?")
            sample_used = lib_result.get("sample_used", "?")
            lines.append(
                f"[라이브러리 검증] 샘플: {sample_used}, "
                f"strandedness: {strand_map.get(detected, detected)}, "
                f"엑손 비율: {exon_pct:.1f}%, 판정: {verdict}"
            )
            if updated:
                lines.append(
                    f"  → strandedness config 자동 업데이트: {before} → {detected}"
                )
            for w in lib_result.get("warnings", []):
                lines.append(f"  ⚠ {w}")
        else:
            lines.append(f"[라이브러리 검증] 오류: {lib_result.get('message', '알 수 없음')}")

        # Overall status
        statuses = [fq_result.get("status"), ss_result.get("status"), lib_result.get("status")]
        if any(s == "error" for s in statuses):
            overall = "warning"  # partial failure — don't block everything
        elif any(s == "warning" for s in statuses):
            overall = "warning"
        else:
            overall = "success"

        return {
            "status": overall,
            "fastq_detection": fq_result,
            "sample_sheet": ss_result,
            "library_validation": lib_result,
            "summary": "\n".join(lines),
        }

    except Exception as e:
        return {"status": "error", "message": f"setup_and_validate failed: {e}"}


# ─────────────────────────────────────────────────────────────────────────────
# DE-GO Pipeline execution tools
# ─────────────────────────────────────────────────────────────────────────────

def run_de_pipeline(
    de_config_file: str,
    de_pipeline_dir: Optional[str] = None,
    cores: int = 8,
    dry_run: bool = True,
    background: bool = False,
) -> Dict[str, Any]:
    """
    Execute the DE-GO Snakemake pipeline.

    Mirrors run_pipeline() but targets the DE pipeline Snakefile.
    All rules use conda: "rna-seq-de-go-analysis" env, activated via --use-conda.

    Args:
        de_config_file: Path to DE config YAML (e.g. configs/config_*.yml)
        de_pipeline_dir: Root of DE pipeline (contains Snakefile). If None,
                         walks up from de_config_file to find Snakefile.
        cores: Number of CPU cores
        dry_run: If True (default), show jobs without running (always blocking)
        background: If True and not dry_run, launch in background with Popen

    Returns:
        dry_run:    {"status": "dry_run_ok"|"dry_run_error", "total_jobs": int, ...}
        background: {"status": "running", "pid": int, "log_file": str, ...}
        blocking:   {"status": "success"|"error", "returncode": int, ...}
    """
    try:
        de_cfg_path = Path(de_config_file)
        if not de_cfg_path.exists():
            return {"status": "error", "message": f"DE config not found: {de_config_file}"}

        # Resolve DE pipeline root (directory containing Snakefile)
        if de_pipeline_dir:
            pipeline_root = Path(de_pipeline_dir).resolve()
        else:
            # Walk up from config file to find Snakefile
            pipeline_root = de_cfg_path.resolve().parent
            for _ in range(5):
                if (pipeline_root / "Snakefile").exists():
                    break
                pipeline_root = pipeline_root.parent
            else:
                return {
                    "status": "error",
                    "message": f"Could not find Snakefile from {de_config_file}. Provide de_pipeline_dir explicitly."
                }

        snakefile = pipeline_root / "Snakefile"
        if not snakefile.exists():
            return {"status": "error", "message": f"Snakefile not found at {snakefile}"}

        # 4-tier snakemake path resolution (same as run_pipeline)
        import shutil as _shutil
        _snakemake_candidates = [
            Path(sys.executable).parent / "snakemake",
            Path(os.environ.get("CONDA_PREFIX", "/nonexistent")) / "bin" / "snakemake",
            Path("/home/ygkim/program/anaconda3/envs/rna-seq-pipeline/bin/snakemake"),
        ]
        snakemake_cmd = next(
            (str(p) for p in _snakemake_candidates if p.exists()),
            _shutil.which("snakemake") or "snakemake",
        )

        cmd = [
            snakemake_cmd,
            "--snakefile", str(snakefile),
            "--configfile", str(de_cfg_path.resolve()),
            "--use-conda",
            "--cores", str(cores),
            "--directory", str(pipeline_root),
        ]

        if dry_run:
            cmd.append("--dry-run")

        # ── Background launch ───────────────────────────────────────────────
        if background and not dry_run:
            log_dir = pipeline_root / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / f"snakemake_{de_cfg_path.stem}.log"
            pid_file = log_dir / f"snakemake_{de_cfg_path.stem}.pid"

            with open(log_file, "w") as lf:
                proc = subprocess.Popen(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    cwd=str(pipeline_root),
                )
            pid_file.write_text(str(proc.pid))

            return {
                "status": "running",
                "pid": proc.pid,
                "log_file": str(log_file),
                "pid_file": str(pid_file),
                "command": " ".join(cmd),
                "de_config_file": str(de_cfg_path),
                "de_pipeline_dir": str(pipeline_root),
                "note": f"DE pipeline launched in background. Use monitor_de_pipeline() to check progress, or: tail -f {log_file}",
            }

        # ── Blocking (dry-run or background=False) ──────────────────────────
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(pipeline_root))
        stdout = result.stdout
        stderr = result.stderr
        combined = stdout + stderr

        if dry_run:
            jobs_by_rule: Dict[str, int] = {}

            stats_match = re.search(r'Job stats:\s*\n', combined)
            if stats_match:
                block = combined[stats_match.end():]
                for line in block.splitlines():
                    if not line.strip() or line.startswith('=') or line.startswith('Building'):
                        break
                    if re.match(r'^[-\s]+$', line) or line.lstrip().startswith('job') or line.lstrip().startswith('rule'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            jobs_by_rule[parts[0]] = int(parts[-1])
                        except ValueError:
                            pass

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
            issues = [
                ln.strip() for ln in combined.splitlines()
                if any(kw in ln.lower() for kw in ('error:', 'exception', 'traceback', 'syntaxerror'))
                and 'missing' not in ln.lower()
            ][:20]

            status = "dry_run_ok" if result.returncode == 0 else "dry_run_error"
            result_dict: Dict[str, Any] = {
                "status": status,
                "returncode": result.returncode,
                "dry_run": True,
                "total_jobs": total_jobs,
                "jobs_by_rule": jobs_by_rule,
                "issues": issues,
                "de_config_file": str(de_cfg_path),
                "de_pipeline_dir": str(pipeline_root),
                "command": " ".join(cmd),
                "note": (
                    "dry_run_ok: Jobs listed above will run. Call run_de_pipeline(dry_run=False, background=True) to execute."
                    if status == "dry_run_ok" else
                    "dry_run_error: Configuration error. Check issues[] and stderr_snippet."
                ),
            }
            if result.returncode != 0 or not jobs_by_rule:
                result_dict["stderr_snippet"] = (stderr or stdout)[:3000]
            return result_dict

        else:
            output_lines = combined.splitlines()
            return {
                "status": "success" if result.returncode == 0 else "error",
                "returncode": result.returncode,
                "dry_run": False,
                "output_tail": "\n".join(output_lines[-50:]),
                "command": " ".join(cmd),
                "de_config_file": str(de_cfg_path),
                "de_pipeline_dir": str(pipeline_root),
            }

    except Exception as e:
        return {"status": "error", "message": f"run_de_pipeline failed: {e}"}


def monitor_de_pipeline(
    de_config_file: Optional[str] = None,
    de_pipeline_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Check DE-GO pipeline execution status by parsing Snakemake logs.

    Simpler than monitor_pipeline() — no STAR logs to parse.

    Returns:
        {
            "status": "not_started"|"running"|"completed"|"failed",
            "progress_pct": float | None,
            "jobs_done": int | None,
            "jobs_total": int | None,
            "current_rule": str | None,
            "error_rules": [...],
            "log_file": str | None,
            "pid": int | None,
            "last_lines": str,
        }
    """
    try:
        # Resolve pipeline root
        pipeline_root: Optional[Path] = None
        if de_pipeline_dir:
            pipeline_root = Path(de_pipeline_dir).resolve()
        elif de_config_file:
            p = Path(de_config_file).resolve().parent
            for _ in range(5):
                if (p / "Snakefile").exists():
                    pipeline_root = p
                    break
                p = p.parent

        # Find log file
        log_path = _find_snakemake_log(de_config_file, pipeline_root)

        if log_path is None:
            return {
                "status": "not_started",
                "progress_pct": None,
                "jobs_done": None,
                "jobs_total": None,
                "current_rule": None,
                "error_rules": [],
                "log_file": None,
                "pid": None,
                "last_lines": "",
                "message": "No log file found. Has run_de_pipeline() been called yet?",
            }

        log_info = _parse_snakemake_log(log_path)

        # Check PID file to determine if process is still alive
        pid: Optional[int] = None
        is_running = False
        if de_config_file and pipeline_root:
            pid_file = pipeline_root / "logs" / f"snakemake_{Path(de_config_file).stem}.pid"
            if pid_file.exists():
                try:
                    pid = int(pid_file.read_text().strip())
                    os.kill(pid, 0)   # signal 0: just probe, no actual signal
                    is_running = True
                except (ValueError, ProcessLookupError, PermissionError):
                    is_running = False

        # Determine status
        if log_info.get("fatal_error") or log_info.get("error_rules"):
            status = "failed"
        elif log_info.get("is_complete"):
            status = "completed"
        elif is_running:
            status = "running"
        else:
            # Process gone but no completion signal → likely failed or interrupted
            status = "failed" if not log_info.get("is_complete") else "completed"

        # Last 20 lines of log for context
        try:
            log_lines = log_path.read_text(errors='ignore').splitlines()
            last_lines = "\n".join(log_lines[-20:])
        except OSError:
            last_lines = ""

        return {
            "status": status,
            "progress_pct": log_info.get("progress_pct"),
            "jobs_done": log_info.get("jobs_done"),
            "jobs_total": log_info.get("jobs_total"),
            "current_rule": log_info.get("current_rule"),
            "error_rules": log_info.get("error_rules", []),
            "log_file": str(log_path),
            "pid": pid,
            "last_lines": last_lines,
        }

    except Exception as e:
        return {"status": "error", "message": f"monitor_de_pipeline failed: {e}"}
