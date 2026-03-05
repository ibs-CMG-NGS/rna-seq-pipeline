#!/usr/bin/env python3
"""
Batch Runner for RNA-seq Pipeline

Runs multiple projects sequentially without LLM involvement.
Each project goes through: validate → run_pipeline → qc_check → bridge (optional)

Checkpoint files (batch_logs/<project_id>/<step>.done) allow resuming interrupted runs.

Usage:
    python scripts/batch_runner.py config/batch/my_batch.yaml [--dry-run] [--resume] [--cores 12]
"""

import argparse
import csv
import logging
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import yaml

# Add project root to sys.path so pipeline_tools can be imported
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.pipeline_tools import (
    validate_input_data,
    run_pipeline,
    read_qc_results,
    run_bridge,
)

STEPS = ["validate", "pipeline", "qc_check", "bridge"]


# ─── Checkpoint helpers ───────────────────────────────────────────────────────

def checkpoint_path(log_dir: Path, step: str) -> Path:
    return log_dir / f"{step}.done"


def checkpoint_exists(log_dir: Path, step: str) -> bool:
    return checkpoint_path(log_dir, step).exists()


def write_checkpoint(log_dir: Path, step: str) -> None:
    checkpoint_path(log_dir, step).write_text(datetime.now().isoformat())


# ─── Per-project runner ───────────────────────────────────────────────────────

def run_project(project_cfg: dict, global_cfg: dict, dry_run: bool, resume: bool,
                project_logger: logging.Logger) -> dict:
    """
    Execute all steps for one project. Returns a result dict with status/metrics.
    """
    project_id = project_cfg["project_id"]
    config_file = project_cfg["config_file"]
    cores = global_cfg.get("cores", 8)
    de_pipeline_dir = global_cfg.get("de_pipeline_dir")

    log_dir = Path("batch_logs") / project_id
    log_dir.mkdir(parents=True, exist_ok=True)

    result = {
        "project_id": project_id,
        "status": "success",
        "failed_step": "",
        "duration_sec": 0,
        "qc_pass": "",
        "mapping_pct": "",
        "error": "",
    }
    start = time.time()

    # ── Step 1: validate ──────────────────────────────────────────────────────
    if checkpoint_exists(log_dir, "validate") and resume:
        project_logger.info("[validate] skipped (checkpoint exists)")
    else:
        project_logger.info("[validate] starting...")
        try:
            r = validate_input_data(config_file)
        except Exception as e:
            project_logger.error(f"[validate] exception: {e}\n{traceback.format_exc()}")
            result.update(status="failed", failed_step="validate", error=str(e))
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        if r.get("status") == "invalid":
            msg = r.get("errors") or r.get("message") or str(r)
            project_logger.error(f"[validate] failed: {msg}")
            result.update(status="failed", failed_step="validate",
                          error=str(msg)[:200])
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        write_checkpoint(log_dir, "validate")
        project_logger.info("[validate] OK")

    # ── Step 2: pipeline ──────────────────────────────────────────────────────
    if checkpoint_exists(log_dir, "pipeline") and resume:
        project_logger.info("[pipeline] skipped (checkpoint exists)")
    else:
        project_logger.info(f"[pipeline] starting (cores={cores}, dry_run={dry_run})...")
        try:
            r = run_pipeline(config_file, cores=cores, dry_run=dry_run, background=False)
        except Exception as e:
            project_logger.error(f"[pipeline] exception: {e}\n{traceback.format_exc()}")
            result.update(status="failed", failed_step="pipeline", error=str(e))
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        if r.get("status") == "error":
            msg = r.get("message", str(r))
            project_logger.error(f"[pipeline] failed: {msg}")
            result.update(status="failed", failed_step="pipeline", error=str(msg))
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        write_checkpoint(log_dir, "pipeline")
        project_logger.info(f"[pipeline] OK (status={r.get('status')})")

    # ── Step 3: QC check ─────────────────────────────────────────────────────
    if checkpoint_exists(log_dir, "qc_check") and resume:
        project_logger.info("[qc_check] skipped (checkpoint exists)")
    else:
        project_logger.info("[qc_check] reading QC results...")
        try:
            r = read_qc_results(config_file)
        except Exception as e:
            project_logger.warning(f"[qc_check] exception (non-fatal): {e}")
            r = {}

        # Extract key metrics for summary
        fastqc = r.get("fastqc_evaluation") or {}
        star = r.get("star_alignment") or {}
        summary = r.get("summary") or {}

        overall_pass = summary.get("overall_pass")
        avg_mapping = star.get("average_uniquely_mapped_pct")

        result["qc_pass"] = str(overall_pass) if overall_pass is not None else ""
        result["mapping_pct"] = f"{avg_mapping:.1f}%" if avg_mapping is not None else ""

        project_logger.info(
            f"[qc_check] overall_pass={overall_pass}, avg_mapping={result['mapping_pct']}"
        )

        write_checkpoint(log_dir, "qc_check")

    # ── Step 4: bridge (optional) ─────────────────────────────────────────────
    if not project_cfg.get("run_bridge", False):
        project_logger.info("[bridge] skipped (run_bridge=false)")
    elif not de_pipeline_dir:
        project_logger.warning("[bridge] skipped (de_pipeline_dir not set in batch config)")
    elif checkpoint_exists(log_dir, "bridge") and resume:
        project_logger.info("[bridge] skipped (checkpoint exists)")
    else:
        project_logger.info(f"[bridge] running → {de_pipeline_dir}")
        try:
            r = run_bridge(config_file, de_pipeline_dir, skip_de=True)
        except Exception as e:
            project_logger.error(f"[bridge] exception: {e}\n{traceback.format_exc()}")
            result.update(status="failed", failed_step="bridge", error=str(e))
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        if r.get("status") == "error":
            msg = r.get("message", str(r))
            project_logger.error(f"[bridge] failed: {msg}")
            result.update(status="failed", failed_step="bridge", error=str(msg))
            result["duration_sec"] = round(time.time() - start, 1)
            return result

        write_checkpoint(log_dir, "bridge")
        project_logger.info("[bridge] OK")

    result["duration_sec"] = round(time.time() - start, 1)
    return result


# ─── Summary CSV ─────────────────────────────────────────────────────────────

SUMMARY_FIELDS = ["project_id", "status", "failed_step", "duration_sec",
                  "qc_pass", "mapping_pct", "error"]


def write_summary(results: list, summary_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, "") for k in SUMMARY_FIELDS})


# ─── Logging setup ───────────────────────────────────────────────────────────

def setup_logger(name: str, log_file: Path, level=logging.INFO) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    logger = logging.getLogger(name)
    if logger.handlers:
        return logger  # already configured
    logger.setLevel(level)
    logger.propagate = False

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq batch runner — sequential multi-project execution"
    )
    parser.add_argument("batch_config", help="Path to batch YAML config")
    parser.add_argument("--dry-run", action="store_true",
                        help="Pass dry_run=True to run_pipeline (Snakemake --dry-run)")
    parser.add_argument("--resume", action="store_true",
                        help="Skip steps that already have a .done checkpoint")
    parser.add_argument("--cores", type=int, default=None,
                        help="Override cores from batch config")
    args = parser.parse_args()

    batch_cfg_path = Path(args.batch_config)
    if not batch_cfg_path.exists():
        print(f"ERROR: batch config not found: {batch_cfg_path}", file=sys.stderr)
        sys.exit(1)

    with open(batch_cfg_path) as f:
        batch_cfg = yaml.safe_load(f)

    global_cfg = {
        "cores": args.cores or batch_cfg.get("cores", 8),
        "de_pipeline_dir": batch_cfg.get("de_pipeline_dir"),
        "pipeline_root": batch_cfg.get("pipeline_root"),
    }

    projects = batch_cfg.get("projects", [])
    if not projects:
        print("ERROR: no projects defined in batch config", file=sys.stderr)
        sys.exit(1)

    # Root logger
    root_log_dir = Path("batch_logs")
    root_logger = setup_logger(
        "batch", root_log_dir / "batch.log"
    )

    root_logger.info("=" * 60)
    root_logger.info(f"Batch run started: {len(projects)} projects")
    root_logger.info(f"  cores={global_cfg['cores']}, dry_run={args.dry_run}, resume={args.resume}")
    root_logger.info("=" * 60)

    all_results = []

    for i, project_cfg in enumerate(projects, 1):
        pid = project_cfg.get("project_id", f"project_{i}")
        root_logger.info(f"\n[{i}/{len(projects)}] {pid}")

        project_log_file = root_log_dir / pid / "runner.log"
        project_logger = setup_logger(f"batch.{pid}", project_log_file)

        try:
            result = run_project(
                project_cfg, global_cfg,
                dry_run=args.dry_run,
                resume=args.resume,
                project_logger=project_logger,
            )
        except Exception as e:
            root_logger.error(f"Unexpected error for {pid}: {e}\n{traceback.format_exc()}")
            result = {
                "project_id": pid,
                "status": "failed",
                "failed_step": "unknown",
                "duration_sec": 0,
                "qc_pass": "",
                "mapping_pct": "",
                "error": str(e),
            }

        all_results.append(result)
        status_icon = "OK" if result["status"] == "success" else "FAIL"
        root_logger.info(
            f"  [{status_icon}] {pid} — {result['duration_sec']}s"
            + (f" — failed at: {result.get('failed_step')}" if result["status"] != "success" else "")
        )

    # Write summary CSV
    summary_path = root_log_dir / "summary.csv"
    write_summary(all_results, summary_path)

    n_ok = sum(1 for r in all_results if r["status"] == "success")
    n_fail = len(all_results) - n_ok

    root_logger.info("\n" + "=" * 60)
    root_logger.info(f"Batch complete: {n_ok} succeeded, {n_fail} failed")
    root_logger.info(f"Summary: {summary_path}")
    root_logger.info("=" * 60)

    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
