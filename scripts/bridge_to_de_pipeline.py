#!/usr/bin/env python3
"""
Bridge: RNA-seq Pipeline → DE/GO Analysis Pipeline

Automatically transfer count matrix to DE/GO pipeline and trigger analysis.
"""

import json
import shutil
import argparse
from pathlib import Path
import subprocess
import sys

class PipelineBridge:
    """Connect RNA-seq preprocessing to DE/GO analysis."""
    
    def __init__(self, 
                 rnaseq_output_dir: Path,
                 de_pipeline_dir: Path,
                 project_id: str):
        self.rnaseq_output = Path(rnaseq_output_dir)
        self.de_pipeline = Path(de_pipeline_dir)
        self.project_id = project_id
        
    def check_rnaseq_completion(self) -> bool:
        """Check if RNA-seq pipeline completed successfully."""
        # Check for counts matrix
        counts_file = self.rnaseq_output / "counts" / "counts_matrix.txt"
        if not counts_file.exists():
            print(f"❌ Counts matrix not found: {counts_file}")
            return False
        
        # Check project summary
        summary_file = self.rnaseq_output / "project_summary.json"
        if not summary_file.exists():
            print(f"⚠️  Project summary not found: {summary_file}")
            return True  # Not critical, continue
        
        # Check QC status
        with open(summary_file) as f:
            summary = json.load(f)
        
        qc_status = summary['qc_summary']
        print(f"\n📊 QC Summary:")
        print(f"  Total samples: {qc_status['total_samples']}")
        print(f"  Passed: {qc_status['passed']}")
        print(f"  Failed: {qc_status['failed']}")
        print(f"  Pass rate: {qc_status['pass_rate']}%")
        
        if qc_status['failed'] > 0:
            print(f"\n⚠️  Warning: {qc_status['failed']} samples failed QC")
            response = input("Continue with DE analysis? (y/n): ")
            if response.lower() != 'y':
                return False
        
        return True
    
    def prepare_de_input(self) -> Path:
        """Copy counts matrix to DE pipeline input directory."""
        # Source files
        counts_file = self.rnaseq_output / "counts" / "counts_matrix.txt"
        
        # Destination
        de_data_dir = self.de_pipeline / "data" / "raw"
        de_data_dir.mkdir(parents=True, exist_ok=True)
        
        dest_counts = de_data_dir / f"{self.project_id}_counts.csv"
        
        # Copy counts
        print(f"\n📋 Copying counts matrix:")
        print(f"  From: {counts_file}")
        print(f"  To: {dest_counts}")
        shutil.copy(counts_file, dest_counts)
        
        return dest_counts
    
    def generate_metadata_template(self) -> Path:
        """Generate metadata template from sample sheet."""
        # Read sample sheet from RNA-seq pipeline
        sample_sheet = self.rnaseq_output.parent / "config" / "samples" / f"{self.project_id}.tsv"
        
        if not sample_sheet.exists():
            print(f"\n⚠️  Sample sheet not found: {sample_sheet}")
            print("You'll need to create metadata manually for DE analysis")
            return None
        
        # TODO: Parse sample sheet and generate metadata template
        # For now, just inform user
        print(f"\n📝 Sample sheet found: {sample_sheet}")
        print("Please create metadata file for DE analysis based on this.")
        
        return sample_sheet
    
    def trigger_de_analysis(self, dry_run: bool = False):
        """Trigger DE/GO analysis pipeline."""
        print(f"\n🚀 Triggering DE/GO analysis...")
        
        cmd = [
            "snakemake",
            "-s", str(self.de_pipeline / "Snakefile"),
            "--cores", "4",
            "--use-conda"
        ]
        
        if dry_run:
            cmd.extend(["--dry-run", "--printshellcmds"])
            print("DRY RUN mode - no actual execution")
        
        print(f"Command: {' '.join(cmd)}")
        
        if not dry_run:
            response = input("\nExecute DE/GO pipeline? (y/n): ")
            if response.lower() != 'y':
                print("Aborted.")
                return
        
        try:
            subprocess.run(cmd, cwd=self.de_pipeline, check=True)
            print("\n✅ DE/GO analysis completed!")
        except subprocess.CalledProcessError as e:
            print(f"\n❌ DE/GO analysis failed: {e}")
            sys.exit(1)
    
    def run(self, dry_run: bool = False, skip_de: bool = False):
        """Execute full bridge workflow."""
        print(f"{'='*60}")
        print(f"RNA-seq → DE/GO Pipeline Bridge")
        print(f"Project: {self.project_id}")
        print(f"{'='*60}")
        
        # Step 1: Check RNA-seq completion
        print("\n[Step 1] Checking RNA-seq pipeline completion...")
        if not self.check_rnaseq_completion():
            print("❌ RNA-seq pipeline not ready. Aborting.")
            sys.exit(1)
        
        print("✅ RNA-seq pipeline completed")
        
        # Step 2: Prepare DE input
        print("\n[Step 2] Preparing DE/GO pipeline input...")
        counts_file = self.prepare_de_input()
        print(f"✅ Counts matrix prepared: {counts_file}")
        
        # Step 3: Generate metadata template
        print("\n[Step 3] Checking sample metadata...")
        self.generate_metadata_template()
        
        # Step 4: Trigger DE analysis (optional)
        if not skip_de:
            self.trigger_de_analysis(dry_run=dry_run)
        else:
            print("\n⏭️  Skipping DE analysis (--skip-de flag)")
            print(f"\nTo run DE analysis manually:")
            print(f"  cd {self.de_pipeline}")
            print(f"  snakemake --cores 4 --use-conda")

def main():
    parser = argparse.ArgumentParser(
        description="Bridge RNA-seq preprocessing to DE/GO analysis"
    )
    parser.add_argument(
        '--rnaseq-output',
        type=Path,
        required=True,
        help='RNA-seq pipeline output directory'
    )
    parser.add_argument(
        '--de-pipeline',
        type=Path,
        required=True,
        help='DE/GO pipeline directory'
    )
    parser.add_argument(
        '--project-id',
        required=True,
        help='Project identifier'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing'
    )
    parser.add_argument(
        '--skip-de',
        action='store_true',
        help='Only prepare files, do not trigger DE analysis'
    )
    
    args = parser.parse_args()
    
    # Validate paths
    if not args.rnaseq_output.exists():
        print(f"ERROR: RNA-seq output not found: {args.rnaseq_output}")
        sys.exit(1)
    
    if not args.de_pipeline.exists():
        print(f"ERROR: DE/GO pipeline not found: {args.de_pipeline}")
        sys.exit(1)
    
    # Execute bridge
    bridge = PipelineBridge(
        args.rnaseq_output,
        args.de_pipeline,
        args.project_id
    )
    
    bridge.run(dry_run=args.dry_run, skip_de=args.skip_de)

if __name__ == '__main__':
    main()
