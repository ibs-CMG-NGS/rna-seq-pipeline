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
import os
import pandas as pd
import yaml
from datetime import datetime

# Import auto-config utilities
try:
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from scripts.utils.auto_config import ensure_bridge_config
    AUTO_CONFIG_AVAILABLE = True
except ImportError:
    AUTO_CONFIG_AVAILABLE = False
    print("⚠️  Warning: auto_config not available, --config will be required")


def load_path_config(config_path: Path) -> dict:
    """Load path configuration from YAML file."""
    if not config_path or not config_path.exists():
        return {}
    
    with open(config_path) as f:
        config = yaml.safe_load(f) or {}
    
    return config

class PipelineBridge:
    """Connect RNA-seq preprocessing to DE/GO analysis."""
    
    def __init__(self, 
                 rnaseq_output_dir: Path,
                 de_pipeline_dir: Path,
                 project_id: str,
                 config: dict = None,
                 assume_yes: bool = False):
        self.rnaseq_output = Path(rnaseq_output_dir)
        self.de_pipeline = Path(de_pipeline_dir)
        self.project_id = project_id
        self.config = config or {}
        self.assume_yes = assume_yes
        
    def check_rnaseq_completion(self) -> bool:
        """Check if RNA-seq pipeline completed successfully."""
        # Check for counts matrix
        counts_relpath = self.config.get('counts_relpath', 'project_summary/counts/counts_matrix_clean.csv')
        counts_file = self.rnaseq_output / counts_relpath
        if not counts_file.exists():
            print(f"❌ Counts matrix not found: {counts_file}")
            return False
        
        # Check project summary
        summary_relpath = self.config.get('summary_relpath', 'project_summary.json')
        summary_file = self.rnaseq_output / summary_relpath
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
            if self.assume_yes:
                print("Continuing with failed samples (--yes flag)")
            else:
                response = input("Continue with DE analysis? (y/n): ")
                if response.lower() != 'y':
                    return False
        
        return True
    
    def prepare_de_input(self) -> Path:
        """Copy counts matrix to DE pipeline input directory."""
        # Source files
        counts_relpath = self.config.get('counts_relpath', 'project_summary/counts/counts_matrix_clean.csv')
        counts_file = self.rnaseq_output / counts_relpath
        
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
        """Generate metadata CSV from sample sheet for DE analysis."""
        # Get sample sheet directory from config or use defaults
        sample_sheet_dir = self.config.get('sample_sheet_dir')
        
        # Convert project_id: try both with hyphens and underscores
        # (e.g., "mouse-chd8" -> try both "mouse-chd8.tsv" and "mouse_chd8.tsv")
        project_ids = [self.project_id, self.project_id.replace('-', '_')]
        
        # Try multiple possible locations for sample sheet
        possible_paths = []
        
        # First priority: config specified directory
        if sample_sheet_dir:
            for pid in project_ids:
                possible_paths.append(Path(sample_sheet_dir) / f"{pid}.tsv")
        
        # Fallback locations
        for pid in project_ids:
            possible_paths.extend([
                Path("/data_3tb/shared/rna-seq-pipeline/config/samples") / f"{pid}.tsv",
                self.rnaseq_output.parent / "config" / "samples" / f"{pid}.tsv",
                Path("/home/ygkim/ngs-pipeline/rna-seq-pipeline/config/samples") / f"{pid}.tsv",
            ])
        
        sample_sheet_path = None
        for path in possible_paths:
            if path.exists():
                sample_sheet_path = path
                break
        
        if not sample_sheet_path:
            print(f"\n⚠️  Sample sheet not found. Tried:")
            for path in possible_paths:
                print(f"    - {path}")
            print("You'll need to create metadata manually for DE analysis")
            return None
        
        print(f"\n📝 Generating metadata from sample sheet...")
        print(f"  Source: {sample_sheet_path}")
        
        # Read sample sheet
        df = pd.read_csv(sample_sheet_path, sep='\t')
        
        # Create metadata for DE analysis (sample_id, condition, replicate)
        metadata_cols = ['sample_id', 'condition', 'replicate']
        
        # Check if required columns exist
        missing_cols = [col for col in metadata_cols if col not in df.columns]
        if missing_cols:
            print(f"⚠️  Missing columns in sample sheet: {missing_cols}")
            print("Using available columns...")
            metadata_cols = [col for col in metadata_cols if col in df.columns]
        
        # Add optional columns if available
        optional_cols = ['tissue', 'sex', 'age', 'treatment', 'time_point']
        for col in optional_cols:
            if col in df.columns:
                metadata_cols.append(col)
        
        # Extract metadata
        metadata_df = df[metadata_cols].copy()
        
        # Save metadata
        de_data_dir = self.de_pipeline / "data" / "raw"
        metadata_path = de_data_dir / f"{self.project_id}_metadata.csv"
        metadata_df.to_csv(metadata_path, index=False)
        
        print(f"  Output: {metadata_path}")
        print(f"  Columns: {', '.join(metadata_cols)}")
        print(f"  Samples: {len(metadata_df)}")
        print(f"  Conditions: {metadata_df['condition'].unique().tolist()}")
        
        return metadata_path
    
    def generate_de_config(self, counts_path: Path, metadata_path: Path) -> Path:
        """Generate config.yml for DE/GO analysis pipeline from template."""
        print(f"\n⚙️  Generating DE/GO pipeline config...")
        
        # Load template config
        template_path = self.de_pipeline / "configs" / "template" / "config.yml"
        if not template_path.exists():
            print(f"⚠️  Template config not found: {template_path}")
            print("Creating config from scratch...")
            config = {}
        else:
            print(f"  Loading template: {template_path}")
            with open(template_path) as f:
                config = yaml.safe_load(f)
        
        # Read metadata to determine species and conditions
        metadata_df = pd.read_csv(metadata_path)
        # Detect condition column (prefer 'condition', fall back to 'group')
        condition_col = 'condition' if 'condition' in metadata_df.columns else 'group'
        conditions = metadata_df[condition_col].unique().tolist()

        # Auto-detect species from project_id or use default
        species = "mouse" if "mouse" in self.project_id.lower() or "chd8" in self.project_id.lower() else "human"

        # Determine pairwise comparisons (all vs first condition as control)
        # Format: [treatment, control] — matches template and Snakefile parsing
        control_condition = conditions[0]  # Assume first is control
        pairwise_comparisons = []
        for cond in conditions[1:]:
            pairwise_comparisons.append([cond, control_condition])

        # Update config with project-specific values
        config['species'] = species
        config['gene_id_type'] = 'ENSEMBL'  # From featureCounts
        config['count_data_path'] = f"data/raw/{counts_path.name}"
        config['metadata_path'] = f"data/raw/{metadata_path.name}"
        config['output_dir'] = f"output/{self.project_id}"

        # Update DE analysis section
        if 'de_analysis' not in config:
            config['de_analysis'] = {}

        config['de_analysis']['pairwise_comparisons'] = pairwise_comparisons
        # Set group_variable to match actual metadata column name
        config['de_analysis']['group_variable'] = condition_col
        config['de_analysis']['design_formula'] = f"~ {condition_col}"
        
        # Add generation metadata as comments (YAML doesn't support inline comments in dicts)
        config_header = f"""# Generated by RNA-seq pipeline bridge
# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# Project: {self.project_id}
# Auto-detected species: {species}
# 
# ⚠️  IMPORTANT: Please review and adjust if needed!
#    - Verify species is correct
#    - Check control condition (currently: {control_condition})
#    - Review pairwise comparisons
#    - Adjust other parameters as needed
#
"""
        
        # Save config
        config_path = self.de_pipeline / "configs" / f"config_{self.project_id}.yml"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(config_path, 'w') as f:
            f.write(config_header)
            yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
        
        print(f"  Output: {config_path}")
        print(f"  Species: {species}")
        print(f"  Conditions: {', '.join(conditions)}")
        print(f"  Control: {control_condition} (assumed)")
        print(f"  Comparisons: {len(pairwise_comparisons)}")
        for comp in pairwise_comparisons:
            print(f"    - {comp[0]}_vs_{comp[1]}")
        print(f"\n  ⚠️  Please review config file before running DE analysis!")
        print(f"     Config: {config_path}")
        
        return config_path
    
    def trigger_de_analysis(self, dry_run: bool = False, config_file: Path = None):
        """Trigger DE/GO analysis pipeline."""
        print(f"\n🚀 Triggering DE/GO analysis...")
        
        if config_file:
            print(f"  Using config: {config_file.name}")
        
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
        
        if config_file:
            print(f"\n⚠️  Note: You may need to update Snakefile line 6:")
            print(f'     CONFIG_FILE = "configs/config_{self.project_id}.yml"')
        
        if not dry_run:
            if self.assume_yes:
                print("\nExecuting DE/GO pipeline (--yes flag)")
            else:
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
        
        # Step 2: Prepare DE input - counts matrix
        print("\n[Step 2] Preparing counts matrix...")
        counts_file = self.prepare_de_input()
        print(f"✅ Counts matrix prepared: {counts_file}")
        
        # Step 3: Generate metadata from sample sheet
        print("\n[Step 3] Generating metadata...")
        metadata_file = self.generate_metadata_template()
        if metadata_file:
            print(f"✅ Metadata generated: {metadata_file}")
        else:
            print("⚠️  Metadata generation failed. Manual creation needed.")
            if not dry_run and not skip_de:
                if self.assume_yes:
                    print("   (Continuing without metadata due to --yes flag)")
                else:
                    response = input("Continue without metadata? (y/n): ")
                    if response.lower() != 'y':
                        sys.exit(1)
            elif skip_de:
                print("   (Continuing in --skip-de mode)")
        
        # Step 4: Generate DE/GO config file
        print("\n[Step 4] Generating DE/GO pipeline config...")
        if metadata_file:
            config_file = self.generate_de_config(counts_file, metadata_file)
            print(f"✅ Config generated: {config_file}")
        else:
            print("⚠️  Config generation skipped (no metadata)")
            config_file = None
        
        # Step 5: Trigger DE analysis (optional)
        if not skip_de and config_file:
            self.trigger_de_analysis(dry_run=dry_run, config_file=config_file)
        else:
            if skip_de:
                print("\n⏭️  Skipping DE analysis (--skip-de flag)")
            else:
                print("\n⏭️  Skipping DE analysis (no config file)")
            
            print(f"\n📋 To run DE analysis manually:")
            print(f"  cd {self.de_pipeline}")
            if config_file:
                print(f"  # Edit Snakefile line 6 to use: config_{self.project_id}.yml")
            print(f"  snakemake --cores 4 --use-conda")

def main():
    parser = argparse.ArgumentParser(
        description="Bridge RNA-seq preprocessing to DE/GO analysis",
        epilog="""
Examples:
  # Auto-config mode (recommended):
  %(prog)s --project-id mouse-chd8 --skip-de --yes
  
  # With explicit config:
  %(prog)s --config config/projects/paths_mouse_chd8.yaml --project-id mouse-chd8 --skip-de --yes
  
  # With explicit paths:
  %(prog)s --rnaseq-output /path/to/output --de-pipeline /path/to/de --project-id mouse-chd8
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--config',
        type=Path,
        help='Path to YAML config file (auto-generated if not provided)'
    )
    parser.add_argument(
        '--rnaseq-output',
        type=Path,
        help='RNA-seq pipeline output directory (auto-detected from config if not provided)'
    )
    parser.add_argument(
        '--de-pipeline',
        type=Path,
        help='DE/GO pipeline directory (auto-detected from config if not provided)'
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
    parser.add_argument(
        '--yes',
        action='store_true',
        help='Assume yes to all prompts (non-interactive mode)'
    )
    
    args = parser.parse_args()
    
    # Auto-generate config if not provided
    if not args.config and AUTO_CONFIG_AVAILABLE:
        print(f"📝 No config provided, attempting auto-generation for {args.project_id}...")
        
        try:
            result = ensure_bridge_config(
                args.project_id,
                rnaseq_output=str(args.rnaseq_output) if args.rnaseq_output else None,
                de_pipeline=str(args.de_pipeline) if args.de_pipeline else None,
                force=False
            )
            
            if result['status'] in ['created', 'exists']:
                args.config = Path(result['config_path'])
                print(f"✅ Using config: {args.config}")
                
                if result['status'] == 'created':
                    print(f"   (Auto-generated from project summary)")
            else:
                print(f"⚠️  Auto-config failed: {result['message']}")
                print(f"   Continuing with CLI arguments...")
        
        except Exception as e:
            print(f"⚠️  Auto-config error: {e}")
            print(f"   Continuing with CLI arguments...")
    
    # Load config file if provided or auto-generated
    config = {}
    if args.config:
        config = load_path_config(args.config)
        if config and not AUTO_CONFIG_AVAILABLE:
            print(f"📁 Loaded config: {args.config}")
    
    # Resolve paths with fallback order: CLI > ENV > config > defaults
    rnaseq_output = args.rnaseq_output or \
                    (Path(os.getenv('RNASEQ_OUTPUT')) if os.getenv('RNASEQ_OUTPUT') else None) or \
                    (Path(config['rnaseq_output']) if 'rnaseq_output' in config else None)
    
    de_pipeline = args.de_pipeline or \
                  (Path(os.getenv('DE_PIPELINE')) if os.getenv('DE_PIPELINE') else None) or \
                  (Path(config['de_pipeline']) if 'de_pipeline' in config else None)
    
    # Validate required paths
    if not rnaseq_output:
        print("ERROR: RNA-seq output directory not specified.")
        print("  Provide via: --rnaseq-output, RNASEQ_OUTPUT env var, or config file")
        sys.exit(1)
    
    if not de_pipeline:
        print("ERROR: DE/GO pipeline directory not specified.")
        print("  Provide via: --de-pipeline, DE_PIPELINE env var, or config file")
        sys.exit(1)
    
    # Convert to Path objects and resolve
    rnaseq_output = Path(rnaseq_output).resolve()
    de_pipeline = Path(de_pipeline).resolve()
    
    # Validate paths exist
    if not rnaseq_output.exists():
        print(f"ERROR: RNA-seq output not found: {rnaseq_output}")
        sys.exit(1)
    
    if not de_pipeline.exists():
        print(f"ERROR: DE/GO pipeline not found: {de_pipeline}")
        sys.exit(1)
    
    # Execute bridge
    bridge = PipelineBridge(
        rnaseq_output,
        de_pipeline,
        args.project_id,
        config=config,
        assume_yes=args.yes
    )
    
    bridge.run(dry_run=args.dry_run, skip_de=args.skip_de)

if __name__ == '__main__':
    main()
