#!/usr/bin/env python3
"""
Sample Sheet Validation Tool

Validates sample sheet CSV files for RNA-seq pipeline to ensure:
- Required columns are present
- Data types are correct
- FASTQ files exist
- No duplicate sample IDs
- Proper formatting

Usage:
    python src/utils/validate_sample_sheet.py config/samples/master.csv
    python src/utils/validate_sample_sheet.py config/samples/master.csv --project H2O2_human_2025
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional


class SampleSheetValidator:
    """Validates sample sheet for RNA-seq pipeline"""
    
    # Required columns (must be present)
    REQUIRED_COLUMNS = [
        'project_id',
        'sample_id',
        'fastq_1',
        'fastq_2'
    ]
    
    # Optional columns (recommended but not required)
    OPTIONAL_COLUMNS = [
        'sample_name',
        'condition',
        'replicate',
        'sequencing_platform',
        'library_type',
        'read_type',
        'species',
        'genome_build',
        'batch',
        'tissue',
        'treatment',
        'time_point',
        'notes'
    ]
    
    # Valid values for specific columns
    VALID_VALUES = {
        'read_type': ['single-end', 'paired-end'],
        'library_type': ['RNA-seq', 'ChIP-seq', 'WGS', 'WES', 'ATAC-seq'],
    }
    
    def __init__(self, sample_sheet_path: str, check_files: bool = True):
        """
        Initialize validator
        
        Args:
            sample_sheet_path: Path to sample sheet CSV
            check_files: Whether to check FASTQ file existence
        """
        self.sample_sheet_path = Path(sample_sheet_path)
        self.check_files = check_files
        self.errors = []
        self.warnings = []
        self.df = None
        
    def validate(self, project_id: Optional[str] = None) -> Tuple[bool, List[str], List[str]]:
        """
        Run all validations
        
        Args:
            project_id: If provided, only validate samples for this project
            
        Returns:
            Tuple of (is_valid, errors, warnings)
        """
        # Read CSV
        if not self._read_csv():
            return False, self.errors, self.warnings
        
        # Filter by project if specified
        if project_id:
            original_count = len(self.df)
            self.df = self.df[self.df['project_id'] == project_id]
            if len(self.df) == 0:
                self.errors.append(f"No samples found for project_id: {project_id}")
                return False, self.errors, self.warnings
            self.warnings.append(f"Filtered to project {project_id}: {len(self.df)}/{original_count} samples")
        
        # Run validations
        self._validate_columns()
        self._validate_required_fields()
        self._validate_duplicates()
        self._validate_data_types()
        self._validate_valid_values()
        
        if self.check_files:
            self._validate_fastq_files()
        
        is_valid = len(self.errors) == 0
        return is_valid, self.errors, self.warnings
    
    def _read_csv(self) -> bool:
        """Read and parse CSV file"""
        if not self.sample_sheet_path.exists():
            self.errors.append(f"Sample sheet not found: {self.sample_sheet_path}")
            return False
        
        try:
            self.df = pd.read_csv(self.sample_sheet_path)
            if len(self.df) == 0:
                self.errors.append("Sample sheet is empty")
                return False
            return True
        except Exception as e:
            self.errors.append(f"Failed to read CSV: {e}")
            return False
    
    def _validate_columns(self):
        """Check required columns are present"""
        missing = set(self.REQUIRED_COLUMNS) - set(self.df.columns)
        if missing:
            self.errors.append(f"Missing required columns: {', '.join(missing)}")
        
        # Warn about unknown columns
        known_columns = set(self.REQUIRED_COLUMNS + self.OPTIONAL_COLUMNS)
        unknown = set(self.df.columns) - known_columns
        if unknown:
            self.warnings.append(f"Unknown columns (will be ignored): {', '.join(unknown)}")
    
    def _validate_required_fields(self):
        """Check required fields are not empty"""
        for col in self.REQUIRED_COLUMNS:
            if col not in self.df.columns:
                continue
            
            null_count = self.df[col].isnull().sum()
            if null_count > 0:
                self.errors.append(f"Column '{col}' has {null_count} empty values")
            
            # Check for whitespace-only values
            if self.df[col].dtype == 'object':
                empty_count = (self.df[col].str.strip() == '').sum()
                if empty_count > 0:
                    self.errors.append(f"Column '{col}' has {empty_count} whitespace-only values")
    
    def _validate_duplicates(self):
        """Check for duplicate sample IDs"""
        if 'sample_id' not in self.df.columns:
            return
        
        duplicates = self.df[self.df.duplicated(subset=['sample_id'], keep=False)]
        if len(duplicates) > 0:
            dup_ids = duplicates['sample_id'].unique()
            self.errors.append(
                f"Duplicate sample_id found: {', '.join(dup_ids)} "
                f"(appears in {len(duplicates)} rows)"
            )
    
    def _validate_data_types(self):
        """Validate data types for specific columns"""
        # Replicate should be numeric
        if 'replicate' in self.df.columns:
            non_numeric = self.df[pd.to_numeric(self.df['replicate'], errors='coerce').isnull()]
            if len(non_numeric) > 0:
                samples = non_numeric['sample_id'].tolist()[:5]
                self.warnings.append(
                    f"Non-numeric replicate values in samples: {', '.join(samples)}"
                    f"{' ...' if len(non_numeric) > 5 else ''}"
                )
        
        # Batch should be numeric
        if 'batch' in self.df.columns:
            non_numeric = self.df[pd.to_numeric(self.df['batch'], errors='coerce').isnull()]
            if len(non_numeric) > 0:
                samples = non_numeric['sample_id'].tolist()[:5]
                self.warnings.append(
                    f"Non-numeric batch values in samples: {', '.join(samples)}"
                    f"{' ...' if len(non_numeric) > 5 else ''}"
                )
    
    def _validate_valid_values(self):
        """Check values against valid lists"""
        for col, valid_values in self.VALID_VALUES.items():
            if col not in self.df.columns:
                continue
            
            invalid = self.df[~self.df[col].isin(valid_values + [None, ''])]
            if len(invalid) > 0:
                invalid_vals = invalid[col].unique()
                self.warnings.append(
                    f"Column '{col}' has invalid values: {', '.join(map(str, invalid_vals))}. "
                    f"Valid values: {', '.join(valid_values)}"
                )
    
    def _validate_fastq_files(self):
        """Check if FASTQ files exist"""
        missing_files = []
        
        for idx, row in self.df.iterrows():
            sample_id = row['sample_id']
            
            # Check R1
            if pd.notna(row['fastq_1']):
                fastq_1 = Path(row['fastq_1'])
                if not fastq_1.exists():
                    missing_files.append(f"{sample_id}: R1 not found - {fastq_1}")
            
            # Check R2 (if paired-end)
            if pd.notna(row['fastq_2']):
                fastq_2 = Path(row['fastq_2'])
                if not fastq_2.exists():
                    missing_files.append(f"{sample_id}: R2 not found - {fastq_2}")
        
        if missing_files:
            if len(missing_files) <= 10:
                for msg in missing_files:
                    self.errors.append(msg)
            else:
                self.errors.append(
                    f"FASTQ files not found for {len(missing_files)} entries. "
                    f"First 10: {'; '.join(missing_files[:10])}"
                )
    
    def get_summary(self) -> Dict:
        """Get validation summary statistics"""
        if self.df is None:
            return {}
        
        summary = {
            'total_samples': len(self.df),
            'projects': self.df['project_id'].nunique() if 'project_id' in self.df.columns else 0,
            'conditions': self.df['condition'].nunique() if 'condition' in self.df.columns else 0,
        }
        
        # Count by project
        if 'project_id' in self.df.columns:
            summary['samples_per_project'] = self.df['project_id'].value_counts().to_dict()
        
        # Count by condition
        if 'condition' in self.df.columns:
            summary['samples_per_condition'] = self.df['condition'].value_counts().to_dict()
        
        return summary


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Validate RNA-seq sample sheet CSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate entire sample sheet
  python src/utils/validate_sample_sheet.py config/samples/master.csv
  
  # Validate specific project only
  python src/utils/validate_sample_sheet.py config/samples/master.csv --project H2O2_human_2025
  
  # Skip FASTQ file existence check (faster)
  python src/utils/validate_sample_sheet.py config/samples/master.csv --no-check-files
        """
    )
    
    parser.add_argument(
        'sample_sheet',
        help='Path to sample sheet CSV file'
    )
    parser.add_argument(
        '--project', '-p',
        help='Only validate samples for this project_id'
    )
    parser.add_argument(
        '--no-check-files',
        action='store_true',
        help='Skip FASTQ file existence check'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Only show errors, suppress warnings and summary'
    )
    
    args = parser.parse_args()
    
    # Run validation
    validator = SampleSheetValidator(
        args.sample_sheet,
        check_files=not args.no_check_files
    )
    
    is_valid, errors, warnings = validator.validate(project_id=args.project)
    
    # Print results
    if not args.quiet:
        print("=" * 80)
        print(f"SAMPLE SHEET VALIDATION: {args.sample_sheet}")
        print("=" * 80)
    
    # Errors
    if errors:
        print(f"\n❌ ERRORS ({len(errors)}):")
        for error in errors:
            print(f"  - {error}")
    
    # Warnings
    if warnings and not args.quiet:
        print(f"\n⚠️  WARNINGS ({len(warnings)}):")
        for warning in warnings:
            print(f"  - {warning}")
    
    # Summary
    if is_valid and not args.quiet:
        summary = validator.get_summary()
        print(f"\n✅ VALIDATION PASSED")
        print(f"\nSummary:")
        print(f"  Total samples: {summary.get('total_samples', 0)}")
        print(f"  Projects: {summary.get('projects', 0)}")
        print(f"  Conditions: {summary.get('conditions', 0)}")
        
        if 'samples_per_project' in summary:
            print(f"\n  Samples per project:")
            for proj, count in summary['samples_per_project'].items():
                print(f"    {proj}: {count}")
    
    print("=" * 80)
    
    # Exit code
    sys.exit(0 if is_valid else 1)


if __name__ == '__main__':
    main()
