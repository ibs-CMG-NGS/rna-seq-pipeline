# Sample Metadata Integration Test Guide

## Overview
Step 4 of sample sheet integration is complete. This guide shows how to test the metadata flow from sample sheet to manifest.json.

## What Was Implemented

### 1. Manifest Generator Changes
- **File**: `scripts/standardization/generate_manifest.py`
- **Changes**:
  - Added `sample_metadata` parameter to `__init__`
  - New method `_extract_metadata()` extracts 13 relevant fields
  - Manifest JSON now includes `"sample_metadata": {...}` section
  - Optional pandas support with `HAS_PANDAS` flag

### 2. Snakefile Changes
- **File**: `Snakefile`
- **Changes**:
  - Added `json` import
  - Updated `generate_manifest` rule with new parameter:
    ```python
    params:
        sample_metadata=lambda wildcards: json.dumps(SAMPLE_METADATA.get(wildcards.sample, {})) if USE_SAMPLE_SHEET else "{}"
    ```
  - Shell command now passes `--sample-metadata '{params.sample_metadata}'`

### 3. Metadata Fields Extracted
The following 13 fields are extracted from sample sheet (if present and non-empty):
1. `sample_name` - Human-readable sample name
2. `condition` - Experimental condition (e.g., Control, H2O2_100uM)
3. `replicate` - Biological replicate number
4. `sequencing_platform` - Platform used (e.g., Illumina NovaSeq)
5. `library_type` - Library preparation type
6. `read_type` - Read type (paired-end, single-end)
7. `species` - Species (e.g., Homo sapiens)
8. `genome_build` - Genome version (e.g., GRCh38)
9. `batch` - Sequencing batch
10. `tissue` - Tissue type (e.g., PFC, HPC)
11. `treatment` - Treatment details
12. `time_point` - Time point in experiment
13. `notes` - Additional notes

## Test Procedure

### Server Setup (On NGS Server)

```bash
# 1. Navigate to pipeline directory
cd ~/ngs_pipeline/rna-seq-pipeline

# 2. Pull latest changes
git pull origin main

# 3. Verify commits
git log --oneline -3
# Should see: "feat: add sample metadata to manifest.json"

# 4. Check sample sheet exists
ls -lh config/samples/master.csv
```

### Test 1: Validate Sample Sheet

```bash
# Activate conda environment
conda activate rna-seq-pipeline

# Validate H2O2_human_2025 project samples
python src/utils/validate_sample_sheet.py \
    config/samples/master.csv \
    --project H2O2_human_2025 \
    --check-files

# Expected output:
# ✅ VALIDATION PASSED
# Total samples: 15
# Conditions: 5 (Control, Day1, Day3, H2O2_100uM, H2O2_200uM)
```

### Test 2: Dry-Run with Sample Sheet Mode

```bash
# Dry-run to see what will be executed
snakemake \
    --configfile config/projects/config_H2O2_human.yaml \
    --config use_sample_sheet=true \
    --cores 1 \
    --dry-run \
    2>&1 | head -50

# Check for:
# - "Sample Loading: Sample Sheet"
# - "Found 15 samples"
# - "Conditions: Control, Day1, Day3, H2O2_100uM, H2O2_200uM"
```

### Test 3: Generate Single Manifest (Manual Test)

```bash
# Test manifest generation for one sample
python scripts/standardization/generate_manifest.py \
    --sample-dir /home/ngs/data/results/H2O2_human_2025/h_RNA_Cont_1/rna-seq \
    --sample-id h_RNA_Cont_1 \
    --project-id H2O2_human_2025 \
    --pipeline-type rna-seq \
    --sample-metadata '{"condition": "Control", "replicate": 1, "tissue": "astrocyte"}'

# Check manifest.json
cat /home/ngs/data/results/H2O2_human_2025/h_RNA_Cont_1/rna-seq/final_outputs/manifest.json | python -m json.tool

# Should see sample_metadata section:
# {
#   "sample_id": "h_RNA_Cont_1",
#   "sample_metadata": {
#     "condition": "Control",
#     "replicate": 1,
#     "tissue": "astrocyte"
#   },
#   ...
# }
```

### Test 4: Full Pipeline Run (One Sample)

```bash
# Run pipeline for one sample to test integration
snakemake \
    --configfile config/projects/config_H2O2_human.yaml \
    --config use_sample_sheet=true \
    --cores 8 \
    --until generate_manifest \
    --forcerun generate_manifest \
    h_RNA_Cont_1

# Check generated manifest
cat /home/ngs/data/results/H2O2_human_2025/h_RNA_Cont_1/rna-seq/final_outputs/manifest.json | python -m json.tool | grep -A 10 sample_metadata
```

### Test 5: All Samples Pipeline Run

```bash
# Run full pipeline for all 15 samples
snakemake \
    --configfile config/projects/config_H2O2_human.yaml \
    --config use_sample_sheet=true \
    --cores 40 \
    --rerun-incomplete

# Monitor progress
tail -f logs/snakemake_*.log
```

## Validation Checks

### 1. Check Manifest Content

```bash
# For each sample, verify metadata is present
for sample in h_RNA_Cont_1 h_RNA_Cont_2 h_RNA_Cont_3 h_RNA_100_1 h_RNA_100_2 h_RNA_100_3; do
    echo "=== $sample ==="
    cat /home/ngs/data/results/H2O2_human_2025/$sample/rna-seq/final_outputs/manifest.json \
        | python -m json.tool \
        | grep -A 5 '"sample_metadata"'
    echo
done
```

### 2. Verify Metadata Fields

```python
# In Python
import json
from pathlib import Path

# Load manifest
manifest_path = Path("/home/ngs/data/results/H2O2_human_2025/h_RNA_Cont_1/rna-seq/final_outputs/manifest.json")
with open(manifest_path) as f:
    manifest = json.load(f)

# Check metadata
print("Sample Metadata:")
print(json.dumps(manifest.get('sample_metadata', {}), indent=2))

# Expected fields for h_RNA_Cont_1:
# {
#   "condition": "Control",
#   "replicate": 1,
#   "species": "Homo sapiens",
#   "genome_build": "GRCh38",
#   ...
# }
```

### 3. Compare Metadata Across Conditions

```bash
# Extract conditions from all manifests
for condition in Control Day1 Day3 H2O2_100uM H2O2_200uM; do
    echo "=== $condition ==="
    find /home/ngs/data/results/H2O2_human_2025 -name "manifest.json" -exec grep -l "\"condition\": \"$condition\"" {} \;
done
```

## Expected Results

### Success Criteria
✅ All 15 samples have manifest.json with `sample_metadata` section  
✅ Metadata contains at least: `condition`, `replicate`  
✅ Empty/NaN fields are excluded from metadata  
✅ Metadata matches sample sheet content  
✅ No errors in manifest generation logs  

### Example Manifest Structure
```json
{
  "sample_id": "h_RNA_Cont_1",
  "project_id": "H2O2_human_2025",
  "pipeline_type": "rna-seq",
  "execution_date": "2025-06-15T10:30:00",
  "status": "completed",
  
  "sample_metadata": {
    "condition": "Control",
    "replicate": 1,
    "species": "Homo sapiens",
    "genome_build": "GRCh38",
    "tissue": "astrocyte"
  },
  
  "final_outputs": [...],
  "qc_metrics": {...},
  "next_steps": [...]
}
```

## Troubleshooting

### Issue 1: Metadata Not in Manifest
**Symptom**: `sample_metadata` is empty `{}`  
**Cause**: Sample sheet not used or metadata not passed  
**Solution**:
```bash
# Verify use_sample_sheet is true
grep use_sample_sheet config/projects/config_H2O2_human.yaml

# Or override in command:
snakemake --config use_sample_sheet=true ...
```

### Issue 2: Pandas Import Error
**Symptom**: Script fails with "No module named 'pandas'"  
**Cause**: Pandas not in conda environment  
**Solution**:
```bash
# Update environment
conda env update -f environment.yaml

# Verify pandas installed
conda list pandas
```

### Issue 3: JSON Decode Error
**Symptom**: "Invalid JSON in --sample-metadata"  
**Cause**: Shell escaping issue  
**Solution**: Use single quotes around JSON:
```bash
--sample-metadata '{"condition": "Control", "replicate": 1}'
```

### Issue 4: Metadata Fields Missing
**Symptom**: Expected fields not in manifest  
**Cause**: Fields are NaN or empty in sample sheet  
**Solution**: Check sample sheet CSV:
```bash
# Find empty fields for sample
grep "^H2O2_human_2025,Ctrl_1," config/samples/master.csv
```

## Next Steps

After successful testing:

1. **Step 5**: Create `SAMPLE_SHEET_GUIDE.md` documentation
   - Schema explanation
   - Field descriptions
   - Usage examples
   - Best practices

2. **Enhancement**: Condition-based QC grouping
   - Group MultiQC reports by condition
   - Separate QC summaries per condition
   - Agent can query by experimental group

3. **Agent Integration**: Query metadata
   - "What condition is this sample?"
   - "Show me all replicates for Control"
   - "Compare QC across conditions"

## Commit Information

- **Branch**: main
- **Commit**: 6ca1aa8
- **Message**: "feat: add sample metadata to manifest.json"
- **Files Changed**:
  - Snakefile (+6 lines)
  - scripts/standardization/generate_manifest.py (+52 lines)
