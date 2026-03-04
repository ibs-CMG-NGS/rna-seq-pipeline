# RNA-seq Analysis Pipeline Guide

A Snakemake-based RNA-seq data analysis pipeline.

## Pipeline Steps

1. **Quality Control (FastQC)** — Raw data quality assessment
2. **Adapter Trimming (Cutadapt)** — Remove adapters and low-quality bases
3. **Alignment (STAR)** — Map reads to the reference genome
4. **Quantification (featureCounts)** — Count reads per gene
5. **QC Report** — Generate a comprehensive HTML report

---

## 🤖 Run via LLM Agent (Recommended)

Control the entire pipeline through natural language — no command memorization required.

```bash
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline

python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen2.5:32b
```

**Example conversation:**
```
You: Find FASTQ files in /data_3tb/shared/chd8-rna-seq-raw-data/fastq/
Agent: 38 samples detected (149 GB, paired-end) ✅

You: Create a new project. ID: my-project, save to /data_3tb/shared/output/my-project/
Agent: config/projects/my-project.yaml created ✅

You: How much disk space and time do I need?
Agent: ~746 GB output, 5-10 hours, disk/RAM sufficient ✅

You: Generate the sample sheet
Agent: wildtype 20, heterozygous 18 → TSV saved ✅

You: Validate the input data
Agent: FASTQ ✅, disk ✅, genome path ⚠️ needs checking

You: Do a dry run
Agent: 270 jobs, 11 rules — no issues ✅

You: Run with 16 cores
Agent: Pipeline started 🚀

You: Show pipeline status
Agent: 45% complete, star_align 17/38 in progress...
```

> **Quick start guide**: `docs/user/en/LLM_AGENT_QUICKSTART_EN.md`

---

## 🔧 Manual Execution

### 1. Set Up the Environment

```bash
conda env create -f environment.yaml
conda activate rna-seq-pipeline
```

### 2. Prepare Data

Place raw FASTQ files in `data/raw/`:
```
data/raw/
├── SampleA_1.fastq.gz   # R1
├── SampleA_2.fastq.gz   # R2
├── SampleB_1.fastq.gz
└── SampleB_2.fastq.gz
```

Edit `samples.tsv` to list all samples:
```tsv
sample
SampleA
SampleB
```

### 3. Configure Parameters

Edit `config/config.yaml`:

```yaml
# Reference files
star_index: "genome/star_index/"
annotation_gtf: "genome/genes.gtf"

# Adapter sequences (Illumina TruSeq)
adapter_r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Quality filters
quality_cutoff: 20
min_read_length: 20

# Threading
star_threads: 8
featurecounts_threads: 4

# QC report
generate_qc_report: true
qc_report_output: "results/qc_report.html"
qc_top_genes: 10
```

### 4. Dry Run (Test)

```bash
snakemake -n
```

### 5. Run the Pipeline

```bash
# Use all available cores
snakemake --cores all

# Specify core count
snakemake --cores 8

# Run up to a specific step
snakemake --cores 8 results/counts/counts_matrix.txt
```

### 6. Regenerate QC Report Only

```bash
snakemake --cores 1 results/qc_report.html --force
```

---

## Performance Tuning

### Resource Requirements per Step

| Step | CPU | Memory | Notes |
|------|-----|--------|-------|
| FastQC | Low (1 thread) | ~500 MB | I/O bound |
| Cutadapt | Medium (4 threads) | ~1-2 GB | CPU bound |
| **STAR align** | **High (8-16 threads)** | **~30 GB** | **Memory intensive** |
| featureCounts | Medium (4-8 threads) | ~2-3 GB | I/O bound |

### Recommended Settings by Server

| Server Specs | star_threads | featurecounts_threads | --jobs | Notes |
|---|---|---|---|---|
| 8 cores, 32 GB RAM | 6 | 4 | 1 | Memory-limited |
| 12 cores, 62 GB RAM | 12 | 8 | 2 | Typical lab server |
| 24 cores, 128 GB RAM | 16 | 12 | 4 | Comfortable |
| 48+ cores, 256 GB+ RAM | 24 | 16 | 8 | Large-scale |

### Snakemake Resource Control

```bash
# Safe default (prevents OOM)
snakemake --cores 16 --jobs 2

# 3 samples in parallel
snakemake --cores 18 --jobs 3 --resources mem_gb=60

# Conservative (one at a time)
snakemake --cores 12 --jobs 1

# Network file system (NFS): add latency-wait
snakemake --cores 16 --jobs 2 --latency-wait 120
```

---

## Output Files

```
results/
├── trimmed/                        # Adapter-trimmed FASTQ
│   ├── {sample}_1.fastq.gz
│   └── {sample}_2.fastq.gz
├── aligned/                        # STAR alignment results
│   └── {sample}/
│       ├── Aligned.sortedByCoord.out.bam
│       └── Log.final.out
├── counts/                         # featureCounts results
│   ├── counts_matrix.txt           # Raw featureCounts output
│   ├── counts_matrix.txt.summary   # Mapping statistics
│   └── counts_matrix_clean.csv     # Clean matrix for DE analysis ⭐
└── qc_report.html                  # QC report (auto-generated)
```

### Downstream Analysis

`results/counts/counts_matrix_clean.csv` is ready for DESeq2 / edgeR / limma-voom:
- Rows: gene IDs
- Columns: sample names (cleaned)
- No metadata columns (Chr, Start, End removed)

```csv
Geneid,SampleA,SampleB,SampleC,...
ENSMUSG00000104478,0,0,0,...
ENSMUSG00000086053,0,2,4,...
```

Point `count_data_path` in your DE pipeline config to this file.

---

## QC Report Contents

The HTML report includes:
- ✂️ Adapter trimming statistics (per sample)
- 🎯 Alignment quality (mapping rate, progress bars)
- 🧮 Gene quantification statistics
- 📊 Expressed gene counts
- 🔝 Top N highly expressed genes (configurable)
- 💾 File size information

---

## Troubleshooting

### Pipeline Fails — Check Logs

```bash
cat logs/cutadapt/{sample}.log
cat logs/star/{sample}.log
cat logs/featurecounts.log

# Re-run a specific step
snakemake --cores 8 --forcerun star_align
```

### MissingOutputException (Network File System)

```
MissingOutputException: Job completed successfully, but some output files are missing.
```

This happens on NFS/SMB mounts due to metadata sync delays.

```bash
snakemake --cores 16 --jobs 2 --latency-wait 60   # mild
snakemake --cores 16 --jobs 2 --latency-wait 180  # severe
```

### Command Not Found (cutadapt, STAR, etc.)

```bash
# Ensure the conda environment is active
conda activate rna-seq-pipeline

# Verify tools are installed
conda list | grep -E "cutadapt|fastqc|star|subread"

# Install missing tools
conda install -c bioconda cutadapt fastqc star subread samtools
```

### FASTQ File Issues

```bash
python3 src/check_fastq.py    # Validate file integrity
python3 src/fix_fastq.py      # Fix non-standard formats (if needed)
```

---

## Adding New Samples

1. Place new FASTQ files in `data/raw/`
2. Add sample names to `samples.tsv`
3. Re-run Snakemake — it will process only the new samples
4. The QC report automatically includes all samples

---

## Next Steps After Analysis

1. Open `results/qc_report.html` in a browser
2. Use `results/counts/counts_matrix_clean.csv` for differential expression analysis
3. Check `docs/user/FASTQC_GUIDE.md` for detailed QC interpretation

---

## References

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [STAR Manual](https://github.com/alexdobin/STAR)
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [featureCounts (Subread)](http://subread.sourceforge.net/)
- [DESeq2](https://bioconductor.org/packages/DESeq2/)
