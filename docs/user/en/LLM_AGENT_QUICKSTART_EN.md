# LLM Agent Quick Start Guide

Control the RNA-seq pipeline through natural language conversation.  
**FASTQ detection → Project setup → Validation → Execution → Monitoring** — all in one chat session.

---

## 1. Prerequisites (one-time setup)

```bash
# Install Ollama (skip if already installed on the server)
curl -fsSL https://ollama.com/install.sh | sh

# Start the Ollama service
ollama serve &

# Download the recommended model (~20 GB, GPU server)
ollama pull qwen2.5:32b

# Install the Python client
conda activate rna-seq-pipeline
pip install ollama
```

> **Shared server users**: Ollama is already running. Run `ollama list` to see available models.

---

## 2. Start the Agent

```bash
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline

# Basic (for running a new pipeline)
python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen2.5:32b

# With an existing project (for QC analysis)
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --interactive \
  --model qwen2.5:32b
```

---

## 3. End-to-End Workflow Example

### Step 1 — Detect FASTQ files
```
You: Find FASTQ files in /data_3tb/shared/chd8-rna-seq-raw-data/fastq/

Agent: Found 76 files, 38 samples (149 GB, paired-end)
       Pattern: Chd8_{tissue}_{age}_{sex}_{replicate}_1/2.fq.gz
```

### Step 2 — Create project config
```
You: Create a new project. ID: mouse-chd8-2026,
     save results to /data_3tb/shared/output/mouse-chd8-2026/

Agent: config/projects/mouse-chd8-2026.yaml created ✅
```

### Step 3 — Estimate resources
```
You: How much resources do I need?

Agent: Input 149 GB → estimated output 746 GB
       Estimated time: 5-10 hours (8 cores)
       RAM: ~30 GB (STAR), disk free 2145 GB ✅
```

### Step 4 — Generate sample sheet
```
You: Generate the sample sheet

Agent: 38 samples classified automatically
       wildtype: 20, heterozygous: 18
       → config/samples/mouse-chd8-2026.tsv saved ✅
```

### Step 5 — Validate input data
```
You: Validate input data for config/projects/mouse-chd8-2026.yaml

Agent: ✅ FASTQ files confirmed
       ✅ Disk space sufficient
       ⚠️  Genome path needs verification: /data_3tb/shared/genome/star_index/
```

### Step 6 — Dry run
```
You: Do a dry run

Agent: 270 jobs across 11 rules
       cutadapt: 38 / fastqc_raw: 76 / star_align: 38 / ...
       No issues ✅ Ready for real execution
```

### Step 7 — Run the pipeline
```
You: Run the pipeline with 16 cores

Agent: Pipeline started 🚀
       snakemake --configfile ... --cores 16
       Check progress with "show pipeline status"
```

### Step 8 — Monitor progress
```
You: Show pipeline status

Agent: Progress: 45% (cutadapt done, star_align running)
       Completed: cutadapt (38/38), fastqc (76/76)
       Running: star_align (17/38)
```

### Step 9 — Review QC results
```
You: How does the QC look?

Agent: All 38 samples passed (100%)
       Average mapping rate: 92.4%, average assignment: 69.9%
       Ready for DE analysis ✅
```

---

## 4. QC Analysis Query Examples

```
# Overall status
You: Show me the QC status
You: How many samples passed QC?

# Condition comparison
You: Compare wildtype vs heterozygous
You: Compare only in HPC tissue
You: Compare male vs female

# Multi-axis analysis
You: What experimental groups (axes) are available?
You: List HPC wildtype male samples

# Failed samples
You: Any samples that failed QC?
You: Show details for sample Chd8_HPC_10M_W_1

# DE analysis
You: Start DE analysis
```

---

## 5. Command-Line Options

| Option | Description | Example |
|--------|-------------|---------|
| `--interactive` | Interactive chat mode (recommended) | — |
| `--message "..."` | Single-query mode | `--message "Show QC status"` |
| `--model` | Model selection | `--model qwen2.5:32b` |
| `--project-summary` | Path to existing project summary | `--project-summary .../project_summary.json` |
| `--rnaseq-output` | Pipeline output directory | `--rnaseq-output /data_3tb/.../output` |
| `--ollama-host` | Remote Ollama server address | `--ollama-host http://gpu-server:11434` |

---

## 6. Model Selection Guide

| Model | VRAM/RAM | Speed | Quality | Tool Calling | Best For |
|-------|----------|-------|---------|--------------|----------|
| **qwen2.5:32b** | 20 GB VRAM | ⚡⚡ | ⭐⭐⭐⭐⭐ | Native (GPU) | **Default — recommended** |
| llama3.1:8b | 8 GB RAM | ⚡⚡⚡ | ⭐⭐⭐⭐ | Native | No GPU available |
| mistral:7b | 8 GB RAM | ⚡⚡⚡⚡ | ⭐⭐⭐ | Text pattern | Quick status checks |
| llama3.1:70b | 42 GB RAM | ⚡ | ⭐⭐⭐⭐⭐ | Native | High-memory server |

---

## 7. Available Tools (17 total)

**QC Analysis (Phase 6)**
- `get_project_status` — Overall QC summary
- `compare_conditions` — Compare experimental conditions
- `get_failed_samples` — List samples that failed QC
- `get_sample_details` — Per-sample metadata and metrics
- `list_conditions` — Show all experimental groups
- `start_de_analysis` — Trigger DE/GO analysis pipeline

**Multi-axis Analysis (Phase 6.1)**
- `get_sample_axes` — List experimental axes (tissue, sex, genotype, ...)
- `compare_by_axis` — Compare QC metrics by any axis with optional filters
- `filter_samples` — Filter samples by multiple axes simultaneously

**Pipeline Execution (Phase 8A)**
- `create_project_config` — Generate project config YAML
- `detect_fastq_files` — Scan directory for FASTQ files
- `validate_input_data` — Pre-flight validation checks
- `run_pipeline` — Execute Snakemake workflow (dry-run or real)

**Monitoring & Resources (Phase 8B)**
- `monitor_pipeline` — Check pipeline execution status
- `create_sample_sheet` — Auto-generate sample metadata TSV
- `estimate_resources` — Estimate runtime, disk, and RAM usage

---

## 8. Troubleshooting

**Ollama not responding**
```bash
ollama serve          # Start the service
ollama list           # List installed models
```

**Model not found**
```bash
ollama pull qwen2.5:32b    # Recommended model (~20 GB)
ollama pull llama3.1:8b    # Lightweight alternative (~5 GB)
```

**Slow responses**
```bash
nvidia-smi            # Verify GPU is detected
# If no GPU, switch to a lighter model: --model llama3.1:8b
```

**Config file path error**  
Both relative paths (from the project root) and absolute paths are accepted.

---

## Further Documentation

- Ollama setup: `docs/user/en/OLLAMA_SETUP_EN.md`
- Full pipeline guide: `docs/user/en/PIPELINE_GUIDE_EN.md`
- Developer reference: `docs/developer/PHASE6_LLM_INTEGRATION.md`
- Phase 8 execution details: `docs/developer/PHASE8_PIPELINE_EXECUTION.md`
