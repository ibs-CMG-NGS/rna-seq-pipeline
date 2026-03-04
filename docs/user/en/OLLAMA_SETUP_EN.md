# Ollama Setup Guide

## Overview

Ollama runs large language models locally, so your genomic data never leaves the server.  
This guide covers installation, model selection, and connecting to the RNA-seq pipeline agent.

## Why Run LLMs Locally?

- **Data security** — No genomic data sent to external APIs
- **Privacy** — Fully compliant with HIPAA/GDPR requirements
- **Cost** — Free to use, no per-token charges
- **Control** — Choose any model, adjust parameters freely

## Prerequisites

- Linux / WSL2 environment
- 8 GB+ RAM (20 GB VRAM recommended for GPU acceleration)
- 20 GB+ free disk space for models

---

## Installation

### 1. Install Ollama

```bash
curl -fsSL https://ollama.com/install.sh | sh

# Verify
ollama --version
```

The service runs at `http://localhost:11434` by default.

### 2. Start the Service

```bash
# Foreground (for testing)
ollama serve

# Background
ollama serve &

# As a system service (persistent across reboots)
sudo systemctl enable ollama
sudo systemctl start ollama
sudo systemctl status ollama
```

### 3. Download Models

**Recommended:**

```bash
# qwen2.5:32b — Recommended (GPU server, native tool calling, excellent quality)
ollama pull qwen2.5:32b

# llama3.1:8b — Lightweight alternative (CPU-only environments)
ollama pull llama3.1:8b

# mistral:7b — Fast, for quick status queries
ollama pull mistral:7b

# llama3.1:70b — Highest quality (requires 42 GB+ RAM)
ollama pull llama3.1:70b
```

**Model selection guide:**
- **qwen2.5:32b** — Default recommendation (GPU, 20 GB VRAM, native tool calling, multilingual)
- **llama3.1:8b** — Best for CPU-only setups (8 GB RAM, native tool calling)
- **mistral:7b** — Fastest for simple queries (8 GB RAM, text-pattern tool calling)
- **llama3.1:70b** — Best quality for complex analysis (42 GB RAM required)

### 4. Install the Python Client

```bash
conda activate rna-seq-pipeline
pip install ollama
```

### 5. Verify Everything Works

```bash
# Check the service is running
curl http://localhost:11434/api/tags

# List installed models
ollama list

# Quick chat test
ollama run qwen2.5:32b "Hello! Are you ready?"
```

---

## Using the RNA-seq Pipeline Agent

### Interactive Mode (Recommended)

```bash
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline

python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --interactive \
  --model qwen2.5:32b
```

**Example session:**
```
You: Show me the QC status

Agent: Your mouse-chd8 project has excellent quality:
  - Total samples: 38
  - Passed QC: 38 (100%)
  - Average mapping rate: 92.4%
  All samples are ready for downstream analysis!

You: Find FASTQ files in /data_3tb/shared/chd8-rna-seq-raw-data/fastq/

Agent: Found 76 files, 38 samples (149 GB, paired-end) ✅

You: Create a new project. ID: my-project-2026,
     save to /data_3tb/shared/output/my-project-2026/

Agent: config/projects/my-project-2026.yaml created ✅
```

### Single Query Mode

```bash
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --message "Which samples failed QC?"
```

### Remote Ollama Server

```bash
python scripts/standardization/llm_agent.py \
  --project-summary /path/to/project_summary.json \
  --rnaseq-output /path/to/output \
  --ollama-host http://gpu-server.local:11434 \
  --interactive
```

---

## Performance Comparison

| Model | Size | RAM/VRAM | Speed | Quality | Tool Calling | Best For |
|-------|------|----------|-------|---------|--------------|----------|
| mistral:7b | 4 GB | 8 GB RAM | Fast | Good | Text pattern | Quick status checks |
| llama3.1:8b | 5 GB | 8 GB RAM | Medium | Great | Native | CPU-only environments |
| **qwen2.5:32b** | 20 GB | 20 GB VRAM | Medium | **Excellent** | **Native (GPU)** | **Default — recommended** |
| llama3.1:70b | 40 GB | 42 GB RAM | Slow | Best | Native | High-memory servers |

---

## GPU Acceleration

Ollama automatically uses NVIDIA GPUs when available:

```bash
# Verify GPU is detected
nvidia-smi

# Monitor GPU usage during inference
watch -n 1 nvidia-smi
```

For the lab server (2× RTX 2080 Ti, 22 GB VRAM total), `qwen2.5:32b` runs fully on GPU.

---

## Advanced: Custom Model (Optional)

Create a Modelfile to tune the model for bioinformatics:

```bash
cat > Modelfile << 'MEOF'
FROM qwen2.5:32b

SYSTEM \"\"\"
You are an expert bioinformatics assistant specializing in RNA-seq analysis.
Focus on clear, actionable insights for genomics researchers.
\"\"\"

PARAMETER temperature 0.3
PARAMETER num_ctx 8192
MEOF

ollama create rna-seq-agent -f Modelfile

# Use the custom model
python scripts/standardization/llm_agent.py \
  --model rna-seq-agent \
  --interactive
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| `Connection refused` | Ollama not running | `ollama serve` |
| `model not found` | Model not downloaded | `ollama pull qwen2.5:32b` |
| Slow responses | No GPU / wrong model | Check `nvidia-smi`; use `llama3.1:8b` on CPU |
| Out of memory | Model too large | Use `llama3.1:8b` (8 GB RAM) instead |
| Tool not called | Wrong model version | Use `qwen2.5:32b` or `llama3.1:8b`+ |

---

## Security Notes

✅ All model inference runs on your local server  
✅ No genomic data sent to external services  
✅ Full audit trail in terminal output  
✅ HIPAA/GDPR compliant for sensitive data  

---

## Resources

- Ollama Documentation: https://ollama.com/
- Model Library: https://ollama.com/library
- Python Client: https://github.com/ollama/ollama-python
- Project Issues: https://github.com/ibs-CMG-NGS/rna-seq-pipeline/issues
