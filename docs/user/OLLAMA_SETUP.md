# Ollama Setup Guide for RNA-seq Pipeline

## Overview

This guide explains how to set up Ollama for local LLM-powered pipeline queries. Ollama allows you to run large language models locally, ensuring **data security** for sensitive genomic data.

## Why Local LLM?

- **Data Security**: Genomic data never leaves your server
- **Privacy**: No external API calls to cloud providers
- **Cost**: Free to use, no per-token charges
- **Control**: Full control over model selection and configuration

## Prerequisites

- Linux/WSL2 environment
- 8GB+ RAM (16GB+ recommended)
- 10GB+ free disk space for models

## Installation

### 1. Install Ollama

```bash
# Install Ollama (Linux/WSL)
curl -fsSL https://ollama.com/install.sh | sh

# Start Ollama service
ollama serve
```

The service will run on `http://localhost:11434` by default.

### 2. Download Models

**권장 모델:**

```bash
# qwen2.5:32b — 권장 (GPU 서버, native tool calling, 최고 품질)
ollama pull qwen2.5:32b

# llama3.1:8b — GPU 없는 환경 또는 경량 사용
ollama pull llama3.1:8b

# mistral:7b — 빠른 상태 확인용
ollama pull mistral:7b

# llama3.1:70b — 고성능 서버 (42GB+ RAM)
ollama pull llama3.1:70b
```

**모델 선택 가이드:**
- **qwen2.5:32b** — 기본 권장 (GPU 서버, 20GB VRAM, native tool calling, 한국어 우수)
- **llama3.1:8b** — GPU 없을 때 (8GB RAM, native tool calling)
- **mistral:7b** — 빠른 단순 쿼리 (8GB RAM, text-pattern tool calling)
- **llama3.1:70b** — 가장 고품질 (42GB RAM 필요)

### 3. Install Python Package

```bash
# Activate your RNA-seq environment
conda activate rna-seq-pipeline

# Install Ollama Python client
pip install ollama
```

### 4. Verify Installation

```bash
# Test Ollama is running
curl http://localhost:11434/api/tags

# Test model is available
ollama list

# Test chat functionality
ollama run llama3.1:8b "Hello, are you working?"
```

## Usage

### Interactive Mode (Recommended)

```bash
# Start interactive chat session
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --interactive

# The agent will start with:
# 🤖 RNA-seq Pipeline Agent
# Project: mouse-chd8
# Model: llama3.1:8b (local)
# Type 'quit' or 'exit' to end session.
# 
# You:
```

**Example conversation:**
```
You: Show me the QC status

Agent: Your mouse-chd8 project has excellent quality:
- Total samples: 38
- Passed QC: 38 (100%)
- Failed QC: 0
- Average mapping rate: 92.4%

All samples are ready for downstream analysis!

You: Compare wildtype and heterozygous conditions

Agent: Comparing wildtype vs heterozygous:

Wildtype (20 samples):
- Average mapping rate: 92.39%
- Average assignment rate: 69.81%

Heterozygous (18 samples):
- Average mapping rate: 92.49%
- Average assignment rate: 69.95%

The conditions show very similar QC metrics, indicating high-quality data
across both groups. This suggests no systematic technical differences.

You: Can you start DE analysis?

Agent: I can help you prepare for DE analysis. Here's what I'll do:

1. Verify all RNA-seq samples passed QC ✓
2. Copy counts matrix to DE pipeline
3. Generate metadata CSV from sample sheet
4. Create DE/GO config file

Would you like me to proceed? (yes/no)
```

### Single Query Mode

```bash
# Quick one-off queries
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --message "Show me failed samples"
```

### Custom Model Selection

```bash
# Use different model
python scripts/standardization/llm_agent.py \
  --project-summary /path/to/project_summary.json \
  --rnaseq-output /path/to/output \
  --model llama3.1:70b \
  --interactive
```

### Custom Ollama Host

```bash
# Connect to remote Ollama server
python scripts/standardization/llm_agent.py \
  --project-summary /path/to/project_summary.json \
  --rnaseq-output /path/to/output \
  --ollama-host http://gpu-server.local:11434 \
  --interactive
```

## Available Commands

The agent can execute these pipeline commands via natural language:

1. **Get project status** - "Show me the QC status", "How's the project doing?"
2. **Compare conditions** - "Compare wildtype vs knockout", "Show condition differences"
3. **Get failed samples** - "Which samples failed?", "Show me QC failures"
4. **Get sample details** - "Tell me about sample S001", "What's the status of S001?"
5. **List conditions** - "What conditions do I have?", "Show me all groups"
6. **Start DE analysis** - "Run DE analysis", "Start downstream analysis"

## Troubleshooting

### Ollama Not Running

```bash
# Check if Ollama is running
systemctl status ollama

# Start Ollama manually
ollama serve

# Check logs
journalctl -u ollama -f
```

### Model Not Found

```bash
# List installed models
ollama list

# Pull missing model
ollama pull llama3.1:8b
```

### Connection Refused

```bash
# Check Ollama is listening
netstat -tuln | grep 11434

# Test connection
curl http://localhost:11434/api/tags
```

### Out of Memory

```bash
# Use smaller model
ollama pull mistral:7b

# Or use quantized version
ollama pull llama3.1:8b-q4_0
```

### Slow Responses

**Tips to improve speed:**
- Use smaller models (mistral:7b instead of llama3.1:70b)
- Ensure Ollama has GPU access (if available)
- Reduce context window: Edit `_chat_ollama()` and set `"num_ctx": 2048`

## Performance Comparison

| Model | Size | RAM/VRAM | Speed | Quality | Tool Calling | Recommended For |
|-------|------|----------|-------|---------|--------------|-----------------|
| mistral:7b | 4GB | 8GB RAM | Fast | Good | Text pattern | Quick status checks |
| llama3.1:8b | 5GB | 8GB RAM | Medium | Great | Native | GPU 없는 환경 |
| **qwen2.5:32b** | 20GB | 20GB VRAM | Medium | **Excellent** | **Native (GPU)** | **기본 권장** |
| llama3.1:70b | 40GB | 42GB RAM | Slow | Best | Native | 고성능 서버 |

## Security Notes

✅ **Data stays local** - All model inference runs on your server  
✅ **No external API calls** - No data sent to cloud providers  
✅ **Full audit trail** - All queries logged locally  
✅ **Compliance-friendly** - Meets HIPAA/GDPR requirements for data handling  

## Advanced Configuration

### Running Ollama as System Service

```bash
# Create systemd service
sudo systemctl enable ollama
sudo systemctl start ollama

# Check status
sudo systemctl status ollama
```

### GPU Acceleration

If you have an NVIDIA GPU:

```bash
# Ollama will automatically use GPU if available
nvidia-smi  # Verify GPU is detected

# Monitor GPU usage during inference
watch -n 1 nvidia-smi
```

### Model Customization

Create custom Modelfile for fine-tuned behavior:

```bash
# Create Modelfile
cat > Modelfile << 'EOF'
FROM llama3.1:8b

SYSTEM """
You are an expert bioinformatics assistant specializing in RNA-seq analysis.
You focus on clear, actionable insights for genomics researchers.
"""

PARAMETER temperature 0.5
PARAMETER num_ctx 4096
EOF

# Create custom model
ollama create genomics-assistant -f Modelfile

# Use custom model
python scripts/standardization/llm_agent.py \
  --model genomics-assistant \
  --interactive
```

## Next Steps

1. **Test with your project**:
   ```bash
   python scripts/standardization/llm_agent.py \
     --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
     --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
     --interactive
   ```

2. **Try different queries** to explore capabilities

3. **Integrate into workflow** - Use in Snakemake or shell scripts:
   ```bash
   # Example: Automated QC check
   response=$(python scripts/standardization/llm_agent.py \
     --project-summary project_summary.json \
     --rnaseq-output /output \
     --message "Are all samples ready for DE analysis?")
   
   echo "$response"
   ```

4. **Provide feedback** - Help improve the agent by reporting issues

## Resources

- Ollama Documentation: https://ollama.com/
- Model Library: https://ollama.com/library
- Python Client: https://github.com/ollama/ollama-python
- Project Issues: https://github.com/ibs-CMG-NGS/rna-seq-pipeline/issues

## Support

For questions or issues:
1. Check troubleshooting section above
2. Review agent logs in terminal output
3. Test Ollama directly: `ollama run llama3.1:8b "test"`
4. Open issue on GitHub with error details
