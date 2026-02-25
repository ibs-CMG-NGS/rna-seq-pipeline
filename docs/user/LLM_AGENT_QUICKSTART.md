# LLM Agent Quick Reference

## Quick Start

```bash
# 1. Install Ollama
curl -fsSL https://ollama.com/install.sh | sh
ollama serve  # In separate terminal

# 2. Download model
ollama pull llama3.1:8b

# 3. Install Python client
conda activate rna-seq-pipeline
pip install ollama

# 4. Start interactive session
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --interactive
```

## Example Queries

### QC Status
```
You: Show me the QC status
You: How many samples passed QC?
You: What's the overall pass rate?
```

### Condition Comparison
```
You: Compare wildtype and heterozygous
You: Show me the differences between conditions
You: Which condition has better quality?
```

### Failed Samples
```
You: Which samples failed?
You: Show me QC failures
You: What issues were found?
```

### Sample Details
```
You: Tell me about sample Ctrl_3
You: What's the status of Ctrl_3?
You: Show details for Ctrl_3
```

### Conditions List
```
You: What conditions do I have?
You: List all experimental groups
You: Show me condition summary
```

### Start DE Analysis
```
You: Run DE analysis
You: Start downstream analysis
You: Can you prepare DE/GO pipeline?
```

## Command Line Options

```bash
# Interactive mode (recommended)
--interactive

# Single query
--message "Your question here"

# Select model
--model llama3.1:8b       # Default
--model llama3.1:70b      # Higher quality, slower
--model mistral:7b        # Faster

# Custom Ollama server
--ollama-host http://gpu-server:11434

# Required paths
--project-summary /path/to/project_summary.json
--rnaseq-output /path/to/output/directory
--de-pipeline /path/to/de-pipeline  # Optional
```

## Model Selection

| Model | Speed | Quality | RAM | Use Case |
|-------|-------|---------|-----|----------|
| mistral:7b | ⚡⚡⚡ | ⭐⭐⭐ | 8GB | Quick queries |
| llama3.1:8b | ⚡⚡ | ⭐⭐⭐⭐ | 8GB | General (default) |
| llama3.1:70b | ⚡ | ⭐⭐⭐⭐⭐ | 48GB | Complex analysis |

## Troubleshooting

### Connection Refused
```bash
# Start Ollama server
ollama serve
```

### Model Not Found
```bash
# Download model
ollama pull llama3.1:8b

# List installed models
ollama list
```

### Slow Responses
```bash
# Use faster model
--model mistral:7b

# Or reduce context window (edit llm_agent.py)
"num_ctx": 2048
```

## Data Security

✅ All data stays on your server  
✅ No external API calls  
✅ HIPAA/GDPR compliant  
✅ Full audit control  

## Full Documentation

- User Guide: `docs/user/OLLAMA_SETUP.md`
- Developer Guide: `docs/developer/PHASE6_LLM_INTEGRATION.md`
- API Documentation: https://ollama.com/
