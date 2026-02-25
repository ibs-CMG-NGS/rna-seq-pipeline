# Phase 5: Ollama Installation & Setup

**Status**: In Progress  
**Created**: 2026-02-26  
**Objective**: Install Ollama and prepare LLM agent for natural language queries

---

## Overview

Ollama is a local LLM runtime that allows us to run models like Llama 3.1 without external API dependencies. This phase sets up Ollama on the server to enable natural language interaction with the RNA-seq pipeline.

**Benefits:**
- 🔒 Privacy: No data sent to external services
- 🚀 Speed: Local inference, no network latency
- 💰 Cost: No API fees
- 🛡️ Security: Complete control over the environment

---

## Prerequisites

**System Requirements:**
- Ubuntu 20.04+ (or other Linux distribution)
- 8GB+ RAM (16GB recommended for llama3.1:8b)
- 10GB+ free disk space
- GPU (optional, but recommended for faster inference)

**Server Info:**
```bash
# Check current system
uname -a
free -h
df -h /home
```

---

## Installation Steps

### Step 1: Install Ollama

**Option A: Official Script (Recommended)**
```bash
# Download and run installer
curl -fsSL https://ollama.com/install.sh | sh
```

**Option B: Manual Download**
```bash
# Download binary
curl -L https://ollama.com/download/ollama-linux-amd64 -o ollama
chmod +x ollama
sudo mv ollama /usr/local/bin/

# Create service user
sudo useradd -r -s /bin/false -m -d /usr/share/ollama ollama
```

**Verify Installation:**
```bash
ollama --version
```

Expected output: `ollama version is X.X.X`

---

### Step 2: Start Ollama Service

**Interactive Mode (Testing):**
```bash
# Start server in foreground
ollama serve
```

**Background Mode (Production):**
```bash
# Using systemd (if installed via script)
sudo systemctl start ollama
sudo systemctl enable ollama
sudo systemctl status ollama

# Check logs
sudo journalctl -u ollama -f
```

**Verify Server Running:**
```bash
# Check if API is responding
curl http://localhost:11434/api/version
```

Expected output: JSON with version info

---

### Step 3: Pull Llama 3.1 Model

**Download Model:**
```bash
# Pull 8B parameter model (~4.7GB)
ollama pull llama3.1:8b
```

**Alternative Models:**
```bash
# Smaller model (faster, less accurate)
ollama pull llama3.1:7b

# Larger model (slower, more accurate)
ollama pull llama3.1:70b  # Requires 64GB+ RAM
```

**Verify Model:**
```bash
# List downloaded models
ollama list
```

Expected output:
```
NAME              ID              SIZE      MODIFIED
llama3.1:8b       a1b2c3d4...     4.7 GB    X hours ago
```

---

### Step 4: Test Ollama

**Simple Query Test:**
```bash
# Test chat interface
ollama run llama3.1:8b "Hello, how are you?"
```

**API Test:**
```bash
# Test via API
curl -X POST http://localhost:11434/api/generate \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.1:8b",
    "prompt": "What is 2+2?",
    "stream": false
  }'
```

**Python Test:**
```python
import requests
import json

response = requests.post(
    'http://localhost:11434/api/generate',
    json={
        'model': 'llama3.1:8b',
        'prompt': 'Say hello',
        'stream': False
    }
)

print(response.json()['response'])
```

---

## Configuration

### Ollama Environment Variables

Create `/etc/ollama/ollama.env` (if using systemd):

```bash
# API settings
OLLAMA_HOST=0.0.0.0:11434  # Listen on all interfaces
OLLAMA_ORIGINS=*            # CORS origins

# Resource limits
OLLAMA_NUM_PARALLEL=1       # Number of parallel requests
OLLAMA_MAX_LOADED_MODELS=1  # Max models in memory

# GPU settings (if available)
OLLAMA_CUDA_VISIBLE_DEVICES=0  # Use first GPU
```

**Apply Configuration:**
```bash
sudo systemctl restart ollama
```

---

### Resource Limits

**For systemd service**, edit `/etc/systemd/system/ollama.service`:

```ini
[Service]
# Memory limit
MemoryMax=12G
MemoryHigh=10G

# CPU limit
CPUQuota=400%  # 4 cores

# Process limit
LimitNOFILE=65536
```

**Apply Limits:**
```bash
sudo systemctl daemon-reload
sudo systemctl restart ollama
```

---

## Integration with RNA-seq Pipeline

### Step 1: Install Python Dependencies

```bash
cd /data_3tb/shared/rna-seq-pipeline
conda activate rna-seq-pipeline

# Install Ollama Python client (if needed)
pip install ollama

# Or use requests (already available)
# No additional dependencies needed
```

---

### Step 2: Test LLM Agent

**Check Agent Configuration:**
```bash
# View agent help
python scripts/standardization/llm_agent.py --help
```

**Run Agent:**
```bash
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8
```

**Expected Startup:**
```
🤖 RNA-seq Pipeline LLM Agent
📊 Project: mouse-chd8
📁 Output: /data_3tb/shared/output/mouse-chd8
🔧 Available tools: 7

Type 'help' for commands, 'quit' to exit.
You> 
```

---

### Step 3: Test Natural Language Queries

**Query 1: Sample Information**
```
You> 현재 프로젝트의 QC 상태 보여줘
```

Expected: Agent calls `get_sample_details()` → displays QC summary

**Query 2: Path Validation**
```
You> 경로 검증해줘
```

Expected: Agent calls `validate_paths()` → shows validation results

**Query 3: DE Preparation**
```
You> mouse-chd8 프로젝트 DE 분석 준비해줘
```

Expected: Agent calls `prepare_de_analysis()` → runs bridge script

**Query 4: Config Check**
```
You> bridge 설정 파일 확인해줘
```

Expected: Agent calls `check_bridge_config()` → shows config status

---

## Troubleshooting

### Issue: Ollama Not Starting

**Symptom:**
```
Error: could not connect to ollama server
```

**Solutions:**

1. **Check if service is running:**
```bash
sudo systemctl status ollama
ps aux | grep ollama
```

2. **Check port availability:**
```bash
sudo lsof -i :11434
netstat -tlnp | grep 11434
```

3. **Check logs:**
```bash
sudo journalctl -u ollama -n 100 --no-pager
```

4. **Manual start for debugging:**
```bash
OLLAMA_DEBUG=1 ollama serve
```

---

### Issue: Model Not Found

**Symptom:**
```
Error: model 'llama3.1:8b' not found
```

**Solutions:**

1. **List available models:**
```bash
ollama list
```

2. **Pull model again:**
```bash
ollama pull llama3.1:8b
```

3. **Check disk space:**
```bash
df -h /usr/share/ollama
```

---

### Issue: Out of Memory

**Symptom:**
```
Error: failed to allocate memory
```

**Solutions:**

1. **Check available RAM:**
```bash
free -h
```

2. **Use smaller model:**
```bash
ollama pull llama3.1:7b
```

3. **Increase swap:**
```bash
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

4. **Limit model memory:**
```bash
# Set in ollama.env
OLLAMA_MAX_LOADED_MODELS=1
```

---

### Issue: Slow Inference

**Symptom:**
Responses take 30+ seconds

**Solutions:**

1. **Check CPU usage:**
```bash
top
htop
```

2. **Use GPU (if available):**
```bash
# Check GPU
nvidia-smi

# Configure Ollama to use GPU
OLLAMA_CUDA_VISIBLE_DEVICES=0 ollama serve
```

3. **Use quantized model:**
```bash
# 4-bit quantization (smaller, faster)
ollama pull llama3.1:8b-q4_0
```

4. **Reduce context size:**
```python
# In agent, limit prompt length
max_tokens = 2048  # Instead of 4096
```

---

### Issue: LLM Agent Not Understanding Queries

**Symptom:**
Agent doesn't call correct tools

**Solutions:**

1. **Check system prompt:**
```bash
# View agent's tool definitions
grep -A 10 "Available tools" scripts/standardization/llm_agent.py
```

2. **Use more explicit queries:**
```
Bad:  "정보 보여줘"
Good: "mouse-chd8 프로젝트의 샘플 QC 상태 보여줘"
```

3. **Check model temperature:**
```python
# Lower temperature for more deterministic responses
temperature = 0.1  # Instead of 0.7
```

4. **Add examples to prompt:**
```python
# In system prompt, add few-shot examples
examples = """
Q: "QC 상태 보여줘" → use get_sample_details()
Q: "DE 분석 준비" → use prepare_de_analysis()
"""
```

---

## Security Considerations

### Network Access

**Restrict to Localhost (Recommended):**
```bash
# In ollama.env
OLLAMA_HOST=127.0.0.1:11434
```

**Allow Remote Access (Use with Caution):**
```bash
# In ollama.env
OLLAMA_HOST=0.0.0.0:11434

# Use firewall to restrict
sudo ufw allow from 192.168.1.0/24 to any port 11434
```

---

### Resource Isolation

**Using systemd (Recommended):**
```ini
[Service]
# Run as dedicated user
User=ollama
Group=ollama

# Restrict filesystem access
PrivateTmp=yes
ProtectSystem=strict
ProtectHome=yes
ReadWritePaths=/usr/share/ollama

# Restrict capabilities
NoNewPrivileges=yes
CapabilityBoundingSet=
```

**Using Docker (Alternative):**
```bash
docker run -d \
  --name ollama \
  --gpus all \
  -v ollama:/root/.ollama \
  -p 127.0.0.1:11434:11434 \
  --memory=12g \
  --cpus=4 \
  ollama/ollama
```

---

### Audit Logging

**Log All LLM Queries:**

Create `/data_3tb/shared/rna-seq-pipeline/scripts/utils/audit_log.py`:

```python
import logging
from datetime import datetime

def setup_audit_log():
    logging.basicConfig(
        filename='/var/log/rna-seq-pipeline/llm_audit.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def log_query(user, query, tool_called, result):
    logging.info(
        f"User={user} | Query={query} | Tool={tool_called} | Result={result}"
    )
```

---

## Performance Benchmarks

### Expected Response Times

| Query Type | Model Size | Hardware | Response Time |
|------------|------------|----------|---------------|
| Simple info | 8B | CPU (8 cores) | 5-10s |
| Tool call | 8B | CPU (8 cores) | 8-15s |
| Complex reasoning | 8B | CPU (8 cores) | 15-30s |
| Simple info | 8B | GPU (RTX 3090) | 1-2s |
| Tool call | 8B | GPU (RTX 3090) | 2-4s |
| Complex reasoning | 8B | GPU (RTX 3090) | 4-8s |

**Optimization Tips:**
- Use GPU for 5-10x speedup
- Use quantized models for 2x speedup
- Cache frequent queries
- Batch requests when possible

---

## Next Steps

After completing Ollama setup:

1. ✅ **Verify Installation**: All tests passing
2. ✅ **Test Simple Queries**: Basic chat working
3. 🔄 **Test Tool Integration**: Agent calling tools correctly
4. ⏳ **Production Deployment**: systemd service configured
5. ⏳ **Monitoring Setup**: Logs and metrics collection

---

## References

- **Ollama Documentation**: https://ollama.com/docs
- **Llama 3.1 Model Card**: https://ollama.com/library/llama3.1
- **LLM Agent Implementation**: `scripts/standardization/llm_agent.py`
- **Security Guide**: `docs/developer/SECURITY_GUIDE.md`
- **E2E Testing**: `docs/developer/PHASE4_E2E_TESTING.md`
