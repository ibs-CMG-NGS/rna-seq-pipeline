# English Documentation

This folder contains English versions of the user documentation.

| File | Description |
|------|-------------|
| [LLM_AGENT_QUICKSTART_EN.md](LLM_AGENT_QUICKSTART_EN.md) | Quick start guide for the LLM pipeline agent |
| [OLLAMA_SETUP_EN.md](OLLAMA_SETUP_EN.md) | Ollama installation and model setup |
| [PIPELINE_GUIDE_EN.md](PIPELINE_GUIDE_EN.md) | Full pipeline guide (manual + agent) |

## Korean documentation

The original Korean documentation is in `docs/user/`:
- `LLM_AGENT_QUICKSTART.md`
- `OLLAMA_SETUP.md`
- `PIPELINE_GUIDE.md`

## Quick Start (TL;DR)

```bash
# 1. Start the agent
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline
python scripts/standardization/llm_agent.py --interactive --model qwen2.5:32b

# 2. Tell it what you want (in English or Korean)
You: Find FASTQ files in /data_3tb/shared/chd8-rna-seq-raw-data/fastq/
You: Create a project called my-project, save to /data_3tb/shared/output/my-project/
You: Run the pipeline with 16 cores
```
