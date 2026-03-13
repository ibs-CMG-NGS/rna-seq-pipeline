#!/usr/bin/env python3
"""
LLM-powered Natural Language Agent for RNA-seq Pipeline Management

Enables conversational interaction with the pipeline:
- "What's the QC status of my samples?"
- "Compare wildtype vs heterozygous"
- "Start DE analysis"

Security: Input validation and sandboxing applied to all tool executions.
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
import subprocess
import sys

# Import security utilities
# Support running from project root OR from scripts/standardization/
_project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_project_root))
try:
    from scripts.utils.security import (
        validate_project_id,
        validate_sample_id,
        validate_path,
        sanitize_command_args,
        SecurityError
    )
    SECURITY_ENABLED = True
except ImportError:
    SECURITY_ENABLED = False

# LLM Integration
try:
    import openai  # OpenAI GPT-4 (cloud)
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False

try:
    import anthropic  # Anthropic Claude (cloud)
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False

try:
    import ollama  # Ollama (local) - RECOMMENDED for data security
    HAS_OLLAMA = True
except ImportError:
    HAS_OLLAMA = False

try:
    from llama_cpp import Llama  # llama.cpp (local)
    HAS_LLAMACPP = True
except ImportError:
    HAS_LLAMACPP = False


class PipelineAgent:
    """Natural language agent for pipeline management."""
    
    def __init__(self,
                 project_summary_path: Path = None,
                 rnaseq_output_dir: Path = None,
                 de_pipeline_dir: Path = None,
                 llm_provider: str = "ollama",
                 api_key: str = None,
                 model: str = None,
                 ollama_host: str = "http://localhost:11434",
                 initial_config: str = None):
        """
        Initialize agent.
        
        Args:
            project_summary_path: Path to project_summary.json (optional — not needed for pipeline execution tools)
            rnaseq_output_dir: RNA-seq pipeline output directory (optional)
            de_pipeline_dir: DE/GO pipeline directory (optional)
            llm_provider: "ollama" (local, recommended), "openai", "anthropic", "llamacpp"
            api_key: API key for cloud LLM providers
            model: Model name (e.g., "llama3.1", "gpt-4", "claude-3-sonnet")
            ollama_host: Ollama server URL (default: http://localhost:11434)
        """
        self.project_summary_path = Path(project_summary_path) if project_summary_path else None
        self.rnaseq_output_dir = Path(rnaseq_output_dir) if rnaseq_output_dir else None
        self.de_pipeline_dir = Path(de_pipeline_dir) if de_pipeline_dir else None
        self.llm_provider = llm_provider
        self.ollama_host = ollama_host
        
        # Load project summary (optional — pipeline execution tools don't require it)
        if self.project_summary_path and self.project_summary_path.exists():
            with open(self.project_summary_path) as f:
                self.project_summary = json.load(f)
            self.project_id = self.project_summary['project_id']
        else:
            self.project_summary = {}
            self.project_id = "unknown"
            if self.project_summary_path:
                print(f"⚠️  project_summary.json not found: {self.project_summary_path}")
                print("   Pipeline execution tools (detect_fastq, run_pipeline, etc.) are still available.")
        
        # Initialize LLM client
        if llm_provider == "ollama" and HAS_OLLAMA:
            self.llm_client = ollama
            self.model = model or "qwen3:14b"  # Default to Qwen3 14B
            print(f"Using local Ollama model: {self.model}")
            
        elif llm_provider == "llamacpp" and HAS_LLAMACPP:
            if not model:
                raise ValueError("Model path required for llama.cpp")
            self.llm_client = Llama(model_path=model)
            self.model = model
            print(f"Using llama.cpp with model: {model}")
            
        elif llm_provider == "openai" and HAS_OPENAI:
            openai.api_key = api_key
            self.llm_client = openai
            self.model = model or "gpt-4"
            print(f"Using OpenAI: {self.model}")
            
        elif llm_provider == "anthropic" and HAS_ANTHROPIC:
            self.llm_client = anthropic.Anthropic(api_key=api_key)
            self.model = model or "claude-3-sonnet-20240229"
            print(f"Using Anthropic: {self.model}")
            
        else:
            raise ValueError(f"LLM provider {llm_provider} not available or not installed. "
                           f"Available: ollama={HAS_OLLAMA}, openai={HAS_OPENAI}, "
                           f"anthropic={HAS_ANTHROPIC}, llamacpp={HAS_LLAMACPP}")
        
        # Define available tools
        self.tools = self._define_tools()

        # Conversation history
        self.conversation_history = []

        # Session state: remember paths used in this session
        self.current_config: Optional[str] = None
        self.current_de_pipeline_dir: Optional[str] = None
        self.current_de_config: Optional[str] = None   # DE pipeline config YAML
        self.current_data_dir: Optional[str] = None
        self.current_project_id: Optional[str] = None
        self.current_results_dir: Optional[str] = None

        # Persistent session file — survives agent restarts
        self.session_file = _project_root / ".agent_session.json"
        self._load_session()

        # Explicit config overrides session file
        if initial_config:
            self.current_config = initial_config
            self._save_session()
    
    def _define_tools(self) -> List[Dict]:
        """Define available tools for LLM function calling."""
        return [
            {
                "name": "switch_project",
                "description": (
                    "현재 세션의 프로젝트를 초기화하고 새 프로젝트로 전환. "
                    "기존 config/data_dir/results_dir 등 세션 상태를 모두 지운다. "
                    "트리거: '새 프로젝트 시작', '다른 프로젝트', '프로젝트 바꿔', '새로 시작', "
                    "'새 분석 시작', '프로젝트 전환', 'new project', 'switch project'. "
                    "새 config 경로를 알고 있으면 new_config에 전달하고, 모르면 비워두면 됨."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "new_config": {
                            "type": "string",
                            "description": "새 프로젝트 config 경로 (알고 있는 경우). 없으면 생략."
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "get_project_status",
                "description": "파이프라인 완료 후 샘플 QC 결과 요약 (맵핑률, 통과/실패 샘플 수 등). 'QC 상태 보여줘', 'QC 결과 어때?', '몇 개 통과했어?', 'show QC summary' → 이 도구 사용. 실행 중 진행률은 monitor_pipeline 사용.",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "compare_conditions",
                "description": "RNA-seq 파이프라인 QC 지표(맵핑률, assignment rate 등)를 조건 간 비교. ⚠️ DE 분석 결과(유전자 발현 변화) 요약에는 절대 사용하지 말 것 — 그것은 read_de_results_summary 사용.",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "get_failed_samples",
                "description": "Get list of samples that failed QC",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "get_sample_details",
                "description": "Get detailed information about a specific sample",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "sample_id": {
                            "type": "string",
                            "description": "Sample identifier"
                        }
                    },
                    "required": ["sample_id"]
                }
            },
            {
                "name": "list_conditions",
                "description": "List all experimental conditions in the project",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "start_de_analysis",
                "description": "Start differential expression and GO enrichment analysis",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "confirm": {
                            "type": "boolean",
                            "description": "User confirmation to proceed"
                        }
                    },
                    "required": ["confirm"]
                }
            },
            {
                "name": "prepare_de_analysis",
                "description": "Prepare DE/GO analysis by copying counts, generating metadata, and creating config. Does not start the actual analysis.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {
                            "type": "string",
                            "description": "Project identifier (optional, uses current project if not specified)"
                        },
                        "exclude_samples": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": "Sample IDs to exclude from counts matrix and metadata (e.g. failed QC samples)."
                        },
                        "config_file": {
                            "type": "string",
                            "description": "Path to the RNA-seq project config YAML (e.g. config/projects/mouse-monSTIM-2026.yaml). The bridge script will read rnaseq_output, project_id, and species directly from this file — do NOT guess rnaseq_output separately when this is provided."
                        },
                        "rnaseq_output": {
                            "type": "string",
                            "description": "RNA-seq pipeline output directory. Only set this if config_file is NOT provided."
                        },
                        "de_pipeline_dir": {
                            "type": "string",
                            "description": "DE/GO analysis pipeline directory (overrides session/config auto-detection)"
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "check_bridge_config",
                "description": "Check if bridge configuration exists for the project, create if missing",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {
                            "type": "string",
                            "description": "Project identifier"
                        },
                        "force_regenerate": {
                            "type": "boolean",
                            "description": "Force regenerate config even if exists"
                        }
                    },
                    "required": ["project_id"]
                }
            },
            {
                "name": "validate_paths",
                "description": "Validate that all required paths exist for DE/GO analysis",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {
                            "type": "string",
                            "description": "Project identifier"
                        }
                    },
                    "required": ["project_id"]
                }
            },
            # Multi-axis analysis tools
            {
                "name": "get_sample_axes",
                "description": "Get all experimental axes (genotype, tissue, sex) and sample counts per group",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "compare_by_axis",
                "description": "Compare QC metrics grouped by a specific axis (e.g. condition, tissue, sex), with optional pre-filtering by other axes. Use this instead of compare_conditions when the user wants to filter by tissue, sex, or any other axis.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "axis": {
                            "type": "string",
                            "description": "Axis to group comparison by. Any samplesheet column is valid: 'condition', 'tissue', 'sex', 'age', etc."
                        },
                        "filters": {
                            "type": "object",
                            "description": "Optional filters to apply before comparison. e.g. {\"tissue\": \"Hippocampus\"} or {\"condition\": \"wildtype\"}",
                            "additionalProperties": {"type": "string"}
                        }
                    },
                    "required": ["axis"]
                }
            },
            {
                "name": "filter_samples",
                "description": "Get list of samples matching multiple axis filters simultaneously",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "filters": {
                            "type": "object",
                            "description": "Axis filters e.g. {\"tissue\": \"HPC\", \"sex\": \"Male\", \"genotype\": \"wildtype\"}",
                            "additionalProperties": {"type": "string"}
                        }
                    },
                    "required": ["filters"]
                }
            },
            # Phase 8A: Pipeline Execution Tools
            {
                "name": "create_project_config",
                "description": "Create a complete project configuration YAML file with all required parameters. Always ask the user for genome paths (genome_fasta, annotation_gtf, star_index) before calling this tool.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id":     {"type": "string", "description": "Project identifier (e.g. 'mouse-chd8')"},
                        "data_dir":       {"type": "string", "description": "Absolute path to directory containing FASTQ files"},
                        "results_dir":    {"type": "string", "description": "Absolute path to base output directory (project subdir created inside)"},
                        "species":        {"type": "string", "description": "Species: 'human', 'mouse', or 'rat'", "default": "human"},
                        "genome_dir":     {"type": "string", "description": "Absolute path to genome reference directory"},
                        "genome_fasta":   {"type": "string", "description": "Absolute path to genome FASTA file"},
                        "annotation_gtf": {"type": "string", "description": "Absolute path to gene annotation GTF file"},
                        "star_index":     {"type": "string", "description": "Absolute path to STAR genome index directory"},
                        "genome_build":   {"type": "string", "description": "Genome build string (e.g. 'GRCh38', 'GRCm38')"},
                        "use_sample_sheet": {"type": "boolean", "description": "Use sample sheet TSV instead of auto-detecting FASTQs", "default": False},
                        "sample_sheet":   {"type": "string", "description": "Path to sample sheet TSV (if use_sample_sheet=true)"},
                        "threads":        {"type": "integer", "description": "CPU threads for alignment", "default": 12},
                        "memory_gb":      {"type": "integer", "description": "RAM limit in GB", "default": 48},
                        "strandedness":   {"type": "integer", "description": "0=unstranded, 1=forward, 2=reverse", "default": 0}
                    },
                    "required": ["project_id", "data_dir", "results_dir"]
                }
            },
            {
                "name": "detect_fastq_files",
                "description": "Scan directory and detect FASTQ files, determine read type and sample count",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "data_dir": {"type": "string", "description": "Directory to scan for FASTQ files"}
                    },
                    "required": ["data_dir"]
                }
            },
            {
                "name": "validate_input_data",
                "description": "Pre-flight validation before pipeline execution (disk space, reference genome, etc.)",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml file"}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "run_pipeline",
                "description": "Execute RNA-seq pipeline with Snakemake. Use dry_run=false to ACTUALLY RUN the pipeline. Use dry_run=true only for preview. When user says '실행해줘', '돌려줘', '시작해줘', 'run', 'execute', 'start' → set dry_run=false. Real runs (dry_run=false) launch in background by default so the agent stays responsive — use monitor_pipeline to check progress.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "cores": {"type": "integer", "description": "Number of cores (e.g. 16)", "default": 8},
                        "dry_run": {"type": "boolean", "description": "false=actually run pipeline, true=preview only. Default false when user says run/execute/start.", "default": False},
                        "background": {"type": "boolean", "description": "true=launch in background and return immediately (default for real runs); false=block until complete (used by batch runner). Always true when dry_run=false.", "default": True}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "validate_library_type",
                "description": (
                    "RSeQC를 이용한 pre-flight 라이브러리 검증. "
                    "FASTQ 서브셋을 STAR로 정렬 후 infer_experiment.py(strandedness 자동 감지)와 "
                    "read_distribution.py(mRNA-seq 여부 확인)를 실행. "
                    "strandedness가 config와 다르면 자동 업데이트. "
                    "메인 파이프라인 실행 전에 호출할 것. "
                    "트리거: '라이브러리 검증', 'strandedness 확인', 'mRNA-seq 맞는지', 'library type'"
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "sample_id": {"type": "string", "description": "Sample ID to use for validation (default: first sample in sheet)"},
                        "n_reads": {"type": "integer", "description": "Number of reads to subset for validation (default 500000)", "default": 500000},
                        "cores": {"type": "integer", "description": "Threads for STAR alignment (default 8)", "default": 8}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "setup_and_validate",
                "description": (
                    "FASTQ 감지 + 샘플시트 생성 + 라이브러리 검증을 한 번에 실행. "
                    "새 프로젝트 시작 시 가장 먼저 호출할 것. "
                    "트리거: '프로젝트 설정해줘', '한번에 설정해줘', 'setup', '설정부터 해줘', "
                    "'처음부터 설정', '자동 설정'"
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "sample_id": {"type": "string", "description": "Sample ID for library validation (default: first sample)"},
                        "n_reads": {"type": "integer", "description": "Reads to subset for validation (default 500000)", "default": 500000},
                        "cores": {"type": "integer", "description": "Threads for STAR (default 8)", "default": 8}
                    },
                    "required": ["config_file"]
                }
            },
            # Phase 8B tools
            {
                "name": "monitor_pipeline",
                "description": "Snakemake 파이프라인이 현재 실행 중일 때 실시간 진행률 확인 (완료된 job 수, 진행 중 rule 등). '파이프라인 진행상황 알려줘', '몇 개 완료됐어?', '얼마나 됐어?', '지금 뭐 하고 있어?', 'pipeline progress' → 이 도구 사용. QC 결과/분석 완료 여부는 get_project_status 사용.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml — project_id and results_dir are read from it automatically"}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "create_sample_sheet",
                "description": "Generate a sample metadata TSV from detected FASTQ files, optionally assigning conditions. If the TSV already exists and overwrite is not explicitly requested, preserves existing file to protect manually revised conditions.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {"type": "string", "description": "Project identifier"},
                        "data_dir": {"type": "string", "description": "Directory containing FASTQ files"},
                        "conditions": {
                            "type": "object",
                            "description": "Optional condition mapping: {condition_name: [sample_id_substrings]}",
                            "additionalProperties": {"type": "array", "items": {"type": "string"}}
                        },
                        "output_path": {"type": "string", "description": "Output TSV path (optional)"},
                        "overwrite": {"type": "boolean", "description": "If false (default), skip regeneration when TSV already exists to preserve manual condition edits. Set true only when user explicitly wants to regenerate.", "default": False}
                    },
                    "required": ["project_id", "data_dir"]
                }
            },
            {
                "name": "estimate_resources",
                "description": "Estimate runtime, disk and RAM requirements before pipeline execution",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "cores": {"type": "integer", "description": "Number of cores to use", "default": 8}
                    },
                    "required": ["config_file"]
                }
            },
            # DE analysis bridge
            {
                "name": "run_bridge",
                "description": (
                    "Transfer RNA-seq counts matrix to DE/GO analysis pipeline. "
                    "Use when user asks to: start DE analysis, connect to DE pipeline, "
                    "'DE 분석 연결해줘', 'DE 파이프라인으로 넘겨줘', 'counts 넘겨줘'. "
                    "By default skip_de=true — only prepares files (counts copy + metadata + DE config). "
                    "Set skip_de=false only if user explicitly asks to RUN DE analysis immediately. "
                    "Use exclude_samples to drop failed/unwanted samples from counts and metadata."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {
                            "type": "string",
                            "description": "Path to RNA-seq project config YAML"
                        },
                        "de_pipeline_dir": {
                            "type": "string",
                            "description": "Absolute path to DE/GO analysis pipeline directory"
                        },
                        "skip_de": {
                            "type": "boolean",
                            "description": "If true (default), only prepare files. If false, also trigger DE Snakemake."
                        },
                        "dry_run": {
                            "type": "boolean",
                            "description": "If true, show what would be done without writing files."
                        },
                        "exclude_samples": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": "Sample IDs to exclude from counts matrix and metadata (e.g. failed QC samples)."
                        }
                    },
                    "required": ["config_file", "de_pipeline_dir"]
                }
            },
            # Phase 8C: Project/Pipeline information reading tools
            {
                "name": "read_project_config",
                "description": "Read and summarise the project config file. Use when user asks about project setup, paths, species, genome build, or any config settings.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "read_sample_sheet",
                "description": "Read the sample sheet to list all samples, conditions, tissues, replicates and FASTQ paths. Use when user asks about samples, conditions, or group breakdown.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml (sample_sheet path is read from it)"}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "read_qc_results",
                "description": "Read QC results: FastQC evaluation, MultiQC stats, STAR mapping rates. Use when user asks about QC status, mapping rates, failed samples, or data quality.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "read_counts",
                "description": "Read the featureCounts expression matrix. Use when user asks about gene expression, counts, specific gene levels, or top expressed genes.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "genes": {"type": "array", "items": {"type": "string"}, "description": "Specific gene IDs/names to look up (optional)"},
                        "top_n": {"type": "integer", "description": "Return top N most variable genes if genes not specified", "default": 20}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "read_pipeline_logs",
                "description": "Read pipeline execution logs. Use when user asks about errors, warnings, a specific sample's log, or pipeline execution details.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "rule": {"type": "string", "description": "Specific rule: cutadapt / star / fastqc / featurecounts (optional)"},
                        "sample_id": {"type": "string", "description": "Specific sample ID (optional)"},
                        "tail_lines": {"type": "integer", "description": "Number of tail lines per log", "default": 50}
                    },
                    "required": ["config_file"]
                }
            },
            {
                "name": "read_de_results_summary",
                "description": (
                    "⭐ DESeq2 DE 분석 결과 요약 — 유의미한 유전자 수(up/down), 상위 유전자 목록, GO BP 경로 반환. "
                    "de_config_file은 생략 가능 (세션의 current_de_config 자동 사용). "
                    "⚠️ 'DE 결과', '발현 변화', '유전자 변했어', 'D1 결과', '차발현' 관련 질문은 반드시 이 도구 사용. "
                    "compare_conditions / get_project_status / read_counts 절대 사용 금지."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "DE pipeline config YAML 경로 (생략 시 세션 current_de_config 사용)"
                        },
                        "pair": {
                            "type": "string",
                            "description": "특정 비교 pair (예: 'D1_vs_Control'). 생략 시 전체 pair 요약."
                        },
                        "top_n": {
                            "type": "integer",
                            "description": "상위 유전자 반환 수 (기본 10)",
                            "default": 10
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "check_analysis_readiness",
                "description": (
                    "DE 분석 실행 전 종합 진단 — counts/metadata 파일 존재, 샘플 ID 일치, "
                    "조건명 일치, 조건별 n 수, counts 통계를 한 번에 점검. "
                    "de_config_file 생략 가능 (세션 current_de_config 자동 사용). "
                    "트리거: 'DE 실행해도 돼?', '분석 준비됐어?', 'DE 전에 확인해줘', '샘플 수 맞아?'"
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "DE pipeline config YAML 경로 (생략 가능)"
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "validate_de_config_conditions",
                "description": (
                    "DE config의 pairwise_comparisons 조건명이 metadata CSV의 실제 condition 값과 일치하는지 검증. "
                    "불일치 감지 + 자동 수정 제안 (대소문자 차이, 1D↔D1 순서 스왑 등). "
                    "run_bridge 완료 후 반드시 호출. "
                    "트리거: 'DE config 검증해줘', '조건 확인해줘', run_bridge 후 자동 검증."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "Path to the DE pipeline config YAML (e.g., configs/config_PROJECT.yml)"
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "apply_de_config_corrections",
                "description": (
                    "validate_de_config_conditions에서 발견된 조건명 불일치를 DE config에 자동 적용. "
                    "반드시 validate_de_config_conditions 실행 후 suggestions 확인 후 호출. "
                    "트리거: '수정해줘', '자동으로 고쳐줘', validate 후 사용자 승인 시."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "Path to the DE pipeline config YAML"
                        },
                        "corrections": {
                            "type": "object",
                            "description": "Correction mapping from validate result: {wrong_name: correct_name}",
                            "additionalProperties": {"type": "string"}
                        }
                    },
                    "required": ["corrections"]
                }
            },
            {
                "name": "run_de_pipeline",
                "description": (
                    "DE-GO 분석 Snakemake 파이프라인 실행. "
                    "⚠️ 반드시 dry_run=True(기본)로 먼저 실행해 작업 목록을 확인 후 사용자 승인을 받고 dry_run=False, background=True로 본 실행. "
                    "트리거: 'DE 분석 실행해줘', 'DE 파이프라인 돌려줘', 'GO 분석 시작해', 'snakemake DE 실행'."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "Path to DE config YAML. Auto-resolved from session if omitted."
                        },
                        "de_pipeline_dir": {
                            "type": "string",
                            "description": "Root dir of DE pipeline (contains Snakefile). Auto-resolved if omitted."
                        },
                        "cores": {"type": "integer", "default": 8},
                        "dry_run": {"type": "boolean", "default": True},
                        "background": {"type": "boolean", "default": False}
                    },
                    "required": []
                }
            },
            {
                "name": "monitor_de_pipeline",
                "description": (
                    "실행 중인 DE-GO 파이프라인 진행상황 확인. 로그 파싱 + PID 체크. "
                    "트리거: 'DE 진행상황 알려줘', 'DE 파이프라인 모니터링', 'DE 얼마나 됐어?', 'DE 완료됐어?'."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "de_config_file": {
                            "type": "string",
                            "description": "Path to DE config YAML. Auto-resolved from session if omitted."
                        },
                        "de_pipeline_dir": {
                            "type": "string",
                            "description": "Root dir of DE pipeline. Auto-resolved if omitted."
                        }
                    },
                    "required": []
                }
            },
            {
                "name": "get_pipeline_status",
                "description": "파이프라인 단계별 완료 여부를 파일시스템으로 확인 (실행 중이 아닐 때도 사용 가능). '어디까지 됐어?', '다음에 뭐 해야 해?', '어떤 단계가 완료됐어?', '파이프라인 단계 확인' → 이 도구 사용. 실행 중인 파이프라인 실시간 진행률은 monitor_pipeline 사용.",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
        ]

    def _load_session(self):
        """Restore session state from disk (survives agent restarts)."""
        if self.session_file.exists():
            try:
                with open(self.session_file) as f:
                    state = json.load(f)
                self.current_config = state.get("current_config") or self.current_config
                self.current_de_pipeline_dir = state.get("current_de_pipeline_dir") or self.current_de_pipeline_dir
                self.current_de_config = state.get("current_de_config") or self.current_de_config
                self.current_data_dir = state.get("current_data_dir") or self.current_data_dir
                self.current_project_id = state.get("current_project_id") or self.current_project_id
                self.current_results_dir = state.get("current_results_dir") or self.current_results_dir
                if self.current_config:
                    print(f"📂 이전 세션 복원: {self.current_config}")
            except Exception:
                pass  # Corrupt session file — start fresh

    def _save_session(self):
        """Persist session state to disk."""
        try:
            state = {
                "current_config": self.current_config,
                "current_de_pipeline_dir": self.current_de_pipeline_dir,
                "current_de_config": self.current_de_config,
                "current_data_dir": self.current_data_dir,
                "current_project_id": self.current_project_id,
                "current_results_dir": self.current_results_dir,
            }
            with open(self.session_file, "w") as f:
                json.dump(state, f, indent=2)
        except Exception:
            pass  # Non-fatal — session just won't persist

    def _resolve_config_file(self, arguments: Dict) -> Dict:
        """
        Auto-correct config_file in arguments:
        1. If LLM passed a placeholder (path/to/..., <...>, None) → use self.current_config
        2. If a real path is provided → update self.current_config for future calls
        """
        cfg = arguments.get('config_file')
        
        # Detect placeholder / hallucinated paths
        is_placeholder = (
            not cfg
            or cfg in ('path/to/config.yaml', '<config_file>', 'config.yaml', 'None', 'null')
            or cfg.startswith('path/to/')
            or cfg.startswith('<')
            or cfg.startswith('/path/')
        )
        
        if is_placeholder:
            if self.current_config:
                arguments = dict(arguments)
                arguments['config_file'] = self.current_config
            # else: leave as-is, tool will return a proper error
        else:
            # Real path provided — remember it for next time
            self.current_config = cfg
            self._save_session()

        return arguments

    def _resolve_de_pipeline_dir(self, arguments: Dict) -> Dict:
        """
        Auto-correct de_pipeline_dir in run_bridge arguments.
        1. If LLM passed a placeholder → use self.current_de_pipeline_dir
        2. If a real path is provided → update self.current_de_pipeline_dir
        """
        de = arguments.get('de_pipeline_dir')
        is_placeholder = (
            not de
            or de in ('None', 'null', '')
            or de.startswith('path/to/')
            or de.startswith('/path/')
            or de.startswith('<')
            or 'de/go/analysis/pipeline' in de
            or 'de-pipeline' in de
            or 'your_de_pipeline' in de
        )
        if is_placeholder:
            if self.current_de_pipeline_dir:
                arguments = dict(arguments)
                arguments['de_pipeline_dir'] = self.current_de_pipeline_dir
            # else: leave as-is — run_bridge will return a clear error
        else:
            self.current_de_pipeline_dir = de
            self._save_session()
        return arguments

    def _trim_tool_result_for_llm(self, result) -> Dict:
        """
        Return a compact version of a tool result for the LLM message.
        Removes large per-sample arrays (per_sample, multiqc_stats,
        fastqc_evaluation) that can overflow the context window, keeping
        only summary-level statistics that the LLM needs to answer the user.
        """
        if not isinstance(result, dict):
            return {"result": result}
        trimmed = {}
        for k, v in result.items():
            if k == "star_alignment" and isinstance(v, dict):
                # Keep summary stats; drop per_sample list
                trimmed[k] = {ek: ev for ek, ev in v.items() if ek != "per_sample"}
                low = v.get("low_mapping_samples", [])
                trimmed[k]["low_mapping_samples"] = low[:5]  # show at most 5
            elif k in ("multiqc_stats", "fastqc_evaluation", "samples"):
                # Replace with a count so the LLM knows data is present
                if isinstance(v, (list, dict)):
                    trimmed[k + "_count"] = len(v)
                else:
                    trimmed[k] = v
            else:
                trimmed[k] = v
        return trimmed

    def _execute_tool(self, tool_name: str, arguments: Dict) -> Any:
        """Execute a tool and return results."""
        
        # Auto-resolve config_file for all tools that use it
        if 'config_file' in arguments or tool_name in (
            'validate_input_data', 'run_pipeline', 'estimate_resources'
        ):
            arguments = self._resolve_config_file(arguments)
        
        if tool_name == "switch_project":
            # Clear all session state to start a new project
            self.current_config = None
            self.current_project_id = None
            self.current_data_dir = None
            self.current_results_dir = None
            self.current_de_pipeline_dir = None
            self.current_de_config = None
            new_cfg = arguments.get("new_config")
            if new_cfg:
                self.current_config = new_cfg
            self._save_session()
            msg = "세션 초기화 완료. 새 프로젝트를 시작할 준비가 됐습니다."
            if new_cfg:
                msg += f" Config: {new_cfg}"
            else:
                msg += " 새 프로젝트 config 경로를 알려주세요 (또는 '프로젝트 설정해줘'로 자동 설정)."
            return {"status": "success", "message": msg}

        elif tool_name == "get_project_status":
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "status"
            ], capture_output=True, text=True)
            return json.loads(result.stdout)
        
        elif tool_name == "compare_conditions":
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "compare"
            ], capture_output=True, text=True)
            return json.loads(result.stdout)
        
        elif tool_name == "get_failed_samples":
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "failed",
                "--format", "json"
            ], capture_output=True, text=True)
            return json.loads(result.stdout) if result.stdout else []
        
        elif tool_name == "get_sample_details":
            sample_id = arguments.get("sample_id")
            
            # Security: Validate sample_id
            if SECURITY_ENABLED:
                try:
                    sample_id = validate_sample_id(sample_id)
                except SecurityError as e:
                    return {"status": "error", "message": f"Security validation failed: {e}"}
            
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "sample",
                "--sample-id", sample_id
            ], capture_output=True, text=True)
            return json.loads(result.stdout)
        
        elif tool_name == "list_conditions":
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "conditions"
            ], capture_output=True, text=True)
            return json.loads(result.stdout)
        
        elif tool_name == "start_de_analysis":
            if not arguments.get("confirm"):
                return {"status": "cancelled", "message": "User confirmation required"}
            
            if not self.de_pipeline_dir:
                return {"status": "error", "message": "DE pipeline directory not configured"}
            
            # Execute bridge script
            result = subprocess.run([
                "python", "scripts/bridge_to_de_pipeline.py",
                "--rnaseq-output", str(self.rnaseq_output_dir),
                "--de-pipeline", str(self.de_pipeline_dir),
                "--project-id", self.project_id
            ], capture_output=True, text=True)
            
            return {
                "status": "started" if result.returncode == 0 else "failed",
                "output": result.stdout,
                "error": result.stderr
            }
        
        elif tool_name == "prepare_de_analysis":
            project_id = arguments.get("project_id", self.project_id)

            # Security: Validate project_id
            if SECURITY_ENABLED:
                try:
                    project_id = validate_project_id(project_id)
                except SecurityError as e:
                    return {"status": "error", "message": f"Security validation failed: {e}"}

            # config_file takes priority: let the bridge script parse paths from it directly
            config_file = arguments.get("config_file")
            rnaseq_output = arguments.get("rnaseq_output") or self.current_results_dir
            de_pipeline_dir = arguments.get("de_pipeline_dir") or self.current_de_pipeline_dir or self.de_pipeline_dir

            # Build command — omit --project-id when config_file is given (bridge parses it)
            cmd = [
                "conda", "run", "-n", "rna-seq-pipeline",
                "python", "scripts/bridge_to_de_pipeline.py",
                "--skip-de",
                "--yes"
            ]
            if config_file:
                cmd += ["--config", str(config_file)]
            else:
                cmd += ["--project-id", project_id]
                if rnaseq_output:
                    cmd += ["--rnaseq-output", str(rnaseq_output)]
            if de_pipeline_dir:
                cmd += ["--de-pipeline", str(de_pipeline_dir)]

            # Run bridge script in preparation mode (--skip-de)
            result = subprocess.run(
                cmd, capture_output=True, text=True, cwd="/data_3tb/shared/rna-seq-pipeline"
            )
            
            return {
                "status": "success" if result.returncode == 0 else "failed",
                "output": result.stdout,
                "error": result.stderr,
                "message": "DE analysis preparation complete" if result.returncode == 0 else "Preparation failed"
            }
        
        elif tool_name == "check_bridge_config":
            project_id = arguments.get("project_id")
            force = arguments.get("force_regenerate", False)
            
            # Security: Validate project_id
            if SECURITY_ENABLED:
                try:
                    project_id = validate_project_id(project_id)
                except SecurityError as e:
                    return {"status": "error", "message": f"Security validation failed: {e}"}
            
            # Import auto_config
            sys.path.insert(0, str(Path(__file__).parent.parent))
            from scripts.utils.auto_config import ensure_bridge_config
            
            result = ensure_bridge_config(project_id, force=force)
            return result
        
        elif tool_name == "validate_paths":
            project_id = arguments.get("project_id")
            
            # Security: Validate project_id
            if SECURITY_ENABLED:
                try:
                    project_id = validate_project_id(project_id)
                except SecurityError as e:
                    return {"status": "error", "message": f"Security validation failed: {e}"}
            
            # Import auto_config
            sys.path.insert(0, str(Path(__file__).parent.parent))
            from scripts.utils.auto_config import ensure_bridge_config, validate_paths
            
            # Get or create config
            config_result = ensure_bridge_config(project_id)
            
            if config_result['status'] == 'error':
                return config_result
            
            # Validate paths
            if 'config' in config_result:
                validation = validate_paths(config_result['config'])
            else:
                # Load existing config
                import yaml
                with open(config_result['config_path']) as f:
                    config = yaml.safe_load(f)
                validation = validate_paths(config)
            
            return {
                "status": "validated",
                "validation": validation,
                "all_valid": all(validation.values()),
                "missing_paths": [k for k, v in validation.items() if not v]
            }
        
        # Phase 8A: Pipeline Execution Tools
        elif tool_name == "get_sample_axes":
            result = subprocess.run([
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "axes"
            ], capture_output=True, text=True, cwd=str(_project_root))
            return json.loads(result.stdout) if result.stdout else {"error": result.stderr}

        elif tool_name == "compare_by_axis":
            axis = arguments.get("axis", "genotype")
            filters = arguments.get("filters", {})
            
            cmd = [
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "compare_axis",
                "--axis", axis
            ]
            for k, v in (filters or {}).items():
                cmd += ["--filter", f"{k}={v}"]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_project_root))
            return json.loads(result.stdout) if result.stdout else {"error": result.stderr}

        elif tool_name == "filter_samples":
            filters = arguments.get("filters", {})
            
            cmd = [
                "python", "scripts/standardization/agent_query.py",
                "--project-summary", str(self.project_summary_path),
                "--query", "filter"
            ]
            for k, v in filters.items():
                cmd += ["--filter", f"{k}={v}"]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_project_root))
            return json.loads(result.stdout) if result.stdout else {"error": result.stderr}

        elif tool_name == "create_project_config":
            try:
                from scripts.utils.pipeline_tools import create_project_config
                result = create_project_config(
                    project_id=arguments['project_id'],
                    data_dir=arguments['data_dir'],
                    results_dir=arguments['results_dir'],
                    species=arguments.get('species', 'human'),
                    genome_dir=arguments.get('genome_dir'),
                    genome_fasta=arguments.get('genome_fasta'),
                    annotation_gtf=arguments.get('annotation_gtf'),
                    star_index=arguments.get('star_index'),
                    genome_build=arguments.get('genome_build'),
                    use_sample_sheet=arguments.get('use_sample_sheet', False),
                    sample_sheet=arguments.get('sample_sheet'),
                    threads=arguments.get('threads', 12),
                    memory_gb=arguments.get('memory_gb', 48),
                    strandedness=arguments.get('strandedness', 0),
                )
                # Auto-capture config path so subsequent tools don't need it repeated
                if result.get('status') == 'success' and result.get('config_path'):
                    self.current_config = result['config_path']
                    self.current_project_id = arguments['project_id']
                    self.current_data_dir = arguments['data_dir']
                    self.current_results_dir = arguments['results_dir']
                    self._save_session()
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}
        
        elif tool_name == "detect_fastq_files":
            try:
                from scripts.utils.pipeline_tools import detect_fastq_files
                result = detect_fastq_files(arguments['data_dir'])
                # Remember the data_dir for subsequent create_project_config calls
                if result.get('status') != 'error':
                    self.current_data_dir = arguments['data_dir']
                    self._save_session()
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}
        
        elif tool_name == "validate_input_data":
            try:
                from scripts.utils.pipeline_tools import validate_input_data
                result = validate_input_data(arguments['config_file'])
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}
        
        elif tool_name == "validate_library_type":
            try:
                from scripts.utils.pipeline_tools import validate_library_type
                cfg_resolved = self._resolve_config_file(arguments)
                return validate_library_type(
                    config_file=cfg_resolved.get('config_file') or arguments['config_file'],
                    n_reads=arguments.get('n_reads', 500_000),
                    cores=arguments.get('cores', 8),
                    sample_id=arguments.get('sample_id'),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "setup_and_validate":
            try:
                from scripts.utils.pipeline_tools import setup_and_validate
                cfg_resolved = self._resolve_config_file(arguments)
                return setup_and_validate(
                    config_file=cfg_resolved.get('config_file') or arguments['config_file'],
                    n_reads=arguments.get('n_reads', 500_000),
                    cores=arguments.get('cores', 8),
                    sample_id=arguments.get('sample_id'),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "run_pipeline":
            try:
                from scripts.utils.pipeline_tools import run_pipeline
                result = run_pipeline(
                    config_file=arguments['config_file'],
                    cores=arguments.get('cores', 8),
                    dry_run=arguments.get('dry_run', True),
                    background=arguments.get('background', True),  # non-blocking by default for agent
                )
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}

        # ── Phase 8B tools ────────────────────────────────────────────────
        elif tool_name == "monitor_pipeline":
            try:
                from scripts.utils.pipeline_tools import monitor_pipeline
                import yaml as _yaml
                config_file = self._resolve_config_file(arguments).get('config_file') or arguments.get('config_file', '')
                # Auto-extract project_id and base_results_dir from config
                project_id = arguments.get('project_id', '')
                base_results_dir = arguments.get('base_results_dir', '')
                if (not project_id or not base_results_dir) and config_file:
                    with open(config_file) as _f:
                        _cfg = _yaml.safe_load(_f)
                    project_id = project_id or _cfg.get('project_id', '')
                    results_dir = _cfg.get('results_dir') or _cfg.get('base_results_dir', '')
                    if not base_results_dir and results_dir:
                        _rp = Path(results_dir)
                        base_results_dir = str(_rp.parent) if _rp.name == project_id else str(_rp)
                return monitor_pipeline(
                    project_id=project_id,
                    base_results_dir=base_results_dir,
                    config_file=config_file,
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "create_sample_sheet":
            try:
                from scripts.utils.pipeline_tools import create_sample_sheet
                return create_sample_sheet(
                    project_id=arguments['project_id'],
                    data_dir=arguments['data_dir'],
                    conditions=arguments.get('conditions'),
                    output_path=arguments.get('output_path'),
                    overwrite=arguments.get('overwrite', False),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "estimate_resources":
            try:
                from scripts.utils.pipeline_tools import estimate_resources
                return estimate_resources(
                    config_file=arguments['config_file'],
                    cores=arguments.get('cores', 8),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        # ── DE analysis bridge ────────────────────────────────────────────
        elif tool_name == "run_bridge":
            try:
                arguments = self._resolve_de_pipeline_dir(arguments)
                from scripts.utils.pipeline_tools import run_bridge
                result = run_bridge(
                    config_file=arguments['config_file'],
                    de_pipeline_dir=arguments['de_pipeline_dir'],
                    skip_de=arguments.get('skip_de', True),
                    dry_run=arguments.get('dry_run', False),
                    exclude_samples=arguments.get('exclude_samples', []),
                )
                # Auto-save DE config path for cognition tools
                if result.get("status") == "success" and result.get("de_config_file"):
                    self.current_de_config = result["de_config_file"]
                    self._save_session()
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}

        # ── Phase 8C: reading tools ───────────────────────────────────────
        elif tool_name == "read_project_config":
            try:
                from scripts.utils.pipeline_tools import read_project_config
                return read_project_config(arguments['config_file'])
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "read_sample_sheet":
            try:
                from scripts.utils.pipeline_tools import read_sample_sheet
                return read_sample_sheet(arguments['config_file'])
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "read_qc_results":
            try:
                from scripts.utils.pipeline_tools import read_qc_results
                return read_qc_results(arguments['config_file'])
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "read_counts":
            try:
                from scripts.utils.pipeline_tools import read_counts
                return read_counts(
                    config_file=arguments['config_file'],
                    genes=arguments.get('genes'),
                    top_n=arguments.get('top_n', 20),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "read_pipeline_logs":
            try:
                from scripts.utils.pipeline_tools import read_pipeline_logs
                return read_pipeline_logs(
                    config_file=arguments['config_file'],
                    rule=arguments.get('rule'),
                    sample_id=arguments.get('sample_id'),
                    tail_lines=arguments.get('tail_lines', 50),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "read_de_results_summary":
            try:
                from scripts.utils.pipeline_tools import read_de_results_summary
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                if not de_cfg:
                    return {"status": "error", "message": "DE config 경로가 없습니다. run_bridge를 먼저 실행하거나 de_config_file 경로를 직접 전달해 주세요."}
                return read_de_results_summary(
                    de_config_file=de_cfg,
                    pair=arguments.get("pair"),
                    top_n=arguments.get("top_n", 10),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "check_analysis_readiness":
            try:
                from scripts.utils.pipeline_tools import check_analysis_readiness
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                if not de_cfg:
                    return {"status": "error", "message": "DE config 경로가 없습니다. run_bridge를 먼저 실행하거나 de_config_file 경로를 직접 전달해 주세요."}
                return check_analysis_readiness(de_config_file=de_cfg)
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "validate_de_config_conditions":
            try:
                from scripts.utils.pipeline_tools import validate_de_config_conditions
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                if not de_cfg:
                    return {"status": "error", "message": "DE config 경로가 없습니다."}
                return validate_de_config_conditions(de_config_file=de_cfg)
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "apply_de_config_corrections":
            try:
                from scripts.utils.pipeline_tools import apply_de_config_corrections
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                if not de_cfg:
                    return {"status": "error", "message": "DE config 경로가 없습니다."}
                return apply_de_config_corrections(
                    de_config_file=de_cfg,
                    corrections=arguments["corrections"],
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "run_de_pipeline":
            try:
                from scripts.utils.pipeline_tools import run_de_pipeline
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                de_dir = arguments.get("de_pipeline_dir") or self.current_de_pipeline_dir
                if not de_cfg:
                    return {"status": "error", "message": "DE config 경로가 없습니다. de_config_file을 지정하거나 run_bridge를 먼저 실행하세요."}
                return run_de_pipeline(
                    de_config_file=de_cfg,
                    de_pipeline_dir=de_dir,
                    cores=arguments.get("cores", 8),
                    dry_run=arguments.get("dry_run", True),
                    background=arguments.get("background", False),
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "monitor_de_pipeline":
            try:
                from scripts.utils.pipeline_tools import monitor_de_pipeline
                de_cfg = arguments.get("de_config_file") or self.current_de_config
                de_dir = arguments.get("de_pipeline_dir") or self.current_de_pipeline_dir
                return monitor_de_pipeline(de_config_file=de_cfg, de_pipeline_dir=de_dir)
            except Exception as e:
                return {"status": "error", "message": str(e)}

        elif tool_name == "get_pipeline_status":
            try:
                from scripts.utils.pipeline_tools import get_pipeline_status
                return get_pipeline_status(
                    config_file=self.current_config,
                    data_dir=self.current_data_dir,
                    results_dir=self.current_results_dir,
                    de_pipeline_dir=self.current_de_pipeline_dir,
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        # ── Aliases for commonly hallucinated tool names ─────────────────
        elif tool_name in ("get_project_status", "check_pipeline_status", "pipeline_status"):
            # Redirect to monitor_pipeline using config_file
            try:
                from scripts.utils.pipeline_tools import monitor_pipeline
                import yaml as _yaml
                config_file = arguments.get('config_file', '')
                with open(config_file) as _f:
                    _cfg = _yaml.safe_load(_f)
                project_id = _cfg.get('project_id', '')
                results_dir = _cfg.get('results_dir') or f"{_cfg.get('base_results_dir','')}/{project_id}"
                base = str(Path(results_dir).parent) if project_id and results_dir.endswith(project_id) else results_dir
                return monitor_pipeline(
                    project_id=project_id,
                    base_results_dir=base,
                    config_file=config_file,
                )
            except Exception as e:
                return {"status": "error", "message": str(e)}

        else:
            available = [t["function"]["name"] for t in self._define_tools()]
            return {
                "error": f"Unknown tool: '{tool_name}'. Available tools: {available}. "
                         f"Please call one of the available tools."
            }
    
    def _build_pipeline_status_block(self) -> str:
        """파이프라인 진행 상태를 파일시스템에서 추론해 compact 텍스트로 반환."""
        try:
            from scripts.utils.pipeline_tools import get_pipeline_status
            status = get_pipeline_status(
                config_file=self.current_config,
                data_dir=self.current_data_dir,
                results_dir=self.current_results_dir,
                de_pipeline_dir=self.current_de_pipeline_dir,
            )
            completed = status.get("completed", [])
            steps = status.get("steps", {})
            done_labels = [steps[s]["label"] for s in completed if s in steps]
            pending_labels = [v["label"] for s, v in steps.items() if not v["done"]]
            lines = [
                "Pipeline progress:",
                f"  {status['summary']}",
            ]
            if done_labels:
                lines.append(f"  완료: {', '.join(done_labels)}")
            if pending_labels:
                lines.append(f"  미완료: {', '.join(pending_labels)}")
            lines.append(f"  → 다음 단계: {status['next_suggestion']}")
            return "\n".join(lines)
        except Exception:
            return ""

    def _build_system_prompt(self, native_tools: bool = False) -> str:
        """Build system prompt with project context and tool usage instructions.
        
        Args:
            native_tools: If True, omit TOOL_CALL text examples (tools are handled
                          by the API layer, not by the LLM's text output).
        """
        tools_description = "\n".join([
            f"- {tool['name']}: {tool['description']}"
            for tool in self.tools
        ])

        base = f"""You are an AI assistant for RNA-seq analysis pipeline management.

Project Context:
- Project ID: {self.project_id}
- Total samples: {self.project_summary.get('qc_summary', {}).get('total_samples', 'unknown')}
- QC pass rate: {self.project_summary.get('qc_summary', {}).get('pass_rate', 'N/A')}%
- Conditions: {', '.join(self.project_summary.get('condition_groups', {}).keys()) or 'not loaded'}

Available Tools:
{tools_description}

FILTERING RULES:
- If user mentions a specific tissue (HPC, PFC, Hippocampus, Prefrontal Cortex) → use compare_by_axis with filters, NOT compare_conditions
- If user mentions a specific sex (Male, Female) → use compare_by_axis with filters
- Available tissue values: "Hippocampus" (for HPC), "Prefrontal Cortex" (for PFC)
- Available condition values: "wildtype", "heterozygous"
- Available sex values: "Male", "Female"

## Reasoning Protocol — Think Before Acting

Before calling any tool, briefly assess the situation:
1. **Intent**: What does the user actually want? (not just the literal words)
2. **State check**: Does the current session have what's needed?
   - If config_file is 'not yet set' and a pipeline tool is requested → ask for config path first
   - If results_dir is missing and QC is requested → explain pipeline hasn't run yet
3. **Multi-step planning**: If the request requires >1 tool, identify the sequence before starting.
   Example: "DE 분석 준비해줘" → run_bridge → THEN validate_de_config_conditions (always)
4. **Anomaly detection**: After tool results arrive, check if the output looks suspicious
   (e.g., 0 samples, condition names that don't match, unusually low mapping rates).
   If something looks wrong, say so before moving on.

Keep this reasoning brief (1-2 sentences) — do NOT show the reasoning to the user unless it changes your action.

## Post-Tool Validation Rules

- After **run_bridge** completes → ALWAYS call `validate_de_config_conditions` automatically.
  Do not wait for user to ask. If mismatch found → show suggestions and ask "수정할까요?"
- After **validate_de_config_conditions** shows mismatches + user confirms → call `apply_de_config_corrections`.
- After **create_sample_sheet** with new conditions → remind user to re-run run_bridge if DE config already exists.
- Before recommending DE pipeline execution → suggest calling `check_analysis_readiness` first if not yet done.
- After DE pipeline completes → proactively offer `read_de_results_summary` ("결과 요약해드릴까요?").

CRITICAL TOOL ROUTING RULES (DO NOT VIOLATE):
- "DE 결과", "차발현", "유전자 변화", "D1 결과", "발현량 변화", "DE 요약" → read_de_results_summary (de_config_file 생략 가능)
- "QC 결과", "맵핑률", "샘플 통과" → get_project_status 또는 compare_conditions
- "파이프라인 진행" → monitor_pipeline
- compare_conditions는 QC 지표 비교 전용. DE 결과에 절대 사용 금지.

IMPORTANT:
1. Call the appropriate tool when action is needed
2. Extract parameters from user's natural language
3. If user's intent is unclear, ask for clarification
4. After tool execution, explain the results conversationally in Korean
5. For ALL pipeline tools (validate_input_data, run_pipeline, estimate_resources, monitor_pipeline):
   ALWAYS use the CURRENT CONFIG FILE: {self.current_config or 'not set — ask user for config file path'}
   Do NOT invent config file paths like "path/to/config.yaml" or "config.yaml"
6. For run_bridge: ALWAYS use the CURRENT DE PIPELINE DIR: {self.current_de_pipeline_dir or 'not set — ask user for DE pipeline directory path'}
   Do NOT invent de_pipeline_dir paths.

Current session state:
- config_file:     {self.current_config or 'not yet set'}
- project_id:      {self.current_project_id or 'not yet set'}
- data_dir (FASTQ): {self.current_data_dir or 'not yet set'}
- results_dir:     {self.current_results_dir or 'not yet set'}
- de_pipeline_dir: {self.current_de_pipeline_dir or 'not yet set'}
- de_config_file:  {self.current_de_config or 'not yet set'}  ← use for read_de_results_summary / check_analysis_readiness

{self._build_pipeline_status_block()}
Be conversational, clear, and actionable. When showing numbers, include units and context.
"""

        if native_tools:
            # For native tool calling, the API handles tool dispatch.
            # Do NOT include TOOL_CALL text examples — they confuse the model
            # into printing raw JSON in its reply instead of calling the API.
            return base

        # Text-pattern fallback: include explicit TOOL_CALL examples
        examples = """
HOW TO USE TOOLS:
When a user asks for information or action, respond with a tool call in JSON format:
```
TOOL_CALL: {"name": "tool_name", "parameters": {"param1": "value1"}}
```

Examples (Analysis Management):
- User: "QC 상태 보여줘", "QC 결과 어때?", "몇 개 통과했어?", "샘플 QC 요약" → TOOL_CALL: {{"name": "get_project_status", "parameters": {{}}}}  ← QC 결과 조회 (파이프라인 완료 후)
- User: "Ctrl_1 샘플 정보 알려줘" → TOOL_CALL: {{"name": "get_sample_details", "parameters": {{"sample_id": "Ctrl_1"}}}}
- User: "DE 분석 준비해줘" → TOOL_CALL: {{"name": "prepare_de_analysis", "parameters": {{"project_id": "{self.project_id}"}}}}

Examples (Multi-axis Analysis):
- User: "어떤 그룹이 있어?" → TOOL_CALL: {{"name": "get_sample_axes", "parameters": {{}}}}
- User: "HPC에서만 wildtype vs heterozygous 비교해줘" → TOOL_CALL: {{"name": "compare_by_axis", "parameters": {{"axis": "condition", "filters": {{"tissue": "Hippocampus"}}}}}}
- User: "Male HPC wildtype 샘플 목록" → TOOL_CALL: {{"name": "filter_samples", "parameters": {{"filters": {{"tissue": "Hippocampus", "sex": "Male", "condition": "wildtype"}}}}}}

Examples (Sample Sheet Creation):
IMPORTANT WORKFLOW — When user asks to create a sample sheet:
1. If you already have the sample list from detect_fastq_files, analyze sample name patterns to infer conditions.
2. Build a conditions mapping: {{"condition_name": ["prefix_or_substring", ...]}}
3. Call create_sample_sheet WITH the conditions parameter populated.
4. If conditions are unclear from names, ask the user: "샘플 이름에서 조건을 자동 추론했습니다: [list]. 맞나요?"
   Common patterns: *WT*/*Wild*/*Ctrl* → wildtype, *KO*/*Het*/*Mut* → mutant, *Tg*/*treat* → treatment, n*Tg* → non_stimulated
- User: "샘플 시트 만들어줘" (after detect shows MonWT1,MonWT2,MonTg1,nMonTg1) →
  TOOL_CALL: {{"name": "create_sample_sheet", "parameters": {{"project_id": "{proj}", "data_dir": "{ddir}", "conditions": {{"wildtype": ["MonWT"], "transgenic": ["MonTg"], "non_stimulated": ["nMonTg"]}}}}}}
- User: "MonWT는 wildtype, MonTg는 transgenic, nMonTg는 non_stimulated야, 샘플 시트 만들어줘" →
  TOOL_CALL: {{"name": "create_sample_sheet", "parameters": {{"project_id": "{proj}", "data_dir": "{ddir}", "conditions": {{"wildtype": ["MonWT"], "transgenic": ["MonTg"], "non_stimulated": ["nMonTg"]}}}}}}
- User: "샘플 시트 만들어줘" (generic, no prior detect) →
  TOOL_CALL: {{"name": "create_sample_sheet", "parameters": {{"project_id": "{proj}", "data_dir": "{ddir}"}}}}

Examples (Project Switching — 현재 세션을 초기화하고 새 프로젝트 시작):
- User: "새 프로젝트 시작할게", "다른 프로젝트 분석하고 싶어", "프로젝트 바꿔줘", "새로 시작" →
  TOOL_CALL: {{"name": "switch_project", "parameters": {{}}}}
- User: "mouse-monSTIM 프로젝트로 바꿔줘" →
  TOOL_CALL: {{"name": "switch_project", "parameters": {{"new_config": "config/projects/mouse-monSTIM-2026.yaml"}}}}
- switch_project 후에는 자동으로 현재 config가 비워지므로, 사용자에게 새 프로젝트 config 경로 또는 'setup' 실행을 안내할 것.

Examples (One-shot Setup — FASTQ 감지 + 샘플시트 + 라이브러리 검증 한번에):
- User: "프로젝트 설정해줘", "한번에 설정해줘", "setup", "설정부터 해줘", "처음부터 설정", "자동 설정" →
  TOOL_CALL: {{"name": "setup_and_validate", "parameters": {{"config_file": "{cfg}"}}}}
- User: "nMonTg4 샘플로 설정해줘" →
  TOOL_CALL: {{"name": "setup_and_validate", "parameters": {{"config_file": "{cfg}", "sample_id": "nMonTg4_S25"}}}}

Examples (Library Validation — run BEFORE main pipeline):
- User: "라이브러리 검증해줘", "strandedness 확인해줘", "mRNA-seq 맞는지 확인해줘", "library type 확인" →
  TOOL_CALL: {{"name": "validate_library_type", "parameters": {{"config_file": "{cfg}"}}}}

Examples (Pipeline Execution):
- User: "/data/raw/ 폴더에서 FASTQ 파일 찾아줘" → TOOL_CALL: {{"name": "detect_fastq_files", "parameters": {{"data_dir": "/data/raw"}}}}
- User: "파이프라인 dry-run 해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": true, "cores": 8}}}}
- User: "16 cores로 실행해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 16, "background": true}}}}
- User: "파이프라인 시작해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 8, "background": true}}}}
- User: "8코어로 돌려줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 8, "background": true}}}}
- User: "리소스 확인해줘" → TOOL_CALL: {{"name": "estimate_resources", "parameters": {{"config_file": "{cfg}", "cores": 8}}}}

Examples (Project/Pipeline Information):
- User: "이 프로젝트 설정이 어떻게 돼 있어?" → TOOL_CALL: {{"name": "read_project_config", "parameters": {{"config_file": "{cfg}"}}}}
- User: "샘플 목록 보여줘", "어떤 조건이 있어?", "wildtype 샘플 몇 개야?" → TOOL_CALL: {{"name": "read_sample_sheet", "parameters": {{"config_file": "{cfg}"}}}}
- User: "QC 결과 어때?", "매핑률 낮은 샘플 있어?", "어떤 샘플이 실패했어?" → TOOL_CALL: {{"name": "read_qc_results", "parameters": {{"config_file": "{cfg}"}}}}
- User: "파이프라인 진행상황 알려줘", "몇 개 완료됐어?", "얼마나 됐어?", "지금 뭐 하고 있어?", "show progress", "진행률 보여줘" → TOOL_CALL: {{"name": "monitor_pipeline", "parameters": {{"config_file": "{cfg}"}}}}  ← 실행 중 실시간 진행률
- User: "어디까지 됐어?", "다음에 뭐 해야 해?", "어떤 단계 완료됐어?", "파이프라인 몇 단계까지 왔어?" → TOOL_CALL: {{"name": "get_pipeline_status", "parameters": {{}}}}  ← 단계 완료 여부 (비실시간)
- User: "특정유전자 발현량 보여줘" → TOOL_CALL: {{"name": "read_counts", "parameters": {{"config_file": "{cfg}", "genes": ["GENE_NAME"]}}}}
- User: "가장 많이 발현된 유전자 보여줘" → TOOL_CALL: {{"name": "read_counts", "parameters": {{"config_file": "{cfg}", "top_n": 20}}}}
- User: "로그 보여줘", "에러 있어?" → TOOL_CALL: {{"name": "read_pipeline_logs", "parameters": {{"config_file": "{cfg}"}}}}
- User: "SAMPLE_ID STAR 로그 보여줘" → TOOL_CALL: {{"name": "read_pipeline_logs", "parameters": {{"config_file": "{cfg}", "rule": "star", "sample_id": "SAMPLE_ID"}}}}
- User: "DE 분석 연결해줘", "DE 파이프라인으로 넘겨줘", "counts 넘겨줘" → TOOL_CALL: {{"name": "run_bridge", "parameters": {{"config_file": "{cfg}", "de_pipeline_dir": "{de}", "skip_de": true}}}}
- User: "MonTg1_S17, MonTg5_S20 제외하고 DE 준비해줘" → TOOL_CALL: {{"name": "run_bridge", "parameters": {{"config_file": "{cfg}", "de_pipeline_dir": "{de}", "skip_de": true, "exclude_samples": ["MonTg1_S17", "MonTg5_S20"]}}}}
- User: "exclude-samples A B, DE 준비, rnaseq-output X, de-pipeline Y" → TOOL_CALL: {{"name": "run_bridge", "parameters": {{"config_file": "{cfg}", "de_pipeline_dir": "Y", "skip_de": true, "exclude_samples": ["A", "B"]}}}}

CRITICAL: If the user mentions "exclude", "제외", "빼고" with sample IDs → ALWAYS include "exclude_samples" in run_bridge parameters. Never omit it.

Examples (DE Config Validation — run_bridge 후 자동 실행):
- [After run_bridge success] → TOOL_CALL: {{"name": "validate_de_config_conditions", "parameters": {{"de_config_file": "<de_config_path from run_bridge result>"}}}}
- User: "DE config 검증해줘", "조건 이름 맞는지 확인해줘" → TOOL_CALL: {{"name": "validate_de_config_conditions", "parameters": {{"de_config_file": "<de_config_path>"}}}}
- [After validate shows mismatches + user says "수정해줘" or "고쳐줘"] → TOOL_CALL: {{"name": "apply_de_config_corrections", "parameters": {{"de_config_file": "<path>", "corrections": {{<suggestions from validate result>}}}}}}
- User: "자동으로 수정해줘" → apply_de_config_corrections with suggestions dict from previous validate result

Examples (DE 실행 전 진단):
- User: "DE 실행해도 돼?", "분석 준비됐어?", "샘플 수 맞아?", "DE 전에 확인해줘" → TOOL_CALL: {{"name": "check_analysis_readiness", "parameters": {{"de_config_file": "<de_config_path>"}}}}
- [Before suggesting snakemake DE run] → proactively call check_analysis_readiness first

Examples (DE-GO 파이프라인 실행):
- User: "DE 분석 실행해줘", "DE 파이프라인 돌려줘", "GO 분석 시작해", "snakemake DE 실행"
  → Step 1: TOOL_CALL: {{"name": "run_de_pipeline", "parameters": {{"dry_run": true}}}}
  → Step 2: dry-run 결과를 사용자에게 보여주고 "실행할까요?" 확인
  → Step 3 (승인 시): TOOL_CALL: {{"name": "run_de_pipeline", "parameters": {{"dry_run": false, "background": true}}}}
- User: "DE dry-run만 해줘", "DE 실행 전 확인해줘" → TOOL_CALL: {{"name": "run_de_pipeline", "parameters": {{"dry_run": true}}}}
- User: "DE 진행상황 알려줘", "DE 파이프라인 모니터링해줘", "DE 얼마나 됐어?", "DE 완료됐어?" → TOOL_CALL: {{"name": "monitor_de_pipeline", "parameters": {{}}}}

Examples (DE 결과 해석):
- User: "DE 결과 어때?", "결과 요약해줘", "어떤 유전자 변했어?", "분석 결과 보여줘" → TOOL_CALL: {{"name": "read_de_results_summary", "parameters": {{"de_config_file": "<de_config_path>"}}}}
- User: "D1 결과만 보여줘", "D3 vs Control 결과 어때?" → TOOL_CALL: {{"name": "read_de_results_summary", "parameters": {{"de_config_file": "<de_config_path>", "pair": "D1_vs_Control"}}}}
- [After DE pipeline completes successfully] → offer: "결과 요약해드릴까요?" and call read_de_results_summary if user agrees

CRITICAL RULE: "실행해줘", "돌려줘", "시작해줘", "run", "execute", "start pipeline" → ALWAYS use run_pipeline with dry_run=false.
"dry-run", "미리보기" → run_pipeline with dry_run=true.
"리소스", "resource", "얼마나 걸려" → estimate_resources (NOT run_pipeline).
"샘플", "조건", "condition", "sample list" → read_sample_sheet.
"QC", "mapping", "매핑", "품질" → read_qc_results.
"발현", "expression", "counts", "gene" → read_counts.
"로그", "log", "에러", "error" → read_pipeline_logs.
"설정", "config", "경로", "path", "어떻게 설정" → read_project_config.
"DE 분석", "DE 연결", "downstream", "counts 넘겨", "de pipeline" → run_bridge (ask for de_pipeline_dir if not given).
"""
        # Substitute session-state placeholders in examples
        _cfg = self.current_config or "config/projects/your_project.yaml"
        _de = self.current_de_pipeline_dir or "/path/to/de_pipeline"
        _proj = self.current_project_id or "your_project_id"
        _ddir = self.current_data_dir or "/path/to/fastq"
        examples = examples.replace("{cfg}", _cfg).replace("{de}", _de)
        examples = examples.replace("{proj}", _proj).replace("{ddir}", _ddir)
        return base + examples

    def _extract_tool_call(self, text: str) -> Optional[Dict]:
        """
        Robustly extract TOOL_CALL JSON from LLM response.
        
        Handles:
        - Single-line:  TOOL_CALL: {"name": ..., "parameters": {...}}
        - Multi-line:   TOOL_CALL:\n{\n  "name": ...\n}
        - Code blocks:  TOOL_CALL:\n```\n{...}\n```
        - Single quotes from LLM (replaced before parsing)
        """
        import re

        # Find the TOOL_CALL: marker and grab everything after it
        idx = text.find("TOOL_CALL:")
        if idx == -1:
            return None

        after = text[idx + len("TOOL_CALL:"):].strip()

        # Strip optional markdown code fences
        after = re.sub(r'^```[a-z]*\n?', '', after).strip()
        after = re.sub(r'```$', '', after).strip()

        # Try to find the outermost {...} block
        brace_start = after.find('{')
        if brace_start == -1:
            return None

        depth = 0
        end = -1
        for i, ch in enumerate(after[brace_start:], start=brace_start):
            if ch == '{':
                depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break

        if end == -1:
            return None

        json_str = after[brace_start:end]

        # Attempt 1: strict JSON
        try:
            return json.loads(json_str)
        except json.JSONDecodeError:
            pass

        # Attempt 2: replace single quotes → double quotes (common LLM mistake)
        try:
            fixed = json_str.replace("'", '"')
            return json.loads(fixed)
        except json.JSONDecodeError:
            pass

        # Attempt 3: fix unquoted true/false/null (Python-style booleans)
        try:
            fixed = re.sub(r'\bTrue\b', 'true', json_str)
            fixed = re.sub(r'\bFalse\b', 'false', fixed)
            fixed = re.sub(r'\bNone\b', 'null', fixed)
            return json.loads(fixed)
        except json.JSONDecodeError:
            pass

        return None

    def chat(self, user_message: str) -> str:
        """
        Process user message and return response.
        
        This is where the LLM integration happens.
        """
        
        if self.llm_provider == "ollama":
            return self._chat_ollama(user_message)
        elif self.llm_provider == "llamacpp":
            return self._chat_llamacpp(user_message)
        elif self.llm_provider == "openai":
            return self._chat_openai(user_message)
        elif self.llm_provider == "anthropic":
            return self._chat_anthropic(user_message)
    
    def _chat_openai(self, user_message: str) -> str:
        """OpenAI GPT-4 implementation."""
        
        messages = [
            {"role": "system", "content": self._build_system_prompt()},
            {"role": "user", "content": user_message}
        ]
        
        # Initial API call with function calling
        response = self.llm_client.ChatCompletion.create(
            model=self.model,
            messages=messages,
            functions=self.tools,
            function_call="auto"
        )
        
        message = response.choices[0].message
        
        # Check if function call is needed
        if message.get("function_call"):
            function_name = message["function_call"]["name"]
            arguments = json.loads(message["function_call"]["arguments"])
            
            # Execute the function
            function_result = self._execute_tool(function_name, arguments)
            
            # Send function result back to GPT
            messages.append(message)
            messages.append({
                "role": "function",
                "name": function_name,
                "content": json.dumps(function_result)
            })
            
            # Get final response
            second_response = self.llm_client.ChatCompletion.create(
                model=self.model,
                messages=messages
            )
            
            return second_response.choices[0].message["content"]
        
        else:
            return message["content"]
    
    def _chat_anthropic(self, user_message: str) -> str:
        """Anthropic Claude implementation."""
        # Similar structure to OpenAI but using Anthropic's API
        # Implementation details depend on Anthropic's function calling API
        pass
    
    # Models known to support Ollama native tool calling
    _NATIVE_TOOL_CALLING_MODELS = {
        "qwen3", "qwen2.5", "qwen2", "llama3.1", "llama3.2", "llama3.3",
        "mistral", "mistral-nemo", "mixtral", "firefunction",
        "command-r", "smollm2", "hermes3", "granite3",
    }

    def _supports_native_tool_calling(self) -> bool:
        """Check if the current model supports Ollama native tool calling."""
        model_base = self.model.split(":")[0].lower()
        return any(model_base.startswith(m) for m in self._NATIVE_TOOL_CALLING_MODELS)

    def _chat_ollama(self, user_message: str) -> str:
        """
        Ollama local LLM implementation.

        Strategy:
        - Models that support native tool calling (qwen2.5, llama3.1, etc.):
          Use Ollama tools= API for reliable structured function calls.
        - Older / smaller models: fall back to TOOL_CALL text-pattern matching.
        """
        if self._supports_native_tool_calling():
            return self._chat_ollama_native_tools(user_message)
        else:
            return self._chat_ollama_text_pattern(user_message)

    def _chat_ollama_native_tools(self, user_message: str) -> str:
        """
        Ollama native tool calling (for models like qwen2.5, llama3.1+).

        Uses the tools= parameter in ollama.chat(), which guarantees
        structured JSON output — no regex parsing needed.
        """
        messages = [
            {"role": "system", "content": self._build_system_prompt(native_tools=True)},
            {"role": "user",   "content": user_message},
        ]
        ollama_tools = self._convert_tools_to_ollama_format()

        try:
            response = ollama.chat(
                model=self.model,
                messages=messages,
                tools=ollama_tools,
                options={"temperature": 0.1, "num_ctx": 8192},
            )

            msg = response["message"]

            # ── Native tool call returned ──────────────────────────────────
            tool_calls = msg.get("tool_calls") or []
            if tool_calls:
                results = []
                for tc in tool_calls:
                    fn    = tc["function"]
                    name  = fn["name"]
                    args  = fn.get("arguments", {})
                    if isinstance(args, str):
                        args = json.loads(args)

                    print(f"\n🔧 Calling tool: {name}")
                    print(f"📋 Parameters: {args}")

                    tool_result = self._execute_tool(name, args)
                    results.append({"tool": name, "result": tool_result})

                    # If tool was unknown, include a correction hint
                    if "error" in tool_result and "Unknown tool" in str(tool_result.get("error", "")):
                        tool_content = json.dumps({
                            "error": tool_result["error"],
                            "hint": "Use one of the available tools listed above."
                        }, ensure_ascii=False)
                    else:
                        # Strip large arrays before sending to LLM to avoid context overflow.
                        # Keep only summary-level fields; per-sample lists can be thousands of tokens.
                        tool_content = json.dumps(
                            self._trim_tool_result_for_llm(tool_result),
                            ensure_ascii=False
                        )

                    # Append tool turn to conversation
                    messages.append({"role": "assistant", "content": None,
                                     "tool_calls": [tc]})
                    messages.append({
                        "role": "tool",
                        "content": tool_content,
                    })

                # Ask LLM to summarise all tool results.
                # Do NOT pass tools= here — we only want a plain text summary,
                # not another round of tool calls or code-block suggestions.
                messages.append({
                    "role": "user",
                    "content": (
                        "위 도구 실행 결과를 사용자에게 친절하게 한국어로 요약해줘. "
                        "코드 블록이나 추가 명령어 예시는 포함하지 말고, "
                        "결과 내용만 자연스럽게 설명해줘."
                    ),
                })
                summary = ollama.chat(
                    model=self.model,
                    messages=messages,
                    options={"temperature": 0.7, "num_ctx": 8192},
                )
                content = summary["message"].get("content", "").strip()
                if not content:
                    # LLM returned empty — build a fallback summary from raw results
                    parts = []
                    for r in results:
                        parts.append(f"[{r['tool']}] {json.dumps(r['result'], ensure_ascii=False)[:500]}")
                    content = "도구 실행 결과:\n" + "\n".join(parts)
                return content

            # ── Plain text reply (no tool needed) ─────────────────────────
            return msg.get("content", "")

        except Exception as e:
            return f"❌ Error communicating with Ollama: {e}"

    def _chat_ollama_text_pattern(self, user_message: str) -> str:
        """
        Fallback: text-pattern TOOL_CALL matching for models without
        native tool calling support.
        """
        messages = [
            {"role": "system", "content": self._build_system_prompt(native_tools=False)},
            {"role": "user",   "content": user_message},
        ]

        try:
            response = ollama.chat(
                model=self.model,
                messages=messages,
                options={"temperature": 0.3, "num_ctx": 4096},
            )

            response_text = response["message"]["content"]

            if "TOOL_CALL:" in response_text:
                import re
                tool_call = self._extract_tool_call(response_text)
                if tool_call:
                    function_name = tool_call["name"]
                    arguments     = tool_call.get("parameters", {})

                    print(f"\n🔧 Calling tool: {function_name}")
                    print(f"📋 Parameters: {arguments}")

                    function_result = self._execute_tool(function_name, arguments)
                    result_str      = json.dumps(function_result, indent=2, ensure_ascii=False)

                    messages.append({"role": "assistant", "content": response_text})
                    messages.append({
                        "role": "user",
                        "content": (
                            f"Tool result:\n{result_str}\n\n"
                            "Please explain this result to the user in a friendly way."
                        ),
                    })

                    explanation = ollama.chat(
                        model=self.model,
                        messages=messages,
                        options={"temperature": 0.7},
                    )
                    return explanation["message"]["content"]
                else:
                    return response_text

            return response_text

        except Exception as e:
            return f"❌ Error communicating with Ollama: {e}"

    
    def _chat_llamacpp(self, user_message: str) -> str:
        """
        llama.cpp implementation with function calling simulation.
        
        llama.cpp doesn't natively support function calling, so we use
        prompt engineering to simulate it.
        """
        
        # Build prompt with tool descriptions
        tools_description = "\n".join([
            f"- {tool['name']}: {tool['description']}"
            for tool in self.tools
        ])
        
        prompt = f"""{self._build_system_prompt()}

Available commands:
{tools_description}

To use a command, respond with:
TOOL: <command_name>
ARGS: <json_arguments>

User: {user_message}
Assistant:"""
        
        try:
            response = self.llm_client(
                prompt,
                max_tokens=512,
                temperature=0.7,
                stop=["User:", "\n\n"]
            )
            
            response_text = response['choices'][0]['text'].strip()
            
            # Parse for tool calls
            if "TOOL:" in response_text and "ARGS:" in response_text:
                lines = response_text.split('\n')
                tool_name = None
                args_json = None
                
                for line in lines:
                    if line.startswith("TOOL:"):
                        tool_name = line.replace("TOOL:", "").strip()
                    elif line.startswith("ARGS:"):
                        args_json = line.replace("ARGS:", "").strip()
                
                if tool_name and args_json:
                    try:
                        arguments = json.loads(args_json)
                        function_result = self._execute_tool(tool_name, arguments)
                        
                        # Generate final response with results
                        final_prompt = f"{prompt}\n\nCommand result:\n{json.dumps(function_result, indent=2)}\n\nProvide a natural language summary:"
                        
                        final_response = self.llm_client(
                            final_prompt,
                            max_tokens=512,
                            temperature=0.7
                        )
                        
                        return final_response['choices'][0]['text'].strip()
                    
                    except json.JSONDecodeError:
                        return f"Error parsing arguments: {args_json}"
            
            return response_text
            
        except Exception as e:
            return f"Error with llama.cpp: {str(e)}"
    
    def _convert_tools_to_ollama_format(self):
        """Convert OpenAI function format to Ollama tools format."""
        ollama_tools = []
        
        for tool in self.tools:
            ollama_tool = {
                "type": "function",
                "function": {
                    "name": tool["name"],
                    "description": tool["description"],
                    "parameters": tool["parameters"]
                }
            }
            ollama_tools.append(ollama_tool)
        
        return ollama_tools


def main():
    parser = argparse.ArgumentParser(
        description="Natural language agent for RNA-seq pipeline"
    )
    parser.add_argument(
        '--project-summary',
        type=Path,
        default=None,
        help='Path to project_summary.json (optional — not required for pipeline execution tools)'
    )
    parser.add_argument(
        '--rnaseq-output',
        type=Path,
        default=None,
        help='RNA-seq pipeline output directory (optional)'
    )
    parser.add_argument(
        '--de-pipeline',
        type=Path,
        help='DE/GO pipeline directory (optional)'
    )
    parser.add_argument(
        '--llm-provider',
        choices=['ollama', 'llamacpp', 'openai', 'anthropic'],
        default='ollama',
        help='LLM provider to use (default: ollama for data security)'
    )
    parser.add_argument(
        '--api-key',
        help='API key for cloud LLM providers (OpenAI/Anthropic)'
    )
    parser.add_argument(
        '--model',
        help='Model name (e.g., llama3.1:8b, gpt-4, claude-3-sonnet)'
    )
    parser.add_argument(
        '--ollama-host',
        default='http://localhost:11434',
        help='Ollama server URL (default: http://localhost:11434)'
    )
    parser.add_argument(
        '--config',
        help='Resume from existing project config YAML (sets current_config for this session)'
    )
    parser.add_argument(
        '--reset-session',
        action='store_true',
        help='Clear saved session state and start fresh'
    )
    parser.add_argument(
        '--interactive',
        action='store_true',
        help='Start interactive chat session'
    )
    parser.add_argument(
        '--message',
        help='Single message to process (non-interactive mode)'
    )

    args = parser.parse_args()

    # Handle session reset before agent init
    session_file = _project_root / ".agent_session.json"
    if args.reset_session and session_file.exists():
        session_file.unlink()
        print("🗑️  세션 초기화 완료")

    # Initialize agent
    try:
        agent = PipelineAgent(
            project_summary_path=args.project_summary,
            rnaseq_output_dir=args.rnaseq_output,
            de_pipeline_dir=args.de_pipeline,
            llm_provider=args.llm_provider,
            api_key=args.api_key,
            model=args.model,
            ollama_host=args.ollama_host
        )
    except Exception as e:
        print(f"❌ Failed to initialize agent: {e}")
        sys.exit(1)

    # --config overrides session file (allows explicit resume)
    if args.config:
        agent.current_config = args.config
        agent._save_session()
        print(f"📂 Config 설정: {args.config}")

    if args.interactive:
        # Interactive mode
        print(f"🤖 RNA-seq Pipeline Agent")
        if agent.project_id != "unknown":
            print(f"Project: {agent.project_id}")
        else:
            print(f"Mode: Pipeline execution (no project summary loaded)")
            print(f"      Use detect_fastq_files, create_project_config, run_pipeline, etc.")
        if agent.current_config:
            print(f"현재 프로젝트: {agent.current_config}")
        print(f"Type 'exit' or 'quit' to end session\n")
        
        while True:
            try:
                user_input = input("You: ")
                if user_input.lower() in ['exit', 'quit', 'bye']:
                    print("Goodbye!")
                    break
                
                response = agent.chat(user_input)
                print(f"\nAgent: {response}\n")
                
            except KeyboardInterrupt:
                print("\nGoodbye!")
                break
            except Exception as e:
                print(f"❌ Error: {e}")
    
    elif args.message:
        # Single message mode
        response = agent.chat(args.message)
        print(response)
    
    else:
        print("Error: Specify --interactive or --message")
        sys.exit(1)


if __name__ == '__main__':
    main()
