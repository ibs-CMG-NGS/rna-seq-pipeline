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
                 ollama_host: str = "http://localhost:11434"):
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
            self.model = model or "qwen2.5:32b"  # Default to Qwen2.5 32B
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
    
    def _define_tools(self) -> List[Dict]:
        """Define available tools for LLM function calling."""
        return [
            {
                "name": "get_project_status",
                "description": "Get overall project QC status and summary",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": []
                }
            },
            {
                "name": "compare_conditions",
                "description": "Compare QC metrics across experimental conditions",
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
                "description": "Execute RNA-seq pipeline with Snakemake. Use dry_run=false to ACTUALLY RUN the pipeline. Use dry_run=true only for preview. When user says '실행해줘', '돌려줘', '시작해줘', 'run', 'execute', 'start' → set dry_run=false.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "cores": {"type": "integer", "description": "Number of cores (e.g. 16)", "default": 8},
                        "dry_run": {"type": "boolean", "description": "false=actually run pipeline, true=preview only. Default false when user says run/execute/start.", "default": False}
                    },
                    "required": ["config_file"]
                }
            },
            # Phase 8B tools
            {
                "name": "monitor_pipeline",
                "description": "Check pipeline execution status and progress by inspecting output files",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {"type": "string", "description": "Project identifier"},
                        "base_results_dir": {"type": "string", "description": "Base results directory"},
                        "config_file": {"type": "string", "description": "Path to config.yaml (optional)"}
                    },
                    "required": ["project_id", "base_results_dir"]
                }
            },
            {
                "name": "create_sample_sheet",
                "description": "Generate a sample metadata TSV from detected FASTQ files, optionally assigning conditions",
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
                        "output_path": {"type": "string", "description": "Output TSV path (optional)"}
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
                    "Set skip_de=false only if user explicitly asks to RUN DE analysis immediately."
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
        ]
    
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
        return arguments

    def _trim_tool_result_for_llm(self, result: Dict) -> Dict:
        """
        Return a compact version of a tool result for the LLM message.
        Removes large per-sample arrays (per_sample, multiqc_stats,
        fastqc_evaluation) that can overflow the context window, keeping
        only summary-level statistics that the LLM needs to answer the user.
        """
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
        
        if tool_name == "get_project_status":
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
            
            # Run bridge script in preparation mode (--skip-de)
            result = subprocess.run([
                "conda", "run", "-n", "rna-seq-pipeline",
                "python", "scripts/bridge_to_de_pipeline.py",
                "--project-id", project_id,
                "--skip-de",
                "--yes"
            ], capture_output=True, text=True, cwd="/data_3tb/shared/rna-seq-pipeline")
            
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
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}
        
        elif tool_name == "detect_fastq_files":
            try:
                from scripts.utils.pipeline_tools import detect_fastq_files
                result = detect_fastq_files(arguments['data_dir'])
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
        
        elif tool_name == "run_pipeline":
            try:
                from scripts.utils.pipeline_tools import run_pipeline
                result = run_pipeline(
                    config_file=arguments['config_file'],
                    cores=arguments.get('cores', 8),
                    dry_run=arguments.get('dry_run', True)
                )
                return result
            except Exception as e:
                return {"status": "error", "message": str(e)}

        # ── Phase 8B tools ────────────────────────────────────────────────
        elif tool_name == "monitor_pipeline":
            try:
                from scripts.utils.pipeline_tools import monitor_pipeline
                return monitor_pipeline(
                    project_id=arguments['project_id'],
                    base_results_dir=arguments['base_results_dir'],
                    config_file=arguments.get('config_file'),
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
                return run_bridge(
                    config_file=arguments['config_file'],
                    de_pipeline_dir=arguments['de_pipeline_dir'],
                    skip_de=arguments.get('skip_de', True),
                    dry_run=arguments.get('dry_run', False),
                )
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

Current session config: {self.current_config or 'not yet set'}
Current DE pipeline dir: {self.current_de_pipeline_dir or 'not yet set'}

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
- User: "QC 상태 보여줘" → TOOL_CALL: {{"name": "get_project_status", "parameters": {{}}}}
- User: "Ctrl_1 샘플 정보 알려줘" → TOOL_CALL: {{"name": "get_sample_details", "parameters": {{"sample_id": "Ctrl_1"}}}}
- User: "DE 분석 준비해줘" → TOOL_CALL: {{"name": "prepare_de_analysis", "parameters": {{"project_id": "{self.project_id}"}}}}

Examples (Multi-axis Analysis):
- User: "어떤 그룹이 있어?" → TOOL_CALL: {{"name": "get_sample_axes", "parameters": {{}}}}
- User: "HPC에서만 wildtype vs heterozygous 비교해줘" → TOOL_CALL: {{"name": "compare_by_axis", "parameters": {{"axis": "condition", "filters": {{"tissue": "Hippocampus"}}}}}}
- User: "Male HPC wildtype 샘플 목록" → TOOL_CALL: {{"name": "filter_samples", "parameters": {{"filters": {{"tissue": "Hippocampus", "sex": "Male", "condition": "wildtype"}}}}}}

Examples (Pipeline Execution):
- User: "/data/raw/ 폴더에서 FASTQ 파일 찾아줘" → TOOL_CALL: {{"name": "detect_fastq_files", "parameters": {{"data_dir": "/data/raw"}}}}
- User: "파이프라인 dry-run 해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": true, "cores": 8}}}}
- User: "16 cores로 실행해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 16}}}}
- User: "파이프라인 시작해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 8}}}}
- User: "8코어로 돌려줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "{cfg}", "dry_run": false, "cores": 8}}}}
- User: "리소스 확인해줘" → TOOL_CALL: {{"name": "estimate_resources", "parameters": {{"config_file": "{cfg}", "cores": 8}}}}

Examples (Project/Pipeline Information):
- User: "이 프로젝트 설정이 어떻게 돼 있어?" → TOOL_CALL: {{"name": "read_project_config", "parameters": {{"config_file": "{cfg}"}}}}
- User: "샘플 목록 보여줘", "어떤 조건이 있어?", "wildtype 샘플 몇 개야?" → TOOL_CALL: {{"name": "read_sample_sheet", "parameters": {{"config_file": "{cfg}"}}}}
- User: "QC 결과 어때?", "매핑률 낮은 샘플 있어?", "어떤 샘플이 실패했어?" → TOOL_CALL: {{"name": "read_qc_results", "parameters": {{"config_file": "{cfg}"}}}}
- User: "현재 상태 보여줘", "진행률 어때?", "얼마나 됐어?", "show status" → TOOL_CALL: {{"name": "monitor_pipeline", "parameters": {{"config_file": "{cfg}"}}}}
- User: "특정유전자 발현량 보여줘" → TOOL_CALL: {{"name": "read_counts", "parameters": {{"config_file": "{cfg}", "genes": ["GENE_NAME"]}}}}
- User: "가장 많이 발현된 유전자 보여줘" → TOOL_CALL: {{"name": "read_counts", "parameters": {{"config_file": "{cfg}", "top_n": 20}}}}
- User: "로그 보여줘", "에러 있어?" → TOOL_CALL: {{"name": "read_pipeline_logs", "parameters": {{"config_file": "{cfg}"}}}}
- User: "SAMPLE_ID STAR 로그 보여줘" → TOOL_CALL: {{"name": "read_pipeline_logs", "parameters": {{"config_file": "{cfg}", "rule": "star", "sample_id": "SAMPLE_ID"}}}}
- User: "DE 분석 연결해줘", "DE 파이프라인으로 넘겨줘", "counts 넘겨줘" → TOOL_CALL: {{"name": "run_bridge", "parameters": {{"config_file": "{cfg}", "de_pipeline_dir": "{de}", "skip_de": true}}}}

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
        examples = examples.replace("{cfg}", _cfg).replace("{de}", _de)
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
        "qwen2.5", "qwen2", "llama3.1", "llama3.2", "llama3.3",
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
        '--interactive',
        action='store_true',
        help='Start interactive chat session'
    )
    parser.add_argument(
        '--message',
        help='Single message to process (non-interactive mode)'
    )
    
    args = parser.parse_args()
    
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
    
    if args.interactive:
        # Interactive mode
        print(f"🤖 RNA-seq Pipeline Agent")
        if agent.project_id != "unknown":
            print(f"Project: {agent.project_id}")
        else:
            print(f"Mode: Pipeline execution (no project summary loaded)")
            print(f"      Use detect_fastq_files, create_project_config, run_pipeline, etc.")
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
