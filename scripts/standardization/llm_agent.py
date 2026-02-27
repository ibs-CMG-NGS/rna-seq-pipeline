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
sys.path.insert(0, str(Path(__file__).parent.parent))
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
    print("⚠️  Warning: Security utilities not available")
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
                 project_summary_path: Path,
                 rnaseq_output_dir: Path,
                 de_pipeline_dir: Path = None,
                 llm_provider: str = "ollama",
                 api_key: str = None,
                 model: str = None,
                 ollama_host: str = "http://localhost:11434"):
        """
        Initialize agent.
        
        Args:
            project_summary_path: Path to project_summary.json
            rnaseq_output_dir: RNA-seq pipeline output directory
            de_pipeline_dir: DE/GO pipeline directory (optional)
            llm_provider: "ollama" (local, recommended), "openai", "anthropic", "llamacpp"
            api_key: API key for cloud LLM providers
            model: Model name (e.g., "llama3.1", "gpt-4", "claude-3-sonnet")
            ollama_host: Ollama server URL (default: http://localhost:11434)
        """
        self.project_summary_path = Path(project_summary_path)
        self.rnaseq_output_dir = Path(rnaseq_output_dir)
        self.de_pipeline_dir = Path(de_pipeline_dir) if de_pipeline_dir else None
        self.llm_provider = llm_provider
        self.ollama_host = ollama_host
        
        # Load project summary
        with open(self.project_summary_path) as f:
            self.project_summary = json.load(f)
        
        self.project_id = self.project_summary['project_id']
        
        # Initialize LLM client
        if llm_provider == "ollama" and HAS_OLLAMA:
            self.llm_client = ollama
            self.model = model or "llama3.1:8b"  # Default to Llama 3.1 8B
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
            # Phase 8A: Pipeline Execution Tools
            {
                "name": "create_project_config",
                "description": "Create new project configuration file from parameters",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "project_id": {"type": "string", "description": "Project identifier"},
                        "data_dir": {"type": "string", "description": "Directory containing FASTQ files"},
                        "results_dir": {"type": "string", "description": "Output directory for results"},
                        "species": {"type": "string", "description": "Species (human/mouse)", "default": "human"}
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
                "description": "Execute RNA-seq pipeline with Snakemake",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "config_file": {"type": "string", "description": "Path to config.yaml"},
                        "cores": {"type": "integer", "description": "Number of cores", "default": 8},
                        "dry_run": {"type": "boolean", "description": "Show what would run without executing", "default": True}
                    },
                    "required": ["config_file"]
                }
            }
        ]
    
    def _execute_tool(self, tool_name: str, arguments: Dict) -> Any:
        """Execute a tool and return results."""
        
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
        elif tool_name == "create_project_config":
            try:
                from scripts.utils.pipeline_tools import create_project_config
                result = create_project_config(
                    project_id=arguments['project_id'],
                    data_dir=arguments['data_dir'],
                    results_dir=arguments['results_dir'],
                    species=arguments.get('species', 'human')
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
        
        else:
            return {"error": f"Unknown tool: {tool_name}"}
    
    def _build_system_prompt(self) -> str:
        """Build system prompt with project context and tool usage instructions."""
        tools_description = "\n".join([
            f"- {tool['name']}: {tool['description']}"
            for tool in self.tools
        ])
        
        return f"""You are an AI assistant for RNA-seq analysis pipeline management.

Project Context:
- Project ID: {self.project_id}
- Total samples: {self.project_summary['qc_summary']['total_samples']}
- QC pass rate: {self.project_summary['qc_summary']['pass_rate']}%
- Conditions: {', '.join(self.project_summary['condition_groups'].keys())}

Available Tools:
{tools_description}

HOW TO USE TOOLS:
When a user asks for information or action, respond with a tool call in JSON format:
```
TOOL_CALL: {{"name": "tool_name", "parameters": {{"param1": "value1"}}}}
```

Examples (Analysis Management):
- User: "QC 상태 보여줘" → TOOL_CALL: {{"name": "get_project_status", "parameters": {{}}}}
- User: "Ctrl_1 샘플 정보 알려줘" → TOOL_CALL: {{"name": "get_sample_details", "parameters": {{"sample_id": "Ctrl_1"}}}}
- User: "경로 검증해줘" → TOOL_CALL: {{"name": "validate_paths", "parameters": {{"project_id": "{self.project_id}"}}}}
- User: "DE 분석 준비해줘" → TOOL_CALL: {{"name": "prepare_de_analysis", "parameters": {{"project_id": "{self.project_id}"}}}}

Examples (Pipeline Execution):
- User: "새 프로젝트 설정 만들어줘. test-2026, 데이터는 /data/raw/" → TOOL_CALL: {{"name": "create_project_config", "parameters": {{"project_id": "test-2026", "data_dir": "/data/raw", "results_dir": "/data/output/test-2026", "species": "human"}}}}
- User: "/data/raw/ 폴더에서 FASTQ 파일 찾아줘" → TOOL_CALL: {{"name": "detect_fastq_files", "parameters": {{"data_dir": "/data/raw"}}}}
- User: "입력 데이터 검증해줘" → TOOL_CALL: {{"name": "validate_input_data", "parameters": {{"config_file": "config/projects/test-2026.yaml"}}}}
- User: "파이프라인 dry-run 해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "config/projects/test-2026.yaml", "dry_run": true}}}}
- User: "파이프라인 실행해줘" → TOOL_CALL: {{"name": "run_pipeline", "parameters": {{"config_file": "config/projects/test-2026.yaml", "dry_run": false, "cores": 8}}}}

IMPORTANT:
1. Always use TOOL_CALL format when action is needed
2. Extract parameters from user's natural language
3. If user's intent is unclear, ask for clarification
4. After tool execution, explain the results conversationally

Be conversational, clear, and actionable. When showing numbers, include units and context.
"""
    
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
    
    def _chat_ollama(self, user_message: str) -> str:
        """
        Ollama local LLM implementation with simplified tool calling.
        
        Uses pattern matching for TOOL_CALL instead of native function calling.
        """
        
        messages = [
            {
                "role": "system",
                "content": self._build_system_prompt()
            },
            {
                "role": "user",
                "content": user_message
            }
        ]
        
        try:
            # Get LLM response
            response = ollama.chat(
                model=self.model,
                messages=messages,
                options={
                    "temperature": 0.3,  # Lower temperature for more consistent tool calling
                    "num_ctx": 4096
                }
            )
            
            response_text = response['message']['content']
            
            # Check if response contains TOOL_CALL
            if "TOOL_CALL:" in response_text:
                # Extract JSON from TOOL_CALL
                import re
                match = re.search(r'TOOL_CALL:\s*({.*?})', response_text, re.DOTALL)
                if match:
                    try:
                        tool_call = json.loads(match.group(1))
                        function_name = tool_call['name']
                        arguments = tool_call.get('parameters', {})
                        
                        # Execute the tool
                        print(f"\n🔧 Calling tool: {function_name}")
                        print(f"📋 Parameters: {arguments}")
                        
                        function_result = self._execute_tool(function_name, arguments)
                        
                        # Format result for user
                        result_str = json.dumps(function_result, indent=2, ensure_ascii=False)
                        
                        # Ask LLM to explain the result
                        messages.append({"role": "assistant", "content": response_text})
                        messages.append({
                            "role": "user",
                            "content": f"Tool result:\n{result_str}\n\nPlease explain this result to the user in a friendly way."
                        })
                        
                        explanation = ollama.chat(
                            model=self.model,
                            messages=messages,
                            options={"temperature": 0.7}
                        )
                        
                        return explanation['message']['content']
                    
                    except json.JSONDecodeError as e:
                        return f"❌ Error parsing tool call: {e}\n\nPlease try rephrasing your request."
                
            # No tool call needed, return response as-is
            return response_text
            
        except Exception as e:
            return f"❌ Error communicating with Ollama: {str(e)}"

    
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
        required=True,
        help='Path to project_summary.json'
    )
    parser.add_argument(
        '--rnaseq-output',
        type=Path,
        required=True,
        help='RNA-seq pipeline output directory'
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
        print(f"Project: {agent.project_id}")
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
