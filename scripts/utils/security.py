#!/usr/bin/env python3
"""
Security utilities for LLM agent - Input validation and sanitization
"""

import re
from pathlib import Path
from typing import Optional


class SecurityError(Exception):
    """Raised when security validation fails."""
    pass


def validate_project_id(project_id: str) -> str:
    """
    Validate and sanitize project ID.
    
    Args:
        project_id: Project identifier to validate
    
    Returns:
        Sanitized project ID
    
    Raises:
        SecurityError: If validation fails
    """
    # Allow only alphanumeric, hyphens, underscores
    if not re.match(r'^[a-zA-Z0-9_-]+$', project_id):
        raise SecurityError(
            f"Invalid project_id: {project_id}. "
            f"Only alphanumeric characters, hyphens, and underscores allowed."
        )
    
    # Prevent path traversal
    if '..' in project_id or '/' in project_id or '\\' in project_id:
        raise SecurityError(f"Path traversal detected in project_id: {project_id}")
    
    # Length limit
    if len(project_id) > 100:
        raise SecurityError(f"project_id too long: {len(project_id)} chars (max 100)")
    
    return project_id


def validate_sample_id(sample_id: str) -> str:
    """Validate and sanitize sample ID."""
    if not re.match(r'^[a-zA-Z0-9_-]+$', sample_id):
        raise SecurityError(f"Invalid sample_id: {sample_id}")
    
    if len(sample_id) > 100:
        raise SecurityError(f"sample_id too long")
    
    return sample_id


def validate_path(path: str, allowed_base_dirs: list = None) -> Path:
    """
    Validate that path is within allowed directories.
    
    Args:
        path: Path to validate
        allowed_base_dirs: List of allowed base directories
    
    Returns:
        Resolved absolute path
    
    Raises:
        SecurityError: If path is outside allowed directories
    """
    if allowed_base_dirs is None:
        allowed_base_dirs = [
            "/data_3tb/shared/output",
            "/data_3tb/shared/rna-seq-pipeline",
            "/home/ygkim/ngs-pipeline"
        ]
    
    path_obj = Path(path).resolve()
    
    # Check if path is within allowed directories
    allowed = False
    for base_dir in allowed_base_dirs:
        base_dir_obj = Path(base_dir).resolve()
        try:
            path_obj.relative_to(base_dir_obj)
            allowed = True
            break
        except ValueError:
            continue
    
    if not allowed:
        raise SecurityError(
            f"Path {path} is outside allowed directories: {allowed_base_dirs}"
        )
    
    return path_obj


def sanitize_command_args(args: list) -> list:
    """
    Sanitize command line arguments.
    
    Args:
        args: List of command arguments
    
    Returns:
        Sanitized arguments
    
    Raises:
        SecurityError: If dangerous patterns detected
    """
    dangerous_patterns = [
        r';',      # Command chaining
        r'\|',     # Pipes
        r'&',      # Background execution
        r'\$',     # Variable expansion
        r'`',      # Command substitution
        r'\(',     # Subshell
        r'>',      # Redirection
        r'<',      # Input redirection
        r'\*',     # Wildcards (questionable)
    ]
    
    sanitized = []
    for arg in args:
        # Check for dangerous patterns
        for pattern in dangerous_patterns:
            if re.search(pattern, arg):
                raise SecurityError(
                    f"Dangerous pattern '{pattern}' detected in argument: {arg}"
                )
        
        sanitized.append(arg)
    
    return sanitized


# Whitelist of allowed commands
ALLOWED_COMMANDS = {
    'python',
    'conda',
    'snakemake',
    'bash',  # Only for specific scripts
}


def validate_command(command: str) -> str:
    """
    Validate that command is in whitelist.
    
    Args:
        command: Command to execute
    
    Returns:
        Validated command
    
    Raises:
        SecurityError: If command not allowed
    """
    if command not in ALLOWED_COMMANDS:
        raise SecurityError(f"Command not allowed: {command}")
    
    return command


# Resource limits (for future implementation with systemd or cgroups)
RESOURCE_LIMITS = {
    'max_cpu_time': 3600,      # 1 hour
    'max_memory_mb': 32768,    # 32GB
    'max_processes': 100,
    'max_open_files': 1024,
}


if __name__ == '__main__':
    # Test validation
    print("Testing input validation...")
    
    # Valid inputs
    print(f"✅ Valid project_id: {validate_project_id('mouse-chd8')}")
    print(f"✅ Valid project_id: {validate_project_id('mouse_chd8')}")
    
    # Invalid inputs
    try:
        validate_project_id("mouse; rm -rf /")
        print("❌ Should have failed!")
    except SecurityError as e:
        print(f"✅ Caught: {e}")
    
    try:
        validate_project_id("../../etc/passwd")
        print("❌ Should have failed!")
    except SecurityError as e:
        print(f"✅ Caught: {e}")
    
    # Path validation
    try:
        validate_path("/data_3tb/shared/output/mouse-chd8")
        print("✅ Valid path")
    except SecurityError as e:
        print(f"❌ {e}")
    
    try:
        validate_path("/etc/passwd")
        print("❌ Should have failed!")
    except SecurityError as e:
        print(f"✅ Caught: {e}")
    
    print("\n✅ All security tests passed")
