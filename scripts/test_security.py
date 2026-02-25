#!/usr/bin/env python3
"""
Security Layer Testing

Tests input validation and security features:
1. Project ID validation (injection prevention)
2. Path validation (traversal prevention)
3. Command sanitization
4. Error handling

Run:
    python scripts/test_security.py
"""

import sys
from pathlib import Path

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.security import (
    validate_project_id,
    validate_sample_id,
    validate_path,
    sanitize_command_args,
    validate_command,
    SecurityError
)


class Colors:
    """Terminal colors for pretty output."""
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    CYAN = '\033[96m'


def print_test(name):
    print(f"\n{Colors.BOLD}{Colors.CYAN}TEST: {name}{Colors.ENDC}")


def print_pass(msg):
    print(f"  {Colors.OKGREEN}✅ PASS: {msg}{Colors.ENDC}")


def print_fail(msg):
    print(f"  {Colors.FAIL}❌ FAIL: {msg}{Colors.ENDC}")


def test_project_id_validation():
    """Test project_id validation with various inputs."""
    print_test("Project ID Validation")
    
    # Valid cases
    valid_ids = [
        "mouse-chd8",
        "human_H2O2",
        "project123",
        "test-project_v2",
    ]
    
    for pid in valid_ids:
        try:
            validate_project_id(pid)
            print_pass(f"Valid ID accepted: '{pid}'")
        except SecurityError as e:
            print_fail(f"Valid ID rejected: '{pid}' - {e}")
            return False
    
    # Invalid cases (should raise SecurityError)
    invalid_ids = [
        ("mouse-chd8; rm -rf /", "command injection"),
        ("../etc/passwd", "path traversal"),
        ("project$(whoami)", "command substitution"),
        ("test`id`", "backtick execution"),
        ("a" * 101, "too long"),
        ("", "empty string"),
        ("project|cat /etc/passwd", "pipe injection"),
    ]
    
    for pid, reason in invalid_ids:
        try:
            validate_project_id(pid)
            print_fail(f"Invalid ID accepted: '{pid}' ({reason})")
            return False
        except SecurityError as e:
            print_pass(f"Invalid ID rejected: '{pid[:50]}...' ({reason})")
    
    return True


def test_sample_id_validation():
    """Test sample_id validation."""
    print_test("Sample ID Validation")
    
    # Valid cases
    valid_ids = [
        "Ctrl_1",
        "Treat-1",
        "Sample123",
        "exp_v2",
    ]
    
    for sid in valid_ids:
        try:
            validate_sample_id(sid)
            print_pass(f"Valid sample ID: '{sid}'")
        except SecurityError as e:
            print_fail(f"Valid sample ID rejected: '{sid}' - {e}")
            return False
    
    # Invalid cases
    invalid_ids = [
        ("Sample; rm -rf /", "command injection"),
        ("../../../etc/passwd", "path traversal"),
        ("sample$(whoami)", "command substitution"),
    ]
    
    for sid, reason in invalid_ids:
        try:
            validate_sample_id(sid)
            print_fail(f"Invalid sample ID accepted: '{sid}' ({reason})")
            return False
        except SecurityError as e:
            print_pass(f"Invalid sample ID rejected: '{sid}' ({reason})")
    
    return True


def test_path_validation():
    """Test path validation with whitelist."""
    print_test("Path Validation")
    
    allowed_dirs = [
        "/data_3tb/shared/output",
        "/home/ygkim/ngs-pipeline",
    ]
    
    # Valid paths
    valid_paths = [
        "/data_3tb/shared/output/mouse-chd8",
        "/home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis",
        "/data_3tb/shared/output/human_H2O2/results",
    ]
    
    for path in valid_paths:
        try:
            validate_path(path, allowed_dirs)
            print_pass(f"Valid path: {path}")
        except SecurityError as e:
            print_fail(f"Valid path rejected: {path} - {e}")
            return False
    
    # Invalid paths
    invalid_paths = [
        ("/etc/passwd", "outside whitelist"),
        ("/tmp/evil", "outside whitelist"),
        ("../../etc/passwd", "relative path"),
        ("/data_3tb/shared/output/../../../etc/passwd", "path traversal"),
    ]
    
    for path, reason in invalid_paths:
        try:
            validate_path(path, allowed_dirs)
            print_fail(f"Invalid path accepted: {path} ({reason})")
            return False
        except SecurityError as e:
            print_pass(f"Invalid path rejected: {path} ({reason})")
    
    return True


def test_command_sanitization():
    """Test command argument sanitization."""
    print_test("Command Argument Sanitization")
    
    # Valid args
    valid_args = [
        ["--project-id", "mouse-chd8"],
        ["--output", "/data_3tb/shared/output"],
        ["--config", "config.yaml"],
    ]
    
    for args in valid_args:
        try:
            sanitize_command_args(args)
            print_pass(f"Valid args: {args}")
        except SecurityError as e:
            print_fail(f"Valid args rejected: {args} - {e}")
            return False
    
    # Invalid args (dangerous patterns)
    invalid_args = [
        (["--project-id", "mouse-chd8; rm -rf /"], "semicolon injection"),
        (["--output", "/tmp/$(whoami)"], "command substitution"),
        (["--config", "config.yaml && cat /etc/passwd"], "command chaining"),
        (["--id", "test`id`"], "backtick execution"),
        (["--path", "/data|grep secret"], "pipe injection"),
    ]
    
    for args, reason in invalid_args:
        try:
            sanitize_command_args(args)
            print_fail(f"Invalid args accepted: {args} ({reason})")
            return False
        except SecurityError as e:
            print_pass(f"Invalid args rejected: {str(args)[:50]}... ({reason})")
    
    return True


def test_command_validation():
    """Test command whitelist validation."""
    print_test("Command Whitelist Validation")
    
    # Valid commands
    valid_commands = [
        "python",
        "conda",
        "git",
        "snakemake",
    ]
    
    for cmd in valid_commands:
        try:
            validate_command(cmd)
            print_pass(f"Valid command: {cmd}")
        except SecurityError as e:
            print_fail(f"Valid command rejected: {cmd} - {e}")
            return False
    
    # Invalid commands
    invalid_commands = [
        ("rm", "dangerous command"),
        ("curl", "not whitelisted"),
        ("wget", "not whitelisted"),
        ("bash", "shell execution"),
    ]
    
    for cmd, reason in invalid_commands:
        try:
            validate_command(cmd)
            print_fail(f"Invalid command accepted: {cmd} ({reason})")
            return False
        except SecurityError as e:
            print_pass(f"Invalid command rejected: {cmd} ({reason})")
    
    return True


def run_all_tests():
    """Run complete security test suite."""
    print(f"\n{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.BOLD}SECURITY LAYER TESTING{Colors.ENDC}")
    print(f"{Colors.BOLD}{'='*60}{Colors.ENDC}")
    
    tests = [
        ("Project ID Validation", test_project_id_validation),
        ("Sample ID Validation", test_sample_id_validation),
        ("Path Validation", test_path_validation),
        ("Command Sanitization", test_command_sanitization),
        ("Command Whitelist", test_command_validation),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print_fail(f"Test crashed: {e}")
            results.append((name, False))
    
    # Summary
    print(f"\n{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.BOLD}TEST SUMMARY{Colors.ENDC}")
    print(f"{Colors.BOLD}{'='*60}{Colors.ENDC}\n")
    
    passed = sum(1 for _, p in results if p)
    total = len(results)
    
    for name, result in results:
        if result:
            print_pass(name)
        else:
            print_fail(name)
    
    print(f"\n{Colors.BOLD}Result: {passed}/{total} tests passed{Colors.ENDC}")
    
    if passed == total:
        print(f"{Colors.OKGREEN}{Colors.BOLD}ALL TESTS PASSED ✅{Colors.ENDC}\n")
        return True
    else:
        print(f"{Colors.FAIL}{Colors.BOLD}SOME TESTS FAILED ❌{Colors.ENDC}\n")
        return False


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
