# Phase 4: E2E Testing Guide

**Status**: In Progress  
**Created**: 2026-02-26  
**Objective**: Validate complete automation flow from natural language to DE/GO preparation

---

## Overview

Phase 4 tests the entire automation pipeline end-to-end:

```
User Natural Language
    ↓
LLM Agent (Ollama)
    ↓
Tool Execution (auto_config.py)
    ↓
Bridge Script (bridge_to_de_pipeline.py)
    ↓
DE/GO Preparation (files ready)
```

---

## Test Scripts

### 1. Security Testing

**File**: `scripts/test_security.py`

**Purpose**: Validate input validation and security features

**Run**:
```bash
cd /data_3tb/shared/rna-seq-pipeline
python scripts/test_security.py
```

**Tests**:
- ✅ Project ID validation (injection prevention)
- ✅ Sample ID validation
- ✅ Path validation (whitelist enforcement)
- ✅ Command sanitization (dangerous pattern blocking)
- ✅ Command whitelist validation

**Expected Output**:
```
TEST SUMMARY
============================================================
✅ PASS: Project ID Validation
✅ PASS: Sample ID Validation
✅ PASS: Path Validation
✅ PASS: Command Sanitization
✅ PASS: Command Whitelist

Result: 5/5 tests passed
ALL TESTS PASSED ✅
```

---

### 2. E2E Automation Testing

**File**: `scripts/test_e2e_automation.py`

**Purpose**: Test complete automation flow with real data

**Run**:
```bash
cd /data_3tb/shared/rna-seq-pipeline

# Full test (creates real files)
python scripts/test_e2e_automation.py --project-id mouse-chd8

# Dry run (skip file creation)
python scripts/test_e2e_automation.py --project-id mouse-chd8 --dry-run
```

**Test Phases**:

#### Phase 1: Auto-Config Generation
- Remove existing config for clean test
- Call `ensure_bridge_config(project_id, force=True)`
- Verify config file created
- Check YAML structure (project_id, rnaseq_output, de_pipeline, counts_relpath)
- Display discovered paths

**Success Criteria**:
- Config file: `config/projects/paths_mouse_chd8.yaml`
- All required keys present
- Paths auto-discovered

#### Phase 2: Path Validation
- Load generated config
- Call `validate_paths(config)`
- Check each required path exists:
  - RNA-seq output directory
  - DE pipeline directory
  - Counts matrix file
  - Sample sheet file

**Success Criteria**:
- All paths exist: True
- Counts file accessible
- Sample sheet found

#### Phase 3: Bridge Execution
- Run bridge script with auto-generated config
- Command: `conda run -n rna-seq-pipeline python scripts/bridge_to_de_pipeline.py --project-id mouse-chd8 --skip-de --yes`
- Capture stdout/stderr
- Parse output for completion messages

**Success Criteria**:
- Return code: 0
- QC summary displayed
- Counts matrix copied
- Metadata generated
- DE config created

#### Phase 4: Output Verification
- Check files in DE pipeline:
  - `data/raw/mouse-chd8_counts.csv`
  - `data/raw/mouse-chd8_metadata.csv`
  - `configs/config_mouse-chd8.yml`
- Verify file contents:
  - Counts matrix has data (rows × columns)
  - Metadata has samples
  - Config has comparisons

**Success Criteria**:
- All 3 files exist
- File sizes > 0
- Content validation passes

---

### 3. LLM Agent Testing (Optional)

**Prerequisite**: Ollama installed and running

**Setup**:
```bash
# Install Ollama (if not installed)
curl -fsSL https://ollama.com/install.sh | sh

# Pull model
ollama pull llama3.1:8b

# Start Ollama server (separate terminal)
ollama serve
```

**Run Agent**:
```bash
cd /data_3tb/shared/rna-seq-pipeline

python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8
```

**Natural Language Test Queries**:

1. **DE Preparation**:
   ```
   User: mouse-chd8 프로젝트 DE 분석 준비해줘
   Expected: prepare_de_analysis() called → bridge script runs → files created
   ```

2. **Path Validation**:
   ```
   User: 경로 검증해줘
   Expected: validate_paths() called → all paths checked → status reported
   ```

3. **Config Check**:
   ```
   User: bridge 설정 확인해줘
   Expected: check_bridge_config() called → config status returned
   ```

4. **Sample Details**:
   ```
   User: 현재 QC 상태 보여줘
   Expected: get_sample_details() called → sample info displayed
   ```

---

## Test Matrix

| Test Type | Test Name | Status | Notes |
|-----------|-----------|--------|-------|
| **Security** | Project ID validation | ✅ Ready | Injection prevention |
| | Sample ID validation | ✅ Ready | Input sanitization |
| | Path validation | ✅ Ready | Whitelist enforcement |
| | Command sanitization | ✅ Ready | Pattern blocking |
| | Command whitelist | ✅ Ready | Allowed commands only |
| **E2E** | Auto-config generation | ✅ Ready | Phase 1 |
| | Path validation | ✅ Ready | Phase 2 |
| | Bridge execution | ✅ Ready | Phase 3 |
| | Output verification | ✅ Ready | Phase 4 |
| **Integration** | LLM agent → prepare_de | ⏳ Pending | Requires Ollama |
| | LLM agent → validate_paths | ⏳ Pending | Requires Ollama |
| | LLM agent → check_config | ⏳ Pending | Requires Ollama |

---

## Expected Workflow

### Scenario 1: Fresh Project (No Config)

**User Action**:
```
사용자: mouse-chd8 프로젝트 DE 분석 준비해줘
```

**System Flow**:
1. LLM agent interprets query → calls `prepare_de_analysis("mouse-chd8")`
2. Tool checks for config → not found
3. Calls `ensure_bridge_config("mouse-chd8", force=False)`
4. Auto-config discovers paths from `project_summary.json`
5. Generates `config/projects/paths_mouse_chd8.yaml`
6. Bridge script runs with auto-generated config
7. Files created in DE pipeline
8. Agent responds: "✅ DE 분석 준비 완료!"

**Verification**:
```bash
# Check config created
ls config/projects/paths_mouse_chd8.yaml

# Check DE files
ls /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis/data/raw/mouse-chd8_*
ls /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis/configs/config_mouse-chd8.yml
```

---

### Scenario 2: Existing Config (Reuse)

**User Action**:
```
사용자: 경로 검증해줘
```

**System Flow**:
1. LLM agent calls `validate_paths("mouse-chd8")`
2. Tool loads existing config
3. Checks each path exists
4. Returns validation results
5. Agent responds with path status

**Expected Output**:
```json
{
  "rnaseq_output": true,
  "de_pipeline": true,
  "counts_matrix": true,
  "sample_sheet": true
}
```

---

### Scenario 3: Security Validation

**Malicious Input**:
```
사용자: mouse-chd8; rm -rf / 프로젝트 정보 보여줘
```

**System Response**:
1. Agent extracts project_id: `"mouse-chd8; rm -rf /"`
2. Tool calls `validate_project_id(project_id)`
3. Security layer detects injection pattern
4. Raises `SecurityError("Invalid project_id: contains dangerous patterns")`
5. Agent catches error
6. Responds: "❌ Security validation failed: Invalid project ID format"

**Expected Behavior**:
- No commands executed
- Error logged
- User notified of validation failure

---

## Success Criteria

### Security Tests
- ✅ All 5 test suites pass
- ✅ Injection attempts blocked
- ✅ Path traversal prevented
- ✅ Dangerous patterns detected

### E2E Tests
- ✅ Config auto-generated (Phase 1)
- ✅ Paths validated (Phase 2)
- ✅ Bridge executed successfully (Phase 3)
- ✅ Output files created (Phase 4)

### Integration Tests (if Ollama available)
- ✅ Natural language queries understood
- ✅ Correct tools called
- ✅ Results returned to user
- ✅ Security validation applied

---

## Troubleshooting

### Test Failure: Config Not Generated

**Symptom**: Phase 1 fails with "Config generation failed"

**Possible Causes**:
1. `project_summary.json` not found
2. Invalid JSON format
3. Missing required paths

**Debug**:
```bash
# Check if project_summary.json exists
ls /data_3tb/shared/output/mouse-chd8/project_summary.json

# Validate JSON
cat /data_3tb/shared/output/mouse-chd8/project_summary.json | python -m json.tool

# Test auto-config manually
cd /data_3tb/shared/rna-seq-pipeline
python scripts/utils/auto_config.py --project-id mouse-chd8
```

---

### Test Failure: Paths Not Validated

**Symptom**: Phase 2 fails with "Missing paths"

**Possible Causes**:
1. RNA-seq output directory doesn't exist
2. Counts file not generated yet
3. Sample sheet filename mismatch

**Debug**:
```bash
# Check RNA-seq output
ls /data_3tb/shared/output/mouse-chd8/

# Check counts file
ls /data_3tb/shared/output/mouse-chd8/project_summary/counts/counts_matrix_clean.csv

# Check sample sheet (try both variants)
ls /data_3tb/shared/rna-seq-pipeline/config/samples/mouse-chd8.tsv
ls /data_3tb/shared/rna-seq-pipeline/config/samples/mouse_chd8.tsv
```

---

### Test Failure: Bridge Execution Error

**Symptom**: Phase 3 fails with non-zero return code

**Possible Causes**:
1. Conda environment not activated
2. Bridge script has bugs
3. Permissions issue

**Debug**:
```bash
# Test conda environment
conda run -n rna-seq-pipeline python --version

# Run bridge manually with verbose output
cd /data_3tb/shared/rna-seq-pipeline
conda run -n rna-seq-pipeline python scripts/bridge_to_de_pipeline.py \
  --project-id mouse-chd8 \
  --skip-de \
  --yes \
  2>&1 | tee bridge_debug.log
```

---

### Test Failure: Output Files Missing

**Symptom**: Phase 4 fails with "Files not found"

**Possible Causes**:
1. Bridge script didn't copy files
2. Wrong output directory
3. Insufficient permissions

**Debug**:
```bash
# Check DE pipeline directory
ls -la /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis/data/raw/

# Check if bridge completed
grep "✅" bridge_debug.log

# Check permissions
ls -ld /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis/data/raw/
```

---

## Next Steps After Testing

### If All Tests Pass ✅

1. **Update Documentation**:
   - Mark Phase 4 as complete
   - Document test results
   - Update README with automation instructions

2. **Production Deployment**:
   - Deploy security hardening (systemd/Docker)
   - Set up Ollama service
   - Configure resource limits
   - Enable audit logging

3. **User Training**:
   - Create user guide for natural language queries
   - Document common workflows
   - Provide troubleshooting tips

### If Tests Fail ❌

1. **Debug Failed Tests**:
   - Run test scripts individually
   - Enable verbose logging
   - Check error messages
   - Review recent changes

2. **Fix Issues**:
   - Update code based on test failures
   - Re-run tests after fixes
   - Document bug fixes

3. **Re-test**:
   - Run security tests first
   - Then E2E tests
   - Finally integration tests

---

## Test Log Template

```markdown
## Test Run: YYYY-MM-DD HH:MM

### Environment
- Server: /data_3tb/shared/rna-seq-pipeline
- Project: mouse-chd8
- Conda Env: rna-seq-pipeline

### Results

#### Security Tests
- [ ] Project ID validation: PASS/FAIL
- [ ] Sample ID validation: PASS/FAIL
- [ ] Path validation: PASS/FAIL
- [ ] Command sanitization: PASS/FAIL
- [ ] Command whitelist: PASS/FAIL

#### E2E Tests
- [ ] Phase 1 (Auto-config): PASS/FAIL
- [ ] Phase 2 (Path validation): PASS/FAIL
- [ ] Phase 3 (Bridge execution): PASS/FAIL
- [ ] Phase 4 (Output verification): PASS/FAIL

#### Integration Tests
- [ ] LLM agent → prepare_de: PASS/FAIL/SKIP
- [ ] LLM agent → validate_paths: PASS/FAIL/SKIP
- [ ] LLM agent → check_config: PASS/FAIL/SKIP

### Issues Encountered
- None / [List issues]

### Notes
- [Any relevant observations]
```

---

## References

- **Security Guide**: `docs/developer/SECURITY_GUIDE.md`
- **LLM Agent Tools**: `docs/developer/PHASE3_LLM_AGENT_TOOLS.md`
- **Auto-Config Examples**: `scripts/AUTOCONFIG_EXAMPLES.py`
- **Bridge Testing**: `scripts/test_bridge_autoconfig.sh`
