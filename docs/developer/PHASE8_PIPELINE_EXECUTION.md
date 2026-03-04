# Phase 8: Pipeline Execution Tools

**Objective**: Enable agent to execute full RNA-seq pipeline from raw FASTQ  
**Status**: ✅ Phase 8A COMPLETE & TESTED | ✅ Phase 8B COMPLETE & TESTED  
**Date**: 2026-03-04 (8A) / 2026-03-05 (8B)

## Testing Results — Phase 8A (2026-03-04)

**Environment**: mouse-chd8, qwen2.5:32b, Snakemake 7+

| Step | Query | Tool Called | Result |
|------|-------|-------------|--------|
| 1 | "/data_3tb/.../fastq/ 폴더에서 FASTQ 파일 찾아줘" | `detect_fastq_files` | ✅ 76파일, 38샘플, 149GB, paired-end 정확히 감지 |
| 2 | "새 프로젝트 만들어줘. ID는 test-2026c ..." | `create_project_config` | ✅ `config/projects/test-2026c.yaml` 생성 성공 |
| 3 | "config/projects/test-2026.yaml 로 입력 데이터 검증해줘" | `validate_input_data` | ✅ 데이터/디스크 OK, 참조게놈 경로 오류 정확히 감지 |
| 4 | "config/projects/test-2026.yaml 로 파이프라인 dry-run 해줘" | `run_pipeline(dry_run=True)` | ✅ 616 jobs, 11개 rule 목록 정확히 파싱 및 설명 |

**Dry-run 상세** (38 samples, Snakemake 7):
```
all: 1 / convert_counts_matrix: 1 / copy_bam_to_standard: 38
cutadapt: 38 / fastqc_raw: 76 / featurecounts_quant: 1
generate_manifest: 38 / generate_qc_summary: 38
index_bam: 38 / multiqc: 1 / star_align: 38  →  total: 270
```

**Bugs fixed during testing**:
| Bug | Fix |
|-----|-----|
| `create_project_config` 응답에 `TOOL_CALL:` JSON 코드블록 노출 | `_build_system_prompt(native_tools=True)` — TOOL_CALL 예시 제거 |
| `create_project_config` 응답 끝에 다음 tool 제안 코드블록 삽입 | summary 요청 user turn에 "코드블록 포함하지 말 것" 명시 |
| `run_pipeline` dry-run 결과를 오류로 오해석 | Snakemake 출력 파싱 → `{total_jobs, jobs_by_rule, note}` 구조체 반환 |
| Snakemake `Job stats:` 파싱 실패 | 실제 출력 형식(`job / count` 헤더) 확인 후 정규식 수정 |
| dry-run "missing output" 줄을 오류로 분류 | `issues` 필터에서 `missing` 제외, `note` 필드로 정상 동작 설명 |

---

## Testing Results — Phase 8B (2026-03-05)

**Environment**: mouse-chd8, qwen2.5:32b, 38 samples

| Step | Query | Tool Called | Result |
|------|-------|-------------|--------|
| 1 | "리소스 얼마나 필요해?" | `estimate_resources` | ✅ 149GB in → 746GB out, 5-10h, RAM 221GB/Disk 2145GB 여유 확인 |
| 2 | "샘플 시트 만들어줘" | `create_sample_sheet` | ✅ 38 샘플, wildtype 20 / heterozygous 18 자동 분류, TSV 생성 |
| 3 | "mouse-chd8 파이프라인 상태 보여줘" | `monitor_pipeline` | ✅ 실행 상태 정확히 반환 (not_started 확인) |

**`estimate_resources` 상세**:
- 입력: 149GB / 예상 출력: 746GB (×5 추정)
- 예상 소요 시간: 5-10시간 (8 cores)
- STAR 메모리: ~30GB (222GB 중 사용)
- 디스크 여유: 2145GB ✅ (충분)

**`create_sample_sheet` 자동 분류 로직**:
- 파일명 `_W_` 또는 `_WT_` 패턴 → `wildtype`
- 파일명 `_H_` 또는 `_HET_` 패턴 → `heterozygous`
- 파일명 `_KO_` 패턴 → `knockout`
- 분류 불가 → `unknown` (수동 지정 필요)

**`monitor_pipeline` 감지 항목**:
```
manifest.json   → generate_manifest 완료 수
*.bam           → star_align 완료 수
logs/cutadapt/  → cutadapt 완료 수
*_fastqc.zip    → fastqc_raw 완료 수
multiqc_report  → multiqc 완료 수
featureCounts   → featurecounts_quant 완료 수
```

**Bugs fixed during Phase 8B testing**:
| Bug | Fix |
|-----|-----|
| `monitor_pipeline` "결과 디렉토리 없음" 반환 (경로 이중 중첩) | `base_results_dir`에 이미 project_id 포함된 경우 fallback: sample dirs 존재 시 base 자체를 `project_dir`로 사용 |

---

## New Tools to Add

### 1. `create_project_config()`
**Purpose**: Generate config.yaml from natural language  
**Parameters**:
- `project_id`: str
- `data_dir`: str (FASTQ location)
- `results_dir`: str (output location)
- `species`: str (human/mouse/...)
- `read_type`: str (paired-end/single-end)

**Example**:
```
You> "새 프로젝트 만들어줘. ID는 test-2026, FASTQ는 /data/raw/test/에 있어"
Agent> TOOL_CALL: {
  "name": "create_project_config",
  "parameters": {
    "project_id": "test-2026",
    "data_dir": "/data/raw/test",
    "results_dir": "/data/output/test-2026",
    "species": "human"  # 추론 또는 질문
  }
}
```

---

### 2. `detect_fastq_files()`
**Purpose**: Scan directory and identify FASTQ files  
**Parameters**:
- `data_dir`: str
- `pattern`: str (optional, e.g., "*_R{1,2}.fastq.gz")

**Returns**:
```json
{
  "total_files": 24,
  "samples": [
    {
      "sample_id": "Sample1",
      "R1": "/data/raw/test/Sample1_R1.fastq.gz",
      "R2": "/data/raw/test/Sample1_R2.fastq.gz",
      "size_gb": 2.5
    }
  ],
  "read_type": "paired-end",
  "estimated_storage": "150 GB"
}
```

---

### 3. `create_sample_sheet()`
**Purpose**: Generate sample metadata TSV  
**Parameters**:
- `project_id`: str
- `samples`: list[str] (from detect_fastq_files)
- `conditions`: dict[str, list[str]] (manual or auto-group)

**Example**:
```
You> "샘플들을 Ctrl과 Treat으로 나눠줘. 1-6은 Ctrl, 7-12는 Treat"
Agent> TOOL_CALL: {
  "name": "create_sample_sheet",
  "parameters": {
    "project_id": "test-2026",
    "samples": ["Sample1", "Sample2", ...],
    "conditions": {
      "Ctrl": ["Sample1", "Sample2", ..., "Sample6"],
      "Treat": ["Sample7", "Sample8", ..., "Sample12"]
    }
  }
}
```

---

### 4. `run_pipeline()`
**Purpose**: Execute Snakemake workflow  
**Parameters**:
- `config_file`: str
- `cores`: int (default: 8)
- `dry_run`: bool (default: True for safety)
- `until_rule`: str (optional, e.g., "fastqc", "multiqc")

**Example**:
```
You> "파이프라인 실행해줘"
Agent> TOOL_CALL: {
  "name": "run_pipeline",
  "parameters": {
    "config_file": "config/projects/test-2026.yaml",
    "cores": 16,
    "dry_run": false
  }
}
```

**Security**: 
- Require explicit confirmation for actual runs
- Show dry-run results first
- Estimate runtime and resource usage

---

### 5. `monitor_pipeline()`
**Purpose**: Check pipeline execution status  
**Parameters**:
- `project_id`: str

**Returns**:
```json
{
  "status": "running",
  "progress": "35%",
  "current_rule": "star_align",
  "current_sample": "Sample8",
  "completed_samples": 7,
  "total_samples": 12,
  "estimated_remaining": "2h 15m",
  "errors": []
}
```

---

### 6. `validate_input_data()`
**Purpose**: Pre-flight checks before pipeline execution  
**Parameters**:
- `config_file`: str

**Checks**:
- FASTQ files exist and readable
- Sufficient disk space
- Reference genome available
- STAR index built
- Sample sheet valid

**Returns**:
```json
{
  "valid": true,
  "warnings": [
    "Disk space: 150GB required, 200GB available (OK)"
  ],
  "errors": [],
  "estimated_runtime": "6-8 hours (12 samples, 16 cores)"
}
```

---

## Workflow Example

### Scenario: Brand New Project

```
You> "새 RNA-seq 프로젝트 시작하고 싶어. FASTQ 파일은 /data/raw/exp123/에 있어"

Agent> 알겠습니다! 먼저 FASTQ 파일들을 확인해볼게요.
🔧 Calling tool: detect_fastq_files
📋 Parameters: {"data_dir": "/data/raw/exp123/"}

Agent> 12개 샘플을 찾았습니다 (paired-end, 총 150GB).
       프로젝트 ID는 뭘로 할까요?

You> "exp123-human-2026"

Agent> 좋아요! 결과는 어디에 저장할까요?

You> "/data/output/exp123/"

Agent> 샘플들의 실험 조건은 어떻게 되나요? 
       예: Control vs Treatment

You> "1-6은 WT, 7-12는 KO야"

Agent> 
🔧 Calling tool: create_project_config
🔧 Calling tool: create_sample_sheet

✅ Config 생성 완료: config/projects/exp123-human-2026.yaml
✅ Sample sheet 생성: config/samples/exp123-human-2026.tsv

데이터 검증 중...
🔧 Calling tool: validate_input_data

✅ 모든 검사 통과
   - FASTQ 파일: 12 samples, 150GB
   - 디스크 공간: 충분 (200GB 사용 가능)
   - Reference: human (GRCh38) 준비됨
   - 예상 실행 시간: 6-8시간 (16 cores)

파이프라인을 실행할까요? (실제 실행 전 dry-run을 보여드릴 수 있어요)

You> "dry-run 먼저 보여줘"

Agent>
🔧 Calling tool: run_pipeline
📋 Parameters: {"dry_run": true}

Dry-run 결과:
- 총 작업: 156 jobs
- FastQC: 24 jobs
- Trimming: 12 jobs  
- Alignment: 12 jobs
- Counting: 12 jobs
- QC report: 1 job

문제없어 보입니다. 실제 실행할까요?

You> "응"

Agent>
🔧 Calling tool: run_pipeline
📋 Parameters: {"dry_run": false, "cores": 16}

✅ 파이프라인 시작!
   진행 상황은 "파이프라인 상태 보여줘" 로 확인할 수 있어요.

[6시간 후...]

You> "파이프라인 상태 보여줘"

Agent>
🔧 Calling tool: monitor_pipeline

✅ 파이프라인 완료! (100%)
   - 12 samples 처리 완료
   - QC pass rate: 100%
   - 결과 위치: /data/output/exp123/

DE 분석도 준비할까요?

You> "응"

Agent>
🔧 Calling tool: prepare_de_analysis

✅ DE/GO 분석 준비 완료!
   다음 명령어로 실행하세요:
   cd /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis
   snakemake --cores 8
```

---

## Implementation Priority

### Phase 8A (Essential)
1. ✅ `create_project_config()` - Config 생성
2. ✅ `detect_fastq_files()` - FASTQ 감지
3. ✅ `validate_input_data()` - 사전 검증
4. ✅ `run_pipeline()` - 파이프라인 실행

### Phase 8B (Nice to Have)
5. ✅ `monitor_pipeline()` - 실행 모니터링
6. ✅ `create_sample_sheet()` - Sample sheet 생성
7. ✅ `estimate_resources()` - 리소스 예측

---

## Technical Considerations

### Config Generation
```python
def create_project_config(project_id, data_dir, results_dir, species="human"):
    """Generate config.yaml from parameters."""
    
    config = {
        'project_id': project_id,
        'use_standard_structure': True,
        'data_dir': data_dir,
        'base_results_dir': results_dir,
        'species': species,
        # ... other defaults
    }
    
    config_path = f"config/projects/{project_id}.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    return {"status": "success", "config_path": config_path}
```

### Pipeline Execution
```python
def run_pipeline(config_file, cores=8, dry_run=True, until_rule=None):
    """Execute Snakemake workflow."""
    
    cmd = [
        "snakemake",
        "--configfile", config_file,
        "--cores", str(cores)
    ]
    
    if dry_run:
        cmd.append("--dry-run")
    
    if until_rule:
        cmd.extend(["--until", until_rule])
    
    # Security: validate config path
    validate_path(config_file, ALLOWED_DIRS)
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    return {
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr
    }
```

### Progress Monitoring
```python
def monitor_pipeline(project_id):
    """Check pipeline execution status."""
    
    # Parse Snakemake logs
    log_file = f"logs/{project_id}/snakemake.log"
    
    # Extract progress from log
    total_jobs = count_total_jobs(log_file)
    completed_jobs = count_completed_jobs(log_file)
    
    progress = (completed_jobs / total_jobs) * 100 if total_jobs > 0 else 0
    
    return {
        "status": "running" if progress < 100 else "completed",
        "progress": f"{progress:.1f}%",
        "completed_jobs": completed_jobs,
        "total_jobs": total_jobs
    }
```

---

## Security & Safety

### Confirmation Required
- Always dry-run before actual execution
- Show resource estimates
- Require explicit user confirmation
- Log all pipeline executions

### Resource Limits
- Max cores: 32 (configurable)
- Disk space check before start
- Timeout for long-running jobs
- Auto-pause if errors exceed threshold

### Input Validation
- Validate all file paths
- Check FASTQ file integrity
- Verify config parameters
- Sanitize project IDs

---

## Testing Strategy

### Unit Tests
- Config generation with various parameters
- FASTQ detection with different naming patterns
- Pipeline execution (dry-run only in tests)

### Integration Tests
- Full workflow: config → detect → validate → run
- Error handling: missing files, invalid config
- Resource estimation accuracy

### E2E Tests
- Real small dataset (2-3 samples)
- Complete pipeline execution
- Verify all outputs generated

---

## Next Steps

1. **Phase 6 완료**: 현재 agent 테스트
2. **Phase 7**: Production 배포 가이드
3. **Phase 8A 구현**: 위의 Essential tools 추가
4. **Phase 8B 테스트**: Small dataset으로 전체 워크플로우
5. **Phase 9**: 데이터 타입 감지 (ATAC-seq, WGS 확장)
6. **Phase 10**: Unified agent 통합

---

**Timeline Estimate**:
- Phase 6-7: 1-2 weeks (현재 agent 안정화)
- Phase 8A: 2-3 weeks (도구 구현 + 테스트)
- Phase 8B: 1-2 weeks (통합 + 검증)
- Phase 9-10: 4-6 weeks (다른 파이프라인 확장)

**Total**: ~3 months to full Option C implementation
