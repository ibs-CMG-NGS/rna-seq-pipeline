# Phase 8: Pipeline Execution Tools

**Objective**: Enable agent to execute full RNA-seq pipeline from raw FASTQ
**Status**: ✅ Phase 8A COMPLETE & TESTED | ✅ Phase 8B COMPLETE & TESTED | ✅ Phase 9 (RSeQC Library Validation) COMPLETE | ✅ Phase 10 (DE-GO Pipeline Execution) COMPLETE
**Date**: 2026-03-04 (8A) / 2026-03-05 (8B) / 2026-03-11 (Phase 9) / 2026-03-13 (Phase 10)

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

---

## Phase 9: RSeQC Library Validation (2026-03-11)

**Objective**: Pre-flight library type verification before pipeline execution
**Status**: ✅ COMPLETE & TESTED (mouse-chd8, mouse-monSTIM-2026)

### New Tools

#### `validate_library_type(config_file, n_reads=500_000, cores=8, sample_id=None)`
RSeQC 기반 라이브러리 검증. 파이프라인 실행 전 자동 호출 권장.

**동작 순서:**
1. GTF→BED12 변환 (`_gtf_to_bed12()` 순수 Python 구현, `{genome_dir}/annotation_rseqc.bed`에 캐시)
2. 대표 샘플 선택 (`sample_id` 지정 또는 sample sheet 첫 번째)
3. FASTQ 서브셋 추출 (head -N lines, gz/uncompressed 모두 지원)
4. STAR 임시 정렬 → BAM index
5. `infer_experiment.py` → strandedness 감지 (0/1/2)
6. `read_distribution.py` → exon_pct, intergenic_pct 계산
7. strandedness가 config와 다르면 YAML 자동 업데이트
8. 임시 파일 정리

**판정 기준:**
- strandedness: failed>0.75 → 0(unstranded), 1++>0.6 → 1(forward), 1+-?>0.6 → 2(reverse)
- mRNA-seq: exon_pct>60% → OK, exon_pct<30% or intergenic>40% → 경고

**주의사항:**
- `gtf_to_gene_bed.py`는 RSeQC 5.x에서 제거됨 → `_gtf_to_bed12()` 순수 Python 대체 구현
- `read_distribution.py` 출력 파싱: col[2](Tag_count) 사용, col[1](Total_bases) 아님

#### `setup_and_validate(config_file, sample_id=None, n_reads=500_000, cores=8)`
detect_fastq_files + create_sample_sheet + validate_library_type 통합.
config에서 data_dir, project_id, sample_sheet 경로를 자동 추출.
반환: `{fastq_detection, sample_sheet, library_validation, summary(Korean)}`

**트리거 예시:** "프로젝트 설정해줘", "한번에 설정해줘", "setup"

### Testing Results (Phase 9)

| 프로젝트 | 샘플 | strandedness 감지 | exon_pct | 판정 | config 업데이트 |
|----------|------|-------------------|----------|------|-----------------|
| mouse-chd8 | Chd8_HPC_10M_W_1 | reverse (2) | ~71% | mRNA-seq ✅ | 0 → 2 ✅ |
| mouse-monSTIM | nMonTg4_S25 | - | 67% | mRNA-seq ✅ | - |
| mouse-monSTIM | MonTg1_S17 | - | 31.7% | suspicious ⚠️ | - (데이터 손상) |

### Bugs Fixed (Phase 9)

| Bug | Fix |
|-----|-----|
| `gtf_to_gene_bed.py` not found (RSeQC 5.x) | `_gtf_to_bed12()` 순수 Python 구현으로 대체 |
| STAR BAM not found (`Path / ""`) | `star_out_prefix = str(star_out_dir) + "/"` trailing slash 추가 |
| `read_distribution.py` 파싱 오류 (exon_pct 27714%) | `parts[2]` (Tag_count) 사용, `parts[1]` (Total_bases) 아님 |
| `fastqc_raw` mv-to-same-file 오류 | `[ src != dst ] && mv || true` 조건부 실행 |
| `monitor_pipeline` project_id 미설정 오류 | `_execute_tool`에서 config YAML 자동 파싱 |
| `run_pipeline` snakemake not found (VS Code) | 다중 후보 탐색 + subprocess PATH 주입 |
| `create_project_config` double path | `if _rp.name == project_id` 검사 추가 |

### VS Code Integration (2026-03-11)

`.vscode/tasks.json` 추가 — `conda run --no-capture-output -n rna-seq-pipeline` 방식:
- **RNA-seq Agent (interactive)** — `Ctrl+Shift+B` 기본 태스크
- **RNA-seq Agent (with existing project)** — config 경로 입력 프롬프트
- **RNA-seq Agent (reset session)** — `.agent_session.json` 초기화 후 시작

현재 모델: `qwen3:14b`

---

---

## Phase 10: DE-GO Pipeline Execution (2026-03-13)

**Objective**: Agent가 DE-GO Snakemake 파이프라인을 직접 실행/모니터링
**Status**: ✅ COMPLETE & TESTED (dry-run: 172 jobs, 16 rules 정상 파싱)

### New Tools

#### `run_de_pipeline(de_config_file, de_pipeline_dir=None, cores=8, dry_run=True, background=False)`
DE-GO Snakemake 파이프라인 실행. `run_pipeline()`과 동일한 패턴.

**핵심 동작:**
- DE Snakefile 위치: `de_pipeline_dir/Snakefile` (또는 `de_config_file`로부터 워크업)
- 명령어: `snakemake --snakefile <sf> --configfile <cfg> --use-conda --cores <N> --directory <dir>`
- `--use-conda` 필수: 모든 DE rules가 `conda: "rna-seq-de-go-analysis"` 사용
- `background=True`: `subprocess.Popen` + `logs/snakemake_<stem>.pid/log` 저장
- dry-run 파싱: `run_pipeline()`과 동일한 `Job stats:` 파서 재사용

**Session auto-resolve:**
- `de_config_file` 미제공 시 `self.current_de_config` 자동 주입
- `de_pipeline_dir` 미제공 시 `self.current_de_pipeline_dir` 자동 주입
- `run_bridge` 성공 후 두 값 모두 세션에 자동 저장됨 → 별도 지정 불필요

**권장 플로우 (agent 내):**
1. `run_de_pipeline(dry_run=True)` → jobs 목록 사용자에게 표시
2. 사용자 승인 후 `run_de_pipeline(dry_run=False, background=True)` → 백그라운드 실행

#### `monitor_de_pipeline(de_config_file=None, de_pipeline_dir=None)`
실행 중인 DE 파이프라인 진행상황 확인. STAR 로그 없어 `monitor_pipeline()`보다 단순.

**동작:**
- `_find_snakemake_log()` 재사용: `logs/snakemake_<stem>.log` → `.snakemake/log/` 순서로 탐색
- `_parse_snakemake_log()` 재사용: progress_pct, current_rule, error_rules, is_complete
- PID 파일 확인 (`os.kill(pid, 0)`): 프로세스 생존 여부로 running/completed/failed 판정
- 로그 마지막 20줄 `last_lines`로 반환 (에러 디버깅용)

**반환값:**
```python
{
    "status": "running"|"completed"|"failed"|"not_started",
    "progress_pct": float,
    "jobs_done": int,
    "jobs_total": int,
    "current_rule": str,
    "error_rules": [...],
    "log_file": str,
    "pid": int | None,
    "last_lines": str  # 최근 20줄
}
```

### Testing Results (Phase 10)

| Test | Result |
|------|--------|
| `run_de_pipeline(dry_run=True)` mouse-h2o2-astrocyte-2026 | ✅ `dry_run_ok`, 172 jobs, 16 rules 정상 파싱 |
| Tool schema `run_de_pipeline` in `_define_tools` | ✅ 확인 |
| Tool schema `monitor_de_pipeline` in `_define_tools` | ✅ 확인 |
| Session auto-resolve (`current_de_config`) | ✅ `run_bridge` 성공 시 자동 저장 |

**dry-run 상세** (4 pairwise comparisons, 16 samples):
```
run_pairwise_de: 4 / go_enrichment: 36 / kegg_enrichment: 12
generate_go_summary_table: 4 / go_barplots: 4
generate_pairwise_qc_plots: 4 / generate_pairwise_volcano: 4
generate_global_pca: 1 / generate_global_qc_plots: 1
aggregate_seqviewer: 1 / export_seqviewer_pair: 4
enrichment_done: 4 / run_omnibus_test: 1 / all: 1
total: 172 jobs
```

### Also Added in Phase 10 (earlier this session)

| Tool | Description |
|------|-------------|
| `validate_de_config_conditions` | DE config `pairwise_comparisons` 조건명 검증 (metadata 대비). 대소문자 차이, 알파벳 순서 오류 자동 감지 |
| `apply_de_config_corrections` | 검증 결과의 `suggestions`를 DE config YAML에 자동 적용 |
| `check_analysis_readiness` | DE 실행 전 5-point 사전 점검 (파일, 샘플 ID, 조건명, n per condition, zero-count %) |
| `read_de_results_summary` | `final_de_results.csv` + GO BP CSV 파싱 → DEG 수, 상위 유전자, GO 결과 요약 (인지 도구) |

### Bugs Fixed in Phase 10

| Bug | Fix |
|-----|-----|
| `run_pairwise_de` 조건명 불일치 (`1D`/`3D` vs `D1`/`D3`, `control` vs `Control`) | `validate_de_config_conditions` + `apply_de_config_corrections` 구현. `run_bridge` 후 자동 검증 |
| `generate_go_summary_table` KEGG 빈 CSV 오류 (`""` 내용) | `05_generate_go_table.R`: `readLines()` pre-check 추가 (GO와 동일 패턴) |
| Agent가 "DE 결과 요약해줘"에 `compare_conditions` 호출 | `required: []` 적용 + `current_de_config` 세션 저장 + system prompt CRITICAL TOOL ROUTING 추가 |
| `create_sample_sheet` overwrite 문제 (사용자 수동 편집 덮어쓰기) | `overwrite=False` 기본값으로 변경 |

---

## New Tools to Add (Archived — originally planned)

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

1. ✅ **Phase 6**: QC 분석 agent 완료
2. ✅ **Phase 8A**: create_project_config, detect_fastq_files, validate_input_data, run_pipeline
3. ✅ **Phase 8B**: monitor_pipeline, create_sample_sheet, estimate_resources
4. ✅ **Phase 9**: validate_library_type, setup_and_validate (RSeQC 라이브러리 검증)
5. ✅ **Phase 10**: run_de_pipeline, monitor_de_pipeline + DE 인지 도구 (조건 검증, 결과 요약)
6. **Phase 11**: ATAC-seq / WGS 지원 확장 (read_distribution 기반 자동 감지)

---

**Timeline Estimate**:
- Phase 6-7: 1-2 weeks (현재 agent 안정화)
- Phase 8A: 2-3 weeks (도구 구현 + 테스트)
- Phase 8B: 1-2 weeks (통합 + 검증)
- Phase 9-10: 4-6 weeks (다른 파이프라인 확장)

**Total**: ~3 months to full Option C implementation
