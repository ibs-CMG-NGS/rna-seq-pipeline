# Pipeline Standardization Guide

**목적**: 모든 NGS 파이프라인에 일관된 데이터 구조, 메타데이터 관리, 쿼리 인터페이스를 적용하기 위한 표준화 가이드

**적용 대상**: RNA-seq, ChIP-seq, ATAC-seq, WGS, WES 등 모든 NGS 분석 파이프라인

**작성일**: 2026-02-25  
**버전**: 1.0  
**참조 구현**: rna-seq-pipeline

---

## 📚 목차

1. [개요](#개요)
2. [재사용 가능한 코드 및 도구](#재사용-가능한-코드-및-도구)
3. [Phase 1: 데이터 구조 표준화](#phase-1-데이터-구조-표준화)
4. [Phase 2: 메타데이터 통합](#phase-2-메타데이터-통합)
5. [Phase 3: 프로젝트 레벨 요약](#phase-3-프로젝트-레벨-요약)
6. [Phase 4: 쿼리 인터페이스](#phase-4-쿼리-인터페이스)
7. [Phase 5: Pipeline 연결 자동화](#phase-5-pipeline-연결-자동화)
8. [Phase 6: LLM 통합 (선택)](#phase-6-llm-통합-선택)
9. [다른 파이프라인 적용 가이드](#다른-파이프라인-적용-가이드)
10. [체크리스트](#체크리스트)

---

## 개요

### 표준화의 필요성

**문제점**:
- 파이프라인마다 다른 출력 구조
- QC 결과를 수동으로 확인해야 함
- 프로젝트 전체 상태 파악이 어려움
- Downstream 분석으로 전환 시 수동 작업 필요
- 메타데이터가 흩어져 있거나 누락됨

**해결책**:
- 표준 디렉토리 구조
- 자동화된 메타데이터 수집
- 프로젝트 레벨 요약 생성
- 명령어/자연어 쿼리 인터페이스
- Pipeline 간 자동 브리지

### 표준화 이점

1. **재현성**: 일관된 구조로 분석 재현이 쉬움
2. **자동화**: QC 상태 자동 수집 및 요약
3. **효율성**: 수동 확인 작업 최소화
4. **통합**: Pipeline 간 seamless 연결
5. **확장성**: 새 프로젝트/샘플 추가가 용이
6. **협업**: 팀원 간 동일한 구조로 작업

### 6단계 표준화 프로세스

```
Phase 1: 데이터 구조화
    ↓
Phase 2: 메타데이터 통합
    ↓
Phase 3: 프로젝트 레벨 요약
    ↓
Phase 4: 쿼리 인터페이스
    ↓
Phase 5: Pipeline 연결
    ↓
Phase 6: LLM 통합 (선택)
```

---

## 재사용 가능한 코드 및 도구

**중요**: 다른 파이프라인에 표준화를 적용할 때는 rna-seq-pipeline의 기존 코드를 **복사하여 커스터마이징**하세요. 새로 작성하지 마세요!

### 파이프라인 독립적 도구 (그대로 사용)

이 도구들은 **모든 파이프라인에서 수정 없이 사용 가능**합니다:

| 도구 | 파일 경로 | 용도 | 사용법 |
|------|-----------|------|--------|
| **샘플 시트 변환기** | `scripts/standardization/convert_sample_sheet.py` | Master CSV → 각 파이프라인 형식 변환 | `python convert_sample_sheet.py master.csv -o config/` |
| **프로젝트 요약 생성기** | `scripts/standardization/generate_project_summary.py` | Manifest 통합 → Project summary | `python generate_project_summary.py --project-dir /output --project-id my-project` |
| **쿼리 에이전트** | `scripts/standardization/agent_query.py` | 명령줄 쿼리 인터페이스 | `python agent_query.py --project-summary summary.json --query status` |
| **LLM 에이전트** | `scripts/standardization/llm_agent.py` | 자연어 쿼리 인터페이스 | `python llm_agent.py --project-summary summary.json --interactive` |
| **Manifest 업데이트** | `scripts/standardization/update_manifests_metadata.py` | 기존 manifest에 metadata 추가 | `python update_manifests_metadata.py --sample-sheet samples.tsv --project-dir /output` |

### 파이프라인 특화 도구 (커스터마이징 필요)

이 도구들은 **파이프라인별 QC 메트릭에 맞춰 수정**이 필요합니다:

| 도구 | RNA-seq 파일 | 커스터마이징 필요 부분 | 적용 방법 |
|------|--------------|----------------------|-----------|
| **Manifest 생성기** | `scripts/standardization/generate_manifest.py` | `extract_qc_metrics()` 함수 | 파이프라인별 로그 파일 파싱 로직 수정 |
| **Manifest 생성기** | `scripts/standardization/generate_manifest.py` | `evaluate_qc_status()` 함수 | 파이프라인별 QC 기준 (threshold) 수정 |
| **Pipeline 브리지** | `scripts/bridge_to_de_pipeline.py` | 전체 클래스 | Downstream 파이프라인에 맞춰 입력 형식 조정 |

### 재사용 코드 구조

```
rna-seq-pipeline/
├── scripts/standardization/
│   ├── convert_sample_sheet.py        ← 그대로 복사
│   ├── generate_manifest.py           ← 함수만 수정
│   ├── update_manifests_metadata.py   ← 그대로 복사
│   ├── generate_project_summary.py    ← 그대로 복사
│   ├── agent_query.py                 ← 그대로 복사
│   └── llm_agent.py                   ← 도구 함수만 추가
└── scripts/
    └── bridge_to_de_pipeline.py       ← 전체 커스터마이징

# 새 파이프라인 적용 시:
chip-seq-pipeline/
├── scripts/standardization/
│   ├── convert_sample_sheet.py        ← [복사] rna-seq에서
│   ├── generate_manifest.py           ← [복사 후 수정] extract_qc_metrics_chipseq()
│   ├── update_manifests_metadata.py   ← [복사] rna-seq에서
│   ├── generate_project_summary.py    ← [복사] rna-seq에서
│   ├── agent_query.py                 ← [복사] rna-seq에서
│   └── llm_agent.py                   ← [복사] rna-seq에서
└── scripts/
    └── bridge_to_downstream.py        ← [새로 작성 또는 템플릿 수정]
```

### 단계별 재사용 가이드

#### Phase 1: 디렉토리 구조
```bash
# 1. RNA-seq 구조를 템플릿으로 복사
cp -r rna-seq-pipeline/scripts/standardization new-pipeline/scripts/

# 2. 수정 없이 바로 사용 가능
```

#### Phase 2: Manifest 생성
```bash
# 1. generate_manifest.py 복사
cp rna-seq-pipeline/scripts/standardization/generate_manifest.py \
   chip-seq-pipeline/scripts/standardization/

# 2. 아래 함수들만 수정:
#    - extract_qc_metrics() → extract_qc_metrics_chipseq()
#    - evaluate_qc_status() → evaluate_qc_status_chipseq()
```

**수정 예시 (ChIP-seq)**:
```python
# RNA-seq 함수 (원본)
def extract_qc_metrics(sample_id: str, output_dir: Path) -> dict:
    """RNA-seq 특화 메트릭"""
    metrics = {}
    
    # STAR 로그 파싱
    star_log = output_dir / "metadata" / "logs" / f"{sample_id}_Log.final.out"
    # ...
    
    return metrics

# ChIP-seq 함수 (수정)
def extract_qc_metrics_chipseq(sample_id: str, output_dir: Path) -> dict:
    """ChIP-seq 특화 메트릭 - RNA-seq 코드 기반"""
    metrics = {}
    
    # FRiP score 파싱 (ChIP-seq 특화)
    frip_file = output_dir / "metadata" / f"{sample_id}_frip.txt"
    if frip_file.exists():
        with open(frip_file) as f:
            metrics["frip"] = float(f.read().strip())
    
    # NSC/RSC 파싱 (ChIP-seq 특화)
    cc_file = output_dir / "metadata" / f"{sample_id}_cc.txt"
    # ...
    
    return metrics

# main() 함수에서 함수 이름만 변경
def main():
    # ...
    manifest["qc_metrics"] = extract_qc_metrics_chipseq(sample_id, output_dir)  # 변경
    # ...
```

#### Phase 3-6: 그대로 사용
```bash
# generate_project_summary.py - 수정 없이 사용
python scripts/standardization/generate_project_summary.py \
    --project-dir /data/chip-seq-output \
    --project-id my-chipseq

# agent_query.py - 수정 없이 사용
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query status

# llm_agent.py - 수정 없이 사용
python scripts/standardization/llm_agent.py \
    --project-summary project_summary.json \
    --interactive
```

### 코드 재사용 체크리스트

**그대로 복사 (0% 수정)**:
- [ ] `convert_sample_sheet.py`
- [ ] `update_manifests_metadata.py`
- [ ] `generate_project_summary.py`
- [ ] `agent_query.py`
- [ ] `llm_agent.py`

**함수 수정 (~20% 수정)**:
- [ ] `generate_manifest.py` → `extract_qc_metrics_XXX()`
- [ ] `generate_manifest.py` → `evaluate_qc_status_XXX()`

**새로 작성 또는 대폭 수정 (~80% 수정)**:
- [ ] `bridge_to_downstream.py` (Downstream이 다르면 새로 작성)

**총 예상 재사용률**: ~85% (7개 중 6개 파일 재사용)

---

## Phase 1: 데이터 구조 표준화

### 목표

모든 파이프라인 출력을 3-tier 구조로 통일:
1. **final_outputs/** - 최종 결과 파일
2. **intermediate/** - 중간 처리 파일
3. **metadata/** - 메타데이터 및 QC 정보

### 표준 디렉토리 구조

```
{project_output_dir}/
├── final_outputs/          # 최종 결과 (downstream 분석용)
│   ├── aligned/            # 정렬 파일 (BAM/CRAM)
│   ├── counts/             # 정량화 결과
│   ├── peaks/              # Peak calling 결과 (ChIP/ATAC)
│   ├── variants/           # Variant calling 결과 (WGS/WES)
│   └── qc_report.html      # 통합 QC 리포트
├── intermediate/           # 중간 파일 (재분석용)
│   ├── trimmed/            # 트리밍된 FASTQ
│   ├── preprocessing/      # 전처리 결과
│   └── temp/               # 임시 파일
└── metadata/               # 메타데이터
    ├── manifest.json       # 샘플별 상세 정보
    ├── qc_metrics.json     # QC 메트릭
    └── processing_log.json # 처리 이력
```

### 구현 방법

#### 1. Snakemake 규칙 수정

```python
# Snakefile 예시
OUTPUT_DIR = config["output_dir"]

rule all:
    input:
        # Final outputs
        f"{OUTPUT_DIR}/final_outputs/qc_report.html",
        f"{OUTPUT_DIR}/metadata/manifest.json"

rule align:
    output:
        bam = f"{OUTPUT_DIR}/final_outputs/aligned/{{sample}}.bam",
        log = f"{OUTPUT_DIR}/metadata/logs/{{sample}}_align.log"
    # ...

rule qc:
    output:
        metrics = f"{OUTPUT_DIR}/metadata/qc_metrics.json"
    # ...
```

#### 2. 디렉토리 자동 생성 스크립트

```python
# scripts/setup_output_structure.py
from pathlib import Path

def setup_output_dirs(base_dir: Path):
    """표준 출력 디렉토리 생성"""
    
    dirs = [
        base_dir / "final_outputs" / "aligned",
        base_dir / "final_outputs" / "counts",
        base_dir / "intermediate" / "trimmed",
        base_dir / "metadata" / "logs"
    ]
    
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        print(f"Created: {d}")

if __name__ == "__main__":
    import sys
    base_dir = Path(sys.argv[1])
    setup_output_dirs(base_dir)
```

사용법:
```bash
python scripts/setup_output_structure.py /path/to/output
```

#### 3. 파일 분류 기준

| 파일 종류 | 저장 위치 | 예시 |
|-----------|-----------|------|
| 최종 분석 결과 | `final_outputs/` | BAM, counts matrix, peaks |
| Downstream 입력 | `final_outputs/` | counts.txt, variants.vcf |
| 중간 처리 파일 | `intermediate/` | Trimmed FASTQ, SAM files |
| 로그 파일 | `metadata/logs/` | STAR logs, cutadapt logs |
| QC 결과 | `metadata/` | FastQC reports, metrics |
| 메타데이터 | `metadata/` | manifest.json, sample info |

### Phase 1 검증

```bash
# 구조가 올바른지 확인
tree -L 3 /path/to/output

# 예상 출력:
# /path/to/output/
# ├── final_outputs
# │   ├── aligned
# │   └── counts
# ├── intermediate
# │   └── trimmed
# └── metadata
#     └── logs
```

---

## Phase 2: 메타데이터 통합

### 목표

샘플 정보, QC 메트릭, 파일 메타데이터를 하나의 JSON manifest로 통합

### Manifest 구조 (manifest.json)

```json
{
  "sample_id": "Sample001",
  "processing_date": "2026-02-25T10:30:00",
  "pipeline_version": "1.0.0",
  
  "sample_metadata": {
    "condition": "wildtype",
    "replicate": 1,
    "tissue": "brain",
    "age": "P30",
    "sex": "male",
    "genotype": "WT",
    "treatment": "control",
    "batch": "batch1"
  },
  
  "qc_metrics": {
    "total_reads": 50000000,
    "mapped_reads": 46000000,
    "mapping_rate": 92.0,
    "duplicate_rate": 15.2,
    "properly_paired_rate": 95.5,
    
    "fastqc_status": "PASS",
    "fastqc_warnings": [],
    "fastqc_failures": [],
    
    "alignment_rate": 92.0,
    "uniquely_mapped_rate": 85.3,
    "multi_mapped_rate": 6.7,
    
    "assignment_rate": 70.5,
    "assigned_reads": 32000000,
    "unassigned_ambiguity": 10.2,
    "unassigned_no_features": 19.3
  },
  
  "final_outputs": {
    "aligned": {
      "path": "final_outputs/aligned/Sample001.bam",
      "md5": "abc123...",
      "size_bytes": 5000000000,
      "format": "BAM"
    },
    "counts": {
      "path": "final_outputs/counts/Sample001.counts.txt",
      "md5": "def456...",
      "size_bytes": 1000000,
      "format": "TXT"
    }
  },
  
  "qc_status": "PASS",
  "issues": []
}
```

### 구현 방법

#### 1. Manifest 생성 스크립트

**🔄 재사용**: `rna-seq-pipeline/scripts/standardization/generate_manifest.py`를 복사하고 아래 함수만 수정하세요.

**원본 파일 위치**: `rna-seq-pipeline/scripts/standardization/generate_manifest.py`

**수정이 필요한 부분**:

1. **`extract_qc_metrics()` 함수** - 파이프라인별 로그 파일 파싱
2. **`evaluate_qc_status()` 함수** - 파이프라인별 QC 기준

**RNA-seq 버전 (참조용)**:
```python
# rna-seq-pipeline/scripts/standardization/generate_manifest.py
def extract_qc_metrics(sample_id: str, output_dir: Path) -> dict:
    """
    RNA-seq QC 메트릭 추출
    
    ✅ 이 함수를 파이프라인에 맞게 수정하세요!
    """
    metrics = {}
    
    # STAR 로그 파싱 (RNA-seq 특화)
    star_log = output_dir / "metadata" / "logs" / f"{sample_id}_Log.final.out"
    if star_log.exists():
        with open(star_log) as f:
            for line in f:
                if "Uniquely mapped reads %" in line:
                    metrics["uniquely_mapped_rate"] = float(line.split("|")[1].strip().rstrip('%'))
                elif "Number of input reads" in line:
                    metrics["total_reads"] = int(line.split("|")[1].strip())
    
    # featureCounts 결과 파싱 (RNA-seq 특화)
    counts_summary = output_dir / "final_outputs" / "counts" / f"{sample_id}.counts.txt.summary"
    if counts_summary.exists():
        with open(counts_summary) as f:
            lines = f.readlines()
            assigned = int(lines[1].split()[1])
            total = sum(int(line.split()[1]) for line in lines[1:])
            metrics["assignment_rate"] = (assigned / total * 100) if total > 0 else 0
    
    return metrics

def evaluate_qc_status(metrics: dict) -> tuple[str, list]:
    """
    QC 메트릭 기반 PASS/WARN/FAIL 판정
    
    ✅ 이 함수의 threshold를 파이프라인에 맞게 조정하세요!
    """
    status = "PASS"
    issues = []
    
    # RNA-seq 기준
    if metrics.get("mapping_rate", 0) < 70:
        status = "FAIL"
        issues.append("Low mapping rate")
    
    if metrics.get("duplicate_rate", 0) > 50:
        status = "WARN"
        issues.append("High duplicate rate")
    
    if metrics.get("assignment_rate", 0) < 50:
        status = "WARN"
        issues.append("Low gene assignment rate")
    
    return status, issues
```

**ChIP-seq 버전 (수정 예시)**:
```python
# chip-seq-pipeline/scripts/standardization/generate_manifest.py
# RNA-seq 파일을 복사한 후 아래 함수들만 수정

def extract_qc_metrics_chipseq(sample_id: str, output_dir: Path) -> dict:
    """
    ChIP-seq QC 메트릭 추출
    
    ✏️ RNA-seq 함수를 복사해서 ChIP-seq 로그 파일로 수정
    """
    metrics = {}
    
    # BWA/Bowtie2 로그 파싱 (ChIP-seq)
    align_log = output_dir / "metadata" / "logs" / f"{sample_id}_align.log"
    if align_log.exists():
        with open(align_log) as f:
            for line in f:
                if "mapped (" in line:
                    pct = line.split("mapped (")[1].split("%")[0]
                    metrics["mapping_rate"] = float(pct)
    
    # FRiP score (ChIP-seq 특화)
    frip_file = output_dir / "metadata" / f"{sample_id}_frip.txt"
    if frip_file.exists():
        with open(frip_file) as f:
            metrics["frip"] = float(f.read().strip())
    
    # NSC/RSC from phantompeakqualtools (ChIP-seq 특화)
    cc_file = output_dir / "metadata" / f"{sample_id}_cc.txt"
    if cc_file.exists():
        with open(cc_file) as f:
            parts = f.read().strip().split('\t')
            metrics["nsc"] = float(parts[8])
            metrics["rsc"] = float(parts[9])
    
    # Peak count (ChIP-seq 특화)
    peaks_file = output_dir / "final_outputs" / "peaks" / f"{sample_id}_peaks.narrowPeak"
    if peaks_file.exists():
        with open(peaks_file) as f:
            metrics["peak_count"] = sum(1 for line in f)
    
    return metrics

def evaluate_qc_status_chipseq(metrics: dict) -> tuple[str, list]:
    """
    ChIP-seq QC 기준
    
    ✏️ RNA-seq 함수를 복사해서 ChIP-seq threshold로 수정
    """
    status = "PASS"
    issues = []
    
    # ChIP-seq 기준
    if metrics.get("frip", 0) < 0.01:  # FRiP < 1%
        status = "FAIL"
        issues.append("Very low FRiP score")
    
    if metrics.get("peak_count", 0) < 1000:
        status = "WARN"
        issues.append("Low peak count")
    
    if metrics.get("nsc", 0) < 1.05:
        status = "WARN"
        issues.append("Low NSC (poor enrichment)")
    
    if metrics.get("rsc", 0) < 0.8:
        status = "WARN"
        issues.append("Low RSC")
    
    return status, issues

# main() 함수에서 함수 이름만 변경
def main():
    # ...
    manifest["qc_metrics"] = extract_qc_metrics_chipseq(sample_id, output_dir)
    manifest["qc_status"], manifest["issues"] = evaluate_qc_status_chipseq(
        manifest["qc_metrics"]
    )
    # ...
```

**적용 방법**:
```bash
# 1. RNA-seq 파일 복사
cp rna-seq-pipeline/scripts/standardization/generate_manifest.py \
   chip-seq-pipeline/scripts/standardization/

# 2. 함수 이름 변경 및 로직 수정
# - extract_qc_metrics() → extract_qc_metrics_chipseq()
# - evaluate_qc_status() → evaluate_qc_status_chipseq()

# 3. main()에서 함수 호출 부분 수정
```

**나머지 코드는 그대로 사용**:
- `calculate_md5()` - 수정 없음
- `load_sample_metadata()` - 수정 없음 (샘플 시트 형식이 같으므로)
- `generate_manifest()` - 함수 호출만 변경
- `main()` - argparse 및 저장 로직 그대로

### 구현 방법 (기존 코드)

```python
# scripts/standardization/generate_manifest.py
import json
from pathlib import Path
from datetime import datetime
import hashlib

### 구현 방법 (기존 코드)

**⚠️ 중요**: 아래 긴 코드를 새로 작성하지 마세요! 대신:

1. **RNA-seq 파일 복사**: `rna-seq-pipeline/scripts/standardization/generate_manifest.py`
2. **함수 2개만 수정**: `extract_qc_metrics()`, `evaluate_qc_status()`
3. **나머지 그대로 사용**: 95%의 코드는 수정 불필요

**전체 코드 보기**: `rna-seq-pipeline/scripts/standardization/generate_manifest.py` (약 200 lines)

#### 2. Snakemake 통합
    """
    status = "PASS"
    issues = []
    
    # RNA-seq 기준 예시
    if metrics.get("mapping_rate", 0) < 70:
        status = "FAIL"
        issues.append("Low mapping rate")
    
    if metrics.get("duplicate_rate", 0) > 50:
        status = "WARN"
        issues.append("High duplicate rate")
    
    if metrics.get("assignment_rate", 0) < 50:
        status = "WARN"
        issues.append("Low gene assignment rate")
    
    return status, issues

def generate_manifest(sample_id: str, 
                     output_dir: Path,
                     sample_sheet: Path,
                     pipeline_version: str = "1.0.0") -> dict:
    """샘플별 manifest 생성"""
    
    manifest = {
        "sample_id": sample_id,
        "processing_date": datetime.now().isoformat(),
        "pipeline_version": pipeline_version,
        "sample_metadata": load_sample_metadata(sample_id, sample_sheet),
        "qc_metrics": extract_qc_metrics(sample_id, output_dir),
        "final_outputs": {}
    }
    
    # 최종 출력 파일 메타데이터
    bam_file = output_dir / "final_outputs" / "aligned" / f"{sample_id}.bam"
    if bam_file.exists():
        manifest["final_outputs"]["aligned"] = {
            "path": str(bam_file.relative_to(output_dir)),
            "md5": calculate_md5(bam_file),
            "size_bytes": bam_file.stat().st_size,
            "format": "BAM"
        }
    
    # QC 상태 평가
    manifest["qc_status"], manifest["issues"] = evaluate_qc_status(
        manifest["qc_metrics"]
    )
    
    return manifest

def main():
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--pipeline-version", default="1.0.0")
    
    args = parser.parse_args()
    
    manifest = generate_manifest(
        args.sample_id,
        args.output_dir,
        args.sample_sheet,
        args.pipeline_version
    )
    
    # 저장
    manifest_dir = args.output_dir / "metadata"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    
    manifest_file = manifest_dir / f"{args.sample_id}_manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"✅ Manifest saved: {manifest_file}")
    print(f"   QC Status: {manifest['qc_status']}")
    if manifest['issues']:
        print(f"   Issues: {', '.join(manifest['issues'])}")

if __name__ == "__main__":
    main()
```

#### 2. Snakemake 통합

```python
# Snakefile에 추가
rule generate_manifest:
    input:
        bam = f"{OUTPUT_DIR}/final_outputs/aligned/{{sample}}.bam",
        counts = f"{OUTPUT_DIR}/final_outputs/counts/{{sample}}.counts.txt",
        star_log = f"{OUTPUT_DIR}/metadata/logs/{{sample}}_Log.final.out"
    output:
        manifest = f"{OUTPUT_DIR}/metadata/{{sample}}_manifest.json"
    params:
        sample_sheet = config["sample_sheet"],
        output_dir = OUTPUT_DIR,
        pipeline_version = "1.0.0"
    shell:
        """
        python scripts/standardization/generate_manifest.py \
            --sample-id {wildcards.sample} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
            --pipeline-version {params.pipeline_version}
        """
```

### Phase 2 검증

```bash
# Manifest 생성 확인
ls metadata/*_manifest.json

# Manifest 내용 확인
cat metadata/Sample001_manifest.json | jq '.qc_status, .issues'

# 모든 샘플의 QC 상태 확인
for f in metadata/*_manifest.json; do
    echo "$(basename $f): $(jq -r '.qc_status' $f)"
done
```

---

## Phase 3: 프로젝트 레벨 요약

### 목표

모든 샘플의 manifest를 통합하여 프로젝트 전체 상태를 한눈에 파악

### Project Summary 구조 (project_summary.json)

```json
{
  "project_id": "mouse-chd8",
  "generated_at": "2026-02-25T15:30:00",
  "pipeline_version": "1.0.0",
  
  "summary": {
    "total_samples": 38,
    "completed_samples": 38,
    "failed_samples": 0,
    "processing_samples": 0
  },
  
  "qc_summary": {
    "total": 38,
    "passed": 38,
    "warned": 0,
    "failed": 0,
    "pass_rate": 100.0
  },
  
  "condition_groups": {
    "wildtype": {
      "count": 20,
      "passed": 20,
      "failed": 0,
      "samples": ["WT_01", "WT_02", "..."]
    },
    "heterozygous": {
      "count": 18,
      "passed": 18,
      "failed": 0,
      "samples": ["HET_01", "HET_02", "..."]
    }
  },
  
  "aggregate_stats": {
    "mapping_rate": {
      "mean": 92.4,
      "std": 1.2,
      "min": 89.5,
      "max": 94.8
    },
    "assignment_rate": {
      "mean": 70.2,
      "std": 2.5,
      "min": 65.3,
      "max": 75.6
    }
  },
  
  "failed_samples_details": [],
  
  "warnings": [
    {
      "sample_id": "Sample003",
      "issue": "High duplicate rate (45%)",
      "severity": "WARN"
    }
  ]
}
```

### 구현 방법

**🔄 재사용**: `rna-seq-pipeline/scripts/standardization/generate_project_summary.py`를 **그대로 복사**하세요!

**원본 파일**: `rna-seq-pipeline/scripts/standardization/generate_project_summary.py` (408 lines)

**수정 필요**: ❌ **없음** - 모든 파이프라인에서 그대로 사용 가능

**이 스크립트가 하는 일**:
1. `metadata/*_manifest.json` 파일들을 모두 수집
2. QC 상태 집계 (total/passed/failed/pass_rate)
3. 조건별 샘플 그룹화
4. 통계 지표 계산 (mean/std/min/max)
5. 문제 샘플 및 경고 식별

**핵심 함수**:
- `collect_manifests()` - Manifest 파일 수집
- `summarize_qc_status()` - QC 상태 요약
- `group_by_condition()` - 조건별 그룹화
- `calculate_aggregate_stats()` - 통계 계산
- `identify_issues()` - 문제 샘플 식별

**전체 코드 보기**: `rna-seq-pipeline/scripts/standardization/generate_project_summary.py`

#### 사용법
    warned = sum(1 for m in manifests if m["qc_status"] == "WARN")
    failed = sum(1 for m in manifests if m["qc_status"] == "FAIL")
    
    return {
        "total": total,
        "passed": passed,
        "warned": warned,
        "failed": failed,
        "pass_rate": round(passed / total * 100, 2) if total > 0 else 0
    }

def group_by_condition(manifests: List[Dict]) -> Dict:
    """조건별 샘플 그룹화"""
    groups = {}
    
    for manifest in manifests:
        condition = manifest["sample_metadata"].get("condition", "Unknown")
        
        if condition not in groups:
            groups[condition] = {
                "count": 0,
                "passed": 0,
                "failed": 0,
                "samples": []
            }
        
        groups[condition]["count"] += 1
        groups[condition]["samples"].append(manifest["sample_id"])
        
        if manifest["qc_status"] == "PASS":
            groups[condition]["passed"] += 1
        elif manifest["qc_status"] == "FAIL":
            groups[condition]["failed"] += 1
    
    return groups

def calculate_aggregate_stats(manifests: List[Dict]) -> Dict:
    """통계 지표 계산"""
    stats = {}
    
    # 수집할 메트릭 목록 (파이프라인별로 커스터마이징)
    metrics_to_aggregate = [
        "mapping_rate",
        "assignment_rate",
        "duplicate_rate",
        "properly_paired_rate"
    ]
    
    for metric in metrics_to_aggregate:
        values = [
            m["qc_metrics"].get(metric)
            for m in manifests
            if metric in m["qc_metrics"] and m["qc_metrics"][metric] is not None
        ]
        
        if values:
            stats[metric] = {
                "mean": round(np.mean(values), 2),
                "std": round(np.std(values), 2),
                "min": round(np.min(values), 2),
                "max": round(np.max(values), 2),
                "median": round(np.median(values), 2)
            }
    
    return stats

def identify_issues(manifests: List[Dict]) -> tuple[List[Dict], List[Dict]]:
    """문제 샘플 및 경고 식별"""
    failed = []
    warnings = []
    
    for manifest in manifests:
        if manifest["qc_status"] == "FAIL":
            failed.append({
                "sample_id": manifest["sample_id"],
                "issues": manifest["issues"],
                "qc_metrics": manifest["qc_metrics"]
            })
        
        elif manifest["qc_status"] == "WARN":
            for issue in manifest["issues"]:
                warnings.append({
                    "sample_id": manifest["sample_id"],
                    "issue": issue,
                    "severity": "WARN"
                })
    
    return failed, warnings

def generate_project_summary(project_dir: Path, project_id: str) -> Dict:
    """프로젝트 요약 생성"""
    
    manifests = collect_manifests(project_dir)
    
    if not manifests:
        raise ValueError(f"No manifest files found in {project_dir}/metadata/")
    
    # 파이프라인 버전 (첫 번째 manifest에서 가져옴)
    pipeline_version = manifests[0].get("pipeline_version", "unknown")
    
    failed_samples, warnings = identify_issues(manifests)
    
    summary = {
        "project_id": project_id,
        "generated_at": datetime.now().isoformat(),
        "pipeline_version": pipeline_version,
        
        "summary": {
            "total_samples": len(manifests),
            "completed_samples": len(manifests),
            "failed_samples": len(failed_samples),
            "processing_samples": 0
        },
        
        "qc_summary": summarize_qc_status(manifests),
        "condition_groups": group_by_condition(manifests),
        "aggregate_stats": calculate_aggregate_stats(manifests),
        "failed_samples_details": failed_samples,
        "warnings": warnings
    }
    
    return summary

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate project-level summary from sample manifests"
    )
    parser.add_argument("--project-dir", type=Path, required=True,
                       help="Project output directory containing metadata/")
    parser.add_argument("--project-id", required=True,
                       help="Project identifier")
    parser.add_argument("--output", type=Path,
                       help="Output file path (default: project_dir/project_summary.json)")
    
    args = parser.parse_args()
    
    # 요약 생성
    summary = generate_project_summary(args.project_dir, args.project_id)
    
    # 저장
    output_file = args.output or (args.project_dir / "project_summary.json")
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # 결과 출력
    print(f"\n✅ Project summary generated: {output_file}")
    print(f"\n📊 Project: {summary['project_id']}")
    print(f"   Total samples: {summary['summary']['total_samples']}")
    print(f"   QC pass rate: {summary['qc_summary']['pass_rate']}%")
    print(f"   Failed samples: {summary['summary']['failed_samples']}")
    
    print(f"\n📁 Conditions:")
    for condition, info in summary['condition_groups'].items():
        print(f"   - {condition}: {info['count']} samples ({info['passed']} passed)")
    
    if summary['failed_samples_details']:
        print(f"\n⚠️  Failed samples:")
        for failed in summary['failed_samples_details']:
            print(f"   - {failed['sample_id']}: {', '.join(failed['issues'])}")
    
    if summary['warnings']:
        print(f"\n⚠️  Warnings: {len(summary['warnings'])}")

if __name__ == "__main__":
    main()
```

#### 2. 사용법

```bash
# 프로젝트 요약 생성
python scripts/standardization/generate_project_summary.py \
    --project-dir /data/output/mouse-chd8 \
    --project-id mouse-chd8

# 출력 예시:
# ✅ Project summary generated: /data/output/mouse-chd8/project_summary.json
#
# 📊 Project: mouse-chd8
#    Total samples: 38
#    QC pass rate: 100.0%
#    Failed samples: 0
#
# 📁 Conditions:
#    - wildtype: 20 samples (20 passed)
#    - heterozygous: 18 samples (18 passed)
```

### Phase 3 검증

```bash
# Project summary 확인
cat project_summary.json | jq '.qc_summary'

# 조건별 요약
cat project_summary.json | jq '.condition_groups'

# 실패 샘플 확인
cat project_summary.json | jq '.failed_samples_details'

# 통계 지표
cat project_summary.json | jq '.aggregate_stats.mapping_rate'
```

---

## Phase 4: 쿼리 인터페이스

### 목표

project_summary.json을 활용한 명령줄 쿼리 도구 제공

### 지원 쿼리 타입

1. **status** - 전체 프로젝트 상태
2. **failed** - 실패한 샘플 목록
3. **condition** - 특정 조건의 샘플 목록
4. **compare** - 조건 간 비교
5. **sample** - 특정 샘플 상세 정보
6. **stats** - 통계 지표
7. **issues** - 모든 문제/경고
8. **conditions** - 조건 목록

### 구현 방법

**🔄 재사용**: `rna-seq-pipeline/scripts/standardization/agent_query.py`를 **그대로 복사**하세요!

**원본 파일**: `rna-seq-pipeline/scripts/standardization/agent_query.py` (408 lines)

**수정 필요**: ❌ **없음** - 모든 파이프라인에서 그대로 사용 가능

**이 스크립트가 하는 일**:
1. `project_summary.json`을 로딩
2. 8가지 쿼리 타입 지원
3. JSON/text 출력 포맷 제공
4. CLI 인터페이스

**핵심 클래스 및 메서드**:
```python
class PipelineQueryAgent:
    def get_overall_status()      # 전체 상태
    def get_failed_samples()      # 실패 샘플
    def get_samples_by_condition() # 조건별 샘플
    def compare_conditions()      # 조건 비교
    def get_sample_details()      # 샘플 상세
    def list_conditions()         # 조건 목록
    def get_aggregate_stats()     # 통계 지표
    def get_all_issues()          # 모든 문제/경고
```

**전체 코드 보기**: `rna-seq-pipeline/scripts/standardization/agent_query.py`

### 사용 예시
        }
    
    def get_failed_samples(self) -> list:
        """실패 샘플 목록"""
        return self.summary['failed_samples_details']
    
    def get_samples_by_condition(self, condition: str) -> dict:
        """특정 조건의 샘플"""
        if condition not in self.summary['condition_groups']:
            available = list(self.summary['condition_groups'].keys())
            raise ValueError(f"Condition '{condition}' not found. Available: {available}")
        
        return self.summary['condition_groups'][condition]
    
    def compare_conditions(self, condition1: str, condition2: str) -> dict:
        """두 조건 비교"""
        c1 = self.get_samples_by_condition(condition1)
        c2 = self.get_samples_by_condition(condition2)
        
        # 각 조건의 메트릭 계산 (manifest에서)
        # ... (구현 생략, 필요시 manifest 재로딩)
        
        return {
            "condition1": {
                "name": condition1,
                "count": c1['count'],
                "passed": c1['passed']
            },
            "condition2": {
                "name": condition2,
                "count": c2['count'],
                "passed": c2['passed']
            }
        }
    
    def get_sample_details(self, sample_id: str, project_dir: Path) -> dict:
        """특정 샘플 상세 정보 (manifest 로딩)"""
        manifest_path = project_dir / "metadata" / f"{sample_id}_manifest.json"
        
        if not manifest_path.exists():
            raise ValueError(f"Sample '{sample_id}' not found")
        
        with open(manifest_path) as f:
            return json.load(f)
    
    def list_conditions(self) -> list:
        """모든 조건 목록"""
        return [
            {
                "condition": cond,
                "count": info['count'],
                "passed": info['passed']
            }
            for cond, info in self.summary['condition_groups'].items()
        ]
    
    def get_aggregate_stats(self) -> dict:
        """통계 지표"""
        return self.summary['aggregate_stats']
    
    def get_all_issues(self) -> dict:
        """모든 문제 및 경고"""
        return {
            "failed": self.summary['failed_samples_details'],
            "warnings": self.summary['warnings']
        }
    
    def query(self, query_type: str, **kwargs) -> dict:
        """쿼리 실행"""
        
        handlers = {
            "status": self.get_overall_status,
            "failed": self.get_failed_samples,
            "condition": lambda: self.get_samples_by_condition(kwargs['condition']),
            "compare": lambda: self.compare_conditions(kwargs['condition1'], kwargs['condition2']),
            "sample": lambda: self.get_sample_details(kwargs['sample_id'], kwargs['project_dir']),
            "stats": self.get_aggregate_stats,
            "issues": self.get_all_issues,
            "conditions": self.list_conditions
        }
        
        if query_type not in handlers:
            raise ValueError(f"Unknown query type: {query_type}")
        
        return handlers[query_type]()

def format_output(result: dict, format_type: str = "json") -> str:
    """출력 포맷팅"""
    if format_type == "json":
        return json.dumps(result, indent=2)
    
    elif format_type == "text":
        # 텍스트 포맷 (사람이 읽기 쉽게)
        lines = []
        
        if isinstance(result, dict):
            for key, value in result.items():
                if isinstance(value, (dict, list)):
                    lines.append(f"{key}:")
                    lines.append(f"  {json.dumps(value, indent=2)}")
                else:
                    lines.append(f"{key}: {value}")
        elif isinstance(result, list):
            for item in result:
                lines.append(json.dumps(item, indent=2))
        
        return "\n".join(lines)
    
    else:
        raise ValueError(f"Unknown format: {format_type}")

def main():
    parser = argparse.ArgumentParser(
        description="Query RNA-seq pipeline status"
    )
    parser.add_argument("--project-summary", type=Path, required=True,
                       help="Path to project_summary.json")
    parser.add_argument("--query", required=True,
                       choices=["status", "failed", "condition", "compare", 
                               "sample", "stats", "issues", "conditions"],
                       help="Query type")
    parser.add_argument("--condition", help="Condition name (for condition/compare)")
    parser.add_argument("--condition1", help="First condition (for compare)")
    parser.add_argument("--condition2", help="Second condition (for compare)")
    parser.add_argument("--sample-id", help="Sample ID (for sample query)")
    parser.add_argument("--project-dir", type=Path, help="Project directory (for sample query)")
    parser.add_argument("--format", choices=["json", "text"], default="json",
                       help="Output format")
    
    args = parser.parse_args()
    
    # 에이전트 초기화
    agent = PipelineQueryAgent(args.project_summary)
    
    # 쿼리 실행
    try:
        kwargs = {}
        if args.condition:
            kwargs['condition'] = args.condition
        if args.condition1:
            kwargs['condition1'] = args.condition1
        if args.condition2:
            kwargs['condition2'] = args.condition2
        if args.sample_id:
            kwargs['sample_id'] = args.sample_id
        if args.project_dir:
            kwargs['project_dir'] = args.project_dir
        
        result = agent.query(args.query, **kwargs)
        print(format_output(result, args.format))
    
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    import sys
    main()
```

### 사용 예시

```bash
# 전체 상태
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query status

# 실패 샘플
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query failed

# 특정 조건
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query condition \
    --condition wildtype

# 조건 비교
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query compare \
    --condition1 wildtype \
    --condition2 knockout

# 샘플 상세 정보
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query sample \
    --sample-id Sample001 \
    --project-dir /data/output/project

# 통계 지표
python scripts/standardization/agent_query.py \
    --project-summary project_summary.json \
    --query stats \
    --format text
```

### Phase 4 검증

```bash
# 모든 쿼리 타입 테스트
for query in status failed conditions stats issues; do
    echo "Testing: $query"
    python scripts/standardization/agent_query.py \
        --project-summary project_summary.json \
        --query $query
done
```

---

## Phase 5: Pipeline 연결 자동화

### 목표

RNA-seq 완료 후 downstream 분석 (DE/GO) 자동 연결

### 구현 방법

**🔧 커스터마이징**: `rna-seq-pipeline/scripts/bridge_to_de_pipeline.py`를 **템플릿으로 사용**하세요.

**원본 파일**: `rna-seq-pipeline/scripts/bridge_to_de_pipeline.py` (315+ lines)

**수정 필요**: ⚠️ **Downstream 파이프라인에 따라 수정** 필요

**이 스크립트가 하는 일**:
1. RNA-seq QC 완료 여부 확인
2. Counts matrix를 DE pipeline으로 복사
3. 메타데이터 CSV 자동 생성
4. DE/GO config 파일 자동 생성
5. DE/GO pipeline 자동 실행 (optional)

**핵심 클래스 및 메서드**:
```python
class PipelineBridge:
    def check_rnaseq_completion()    # QC 상태 확인
    def prepare_de_input()           # Counts matrix 복사
    def generate_metadata_template() # Metadata CSV 생성
    def generate_de_config()         # Config 자동 생성
    def trigger_de_analysis()        # DE pipeline 실행
    def run()                        # 전체 프로세스
```

**수정이 필요한 부분**:
- `prepare_de_input()`: Downstream 파이프라인의 입력 형식에 맞춰 수정
- `generate_metadata_template()`: Downstream이 요구하는 메타데이터 형식으로 수정
- `generate_de_config()`: Downstream config 템플릿에 맞춰 수정

**전체 코드 보기**: `rna-seq-pipeline/scripts/bridge_to_de_pipeline.py`

### Bridge 스크립트 커스터마이징 예시

**RNA-seq → DE/GO (원본)**:
    
    def check_rnaseq_completion(self) -> bool:
        """RNA-seq 완료 여부 확인"""
        qc_pass_rate = self.summary['qc_summary']['pass_rate']
        failed_count = self.summary['summary']['failed_samples']
        
        print(f"📊 RNA-seq QC Summary:")
        print(f"   Pass rate: {qc_pass_rate}%")
        print(f"   Failed samples: {failed_count}")
        
        if failed_count > 0:
            print(f"\n⚠️  Warning: {failed_count} samples failed QC")
            response = input("Continue anyway? (yes/no): ")
            if response.lower() != 'yes':
                return False
        
        return True
    
    def prepare_de_input(self) -> Path:
        """Counts matrix 복사"""
        
        # Source: RNA-seq 출력
        source_counts = self.rnaseq_output / "final_outputs" / "counts" / "counts_matrix.txt"
        
        if not source_counts.exists():
            raise FileNotFoundError(f"Counts matrix not found: {source_counts}")
        
        # Destination: DE pipeline 입력
        dest_dir = self.de_pipeline_dir / "data" / "raw"
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        dest_counts = dest_dir / f"{self.project_id}_counts.txt"
        shutil.copy2(source_counts, dest_counts)
        
        print(f"✅ Copied counts matrix:")
        print(f"   From: {source_counts}")
        print(f"   To: {dest_counts}")
        
        return dest_counts
    
    def generate_metadata_template(self) -> Path:
        """메타데이터 CSV 생성"""
        
        # Sample sheet 로딩
        sample_sheet = self.rnaseq_output / "samples.tsv"  # 경로는 프로젝트에 따라 조정
        df = pd.read_csv(sample_sheet, sep='\t')
        
        # DE 분석용 메타데이터 형식으로 변환
        metadata = df[['sample_id', 'condition', 'replicate']].copy()
        metadata.columns = ['sample', 'condition', 'replicate']
        
        # 저장
        dest_dir = self.de_pipeline_dir / "data" / "metadata"
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        metadata_file = dest_dir / f"{self.project_id}_metadata.csv"
        metadata.to_csv(metadata_file, index=False)
        
        print(f"✅ Generated metadata:")
        print(f"   File: {metadata_file}")
        print(f"   Samples: {len(metadata)}")
        
        return metadata_file
    
    def generate_de_config(self) -> Path:
        """DE/GO 설정 파일 자동 생성"""
        
        # 템플릿 로딩
        template_path = self.de_pipeline_dir / "configs" / "template" / "config.yml"
        with open(template_path) as f:
            config = yaml.safe_load(f)
        
        # 프로젝트별 설정 업데이트
        config['project_id'] = self.project_id
        config['counts_file'] = f"data/raw/{self.project_id}_counts.txt"
        config['metadata_file'] = f"data/metadata/{self.project_id}_metadata.csv"
        
        # 조건 자동 감지
        conditions = list(self.summary['condition_groups'].keys())
        config['conditions'] = conditions
        
        # 대조군 자동 선택 (일반적으로 첫 번째 조건)
        if "wildtype" in conditions:
            config['control_condition'] = "wildtype"
        elif "control" in conditions:
            config['control_condition'] = "control"
        else:
            config['control_condition'] = conditions[0]
            print(f"⚠️  Auto-selected control: {conditions[0]}")
        
        # 비교 조합 자동 생성
        comparisons = []
        control = config['control_condition']
        for cond in conditions:
            if cond != control:
                comparisons.append({
                    "name": f"{cond}_vs_{control}",
                    "condition": cond,
                    "control": control
                })
        config['comparisons'] = comparisons
        
        # 저장
        config_dir = self.de_pipeline_dir / "configs"
        config_file = config_dir / f"config_{self.project_id}.yml"
        
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        
        print(f"✅ Generated DE config:")
        print(f"   File: {config_file}")
        print(f"   Comparisons: {len(comparisons)}")
        for comp in comparisons:
            print(f"      - {comp['name']}")
        
        return config_file
    
    def trigger_de_analysis(self, config_file: Path, dry_run: bool = False):
        """DE/GO 파이프라인 실행"""
        
        cmd = [
            "snakemake",
            "--configfile", str(config_file),
            "--cores", "8"
        ]
        
        if dry_run:
            cmd.append("--dry-run")
        
        print(f"\n{'🔍 Dry-run' if dry_run else '🚀 Starting'} DE/GO pipeline...")
        print(f"Command: {' '.join(cmd)}")
        
        if not dry_run:
            response = input("\nProceed with DE analysis? (yes/no): ")
            if response.lower() != 'yes':
                print("Cancelled.")
                return
        
        subprocess.run(cmd, cwd=self.de_pipeline_dir, check=True)
        
        print("\n✅ DE/GO pipeline completed!")
    
    def run(self, dry_run: bool = False, skip_trigger: bool = False):
        """전체 브리지 프로세스 실행"""
        
        print(f"\n{'='*60}")
        print(f"RNA-seq → DE/GO Pipeline Bridge")
        print(f"Project: {self.project_id}")
        print(f"{'='*60}\n")
        
        # Step 1: RNA-seq 완료 확인
        print("Step 1: Checking RNA-seq completion...")
        if not self.check_rnaseq_completion():
            print("❌ RNA-seq not ready. Exiting.")
            return
        
        # Step 2: Counts matrix 복사
        print("\nStep 2: Preparing DE input...")
        self.prepare_de_input()
        
        # Step 3: 메타데이터 생성
        print("\nStep 3: Generating metadata...")
        self.generate_metadata_template()
        
        # Step 4: Config 생성
        print("\nStep 4: Generating DE config...")
        config_file = self.generate_de_config()
        
        # Step 5: DE 분석 실행
        if not skip_trigger:
            print("\nStep 5: Triggering DE analysis...")
            self.trigger_de_analysis(config_file, dry_run=dry_run)
        else:
            print("\n⏭️  Skipping DE analysis trigger (--skip-trigger)")
            print(f"   Run manually: cd {self.de_pipeline_dir} && snakemake --configfile {config_file}")
        
        print(f"\n{'='*60}")
        print("✅ Bridge process completed!")
        print(f"{'='*60}\n")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Bridge RNA-seq pipeline to DE/GO analysis"
    )
    parser.add_argument("--rnaseq-output", type=Path, required=True,
                       help="RNA-seq output directory")
    parser.add_argument("--de-pipeline", type=Path, required=True,
                       help="DE/GO pipeline directory")
    parser.add_argument("--project-id", required=True,
                       help="Project identifier")
    parser.add_argument("--dry-run", action="store_true",
                       help="Dry run (don't actually run DE pipeline)")
    parser.add_argument("--skip-trigger", action="store_true",
                       help="Skip triggering DE pipeline (setup only)")
    
    args = parser.parse_args()
    
    bridge = PipelineBridge(
        rnaseq_output=args.rnaseq_output,
        de_pipeline_dir=args.de_pipeline,
        project_id=args.project_id
    )
    
    bridge.run(dry_run=args.dry_run, skip_trigger=args.skip_trigger)

if __name__ == "__main__":
    main()
```

### 사용 예시

```bash
# Setup만 (DE 실행 안 함)
python scripts/bridge_to_de_pipeline.py \
    --rnaseq-output /data/output/mouse-chd8 \
    --de-pipeline /home/user/de-go-pipeline \
    --project-id mouse-chd8 \
    --skip-trigger

# Dry-run으로 테스트
python scripts/bridge_to_de_pipeline.py \
    --rnaseq-output /data/output/mouse-chd8 \
    --de-pipeline /home/user/de-go-pipeline \
    --project-id mouse-chd8 \
    --dry-run

# 전체 실행
python scripts/bridge_to_de_pipeline.py \
    --rnaseq-output /data/output/mouse-chd8 \
    --de-pipeline /home/user/de-go-pipeline \
    --project-id mouse-chd8
```

### Phase 5 검증

```bash
# DE pipeline 입력 파일 확인
ls -lh de-pipeline/data/raw/*_counts.txt
ls -lh de-pipeline/data/metadata/*_metadata.csv
ls -lh de-pipeline/configs/config_*.yml

# Config 내용 확인
cat de-pipeline/configs/config_mouse-chd8.yml
```

---

## Phase 6: LLM 통합 (선택)

### 목표

자연어로 파이프라인과 상호작용 (데이터 보안 고려, 로컬 LLM 우선)

### 구현 개요

**참조**: RNA-seq 파이프라인의 완전한 구현 참조
- `scripts/standardization/llm_agent.py` (582 lines)
- `docs/user/OLLAMA_SETUP.md` (설정 가이드)
- `docs/user/LLM_AGENT_QUICKSTART.md` (빠른 참조)
- `docs/developer/PHASE6_LLM_INTEGRATION.md` (기술 문서)

### 핵심 개념

```python
class PipelineAgent:
    """자연어 파이프라인 에이전트"""
    
    def __init__(self, llm_provider="ollama"):
        # Ollama (로컬, 보안) 또는 OpenAI (클라우드)
        self.llm_client = initialize_llm(llm_provider)
        self.tools = self._define_tools()
    
    def _define_tools(self):
        """에이전트가 사용할 도구 정의"""
        return [
            {
                "name": "get_project_status",
                "description": "Get overall QC status and pass rates",
                "function": self._get_status
            },
            {
                "name": "compare_conditions",
                "description": "Compare QC metrics between conditions",
                "function": self._compare_conditions
            },
            # ... 추가 도구
        ]
    
    def chat(self, user_message: str) -> str:
        """자연어 쿼리 처리"""
        # LLM이 적절한 도구 선택
        # 도구 실행
        # 결과를 자연어로 변환하여 반환
```

### 사용 예시

```bash
# 대화형 모드
python scripts/standardization/llm_agent.py \
    --project-summary project_summary.json \
    --rnaseq-output /output \
    --interactive

# 예시 대화:
You: Show me the QC status
Agent: Your project has 38 samples with 100% pass rate...

You: Compare wildtype and knockout
Agent: Comparing conditions:
       - wildtype: 92.4% mapping rate
       - knockout: 91.8% mapping rate
       ...

You: Can you start DE analysis?
Agent: I'll prepare DE analysis. This will:
       1. Copy counts matrix
       2. Generate metadata
       3. Create config file
       Proceed? (yes/no)
```

### 보안 고려사항

**로컬 LLM (Ollama) 사용 권장**:
- ✅ 모든 데이터가 로컬 서버에 유지
- ✅ 외부 API 호출 없음
- ✅ HIPAA/GDPR 준수
- ✅ 완전한 감사 추적

**구현 시 주의사항**:
1. 민감한 샘플 ID를 LLM에 직접 노출하지 않도록 주의
2. 프로덕션 환경에서는 Ollama 사용 강제
3. 로그에 모든 쿼리 기록
4. 사용자 인증 및 권한 관리 추가

### Phase 6 설정

```bash
# 1. Ollama 설치
curl -fsSL https://ollama.com/install.sh | sh
ollama serve

# 2. 모델 다운로드
ollama pull llama3.1:8b

# 3. Python 패키지
pip install ollama

# 4. 테스트
python scripts/standardization/llm_agent.py \
    --project-summary project_summary.json \
    --message "Show QC status"
```

---

## 다른 파이프라인 적용 가이드

### ChIP-seq 파이프라인 표준화

**Phase 1-2 커스터마이징**:

```python
# ChIP-seq 특화 QC 메트릭
def extract_qc_metrics_chipseq(sample_id: str, output_dir: Path) -> dict:
    metrics = {}
    
    # MACS2 peak calling 결과
    peaks_file = output_dir / "final_outputs" / "peaks" / f"{sample_id}_peaks.narrowPeak"
    if peaks_file.exists():
        with open(peaks_file) as f:
            peak_count = sum(1 for line in f)
        metrics["peak_count"] = peak_count
    
    # Fraction of Reads in Peaks (FRiP)
    frip_file = output_dir / "metadata" / f"{sample_id}_frip.txt"
    if frip_file.exists():
        with open(frip_file) as f:
            metrics["frip"] = float(f.read().strip())
    
    # Cross-correlation (phantompeakqualtools)
    cc_file = output_dir / "metadata" / f"{sample_id}_cc.txt"
    if cc_file.exists():
        with open(cc_file) as f:
            parts = f.read().strip().split('\t')
            metrics["nsc"] = float(parts[8])  # Normalized Strand Coefficient
            metrics["rsc"] = float(parts[9])  # Relative Strand Coefficient
    
    return metrics

def evaluate_qc_status_chipseq(metrics: dict) -> tuple:
    """ChIP-seq QC 기준"""
    status = "PASS"
    issues = []
    
    if metrics.get("frip", 0) < 0.01:  # FRiP < 1%
        status = "FAIL"
        issues.append("Very low FRiP score")
    
    if metrics.get("peak_count", 0) < 1000:
        status = "WARN"
        issues.append("Low peak count")
    
    if metrics.get("nsc", 0) < 1.05:
        status = "WARN"
        issues.append("Low NSC (poor enrichment)")
    
    return status, issues
```

### ATAC-seq 파이프라인 표준화

```python
# ATAC-seq 특화 메트릭
def extract_qc_metrics_atacseq(sample_id: str, output_dir: Path) -> dict:
    metrics = {}
    
    # TSS enrichment
    tss_file = output_dir / "metadata" / f"{sample_id}_tss_enrich.txt"
    if tss_file.exists():
        with open(tss_file) as f:
            metrics["tss_enrichment"] = float(f.read().strip())
    
    # Fragment size distribution
    frag_file = output_dir / "metadata" / f"{sample_id}_fragment_sizes.txt"
    if frag_file.exists():
        # Nucleosome-free vs mono-nucleosome 비율 계산
        # ... (구현)
        pass
    
    # Mitochondrial reads %
    mito_file = output_dir / "metadata" / f"{sample_id}_mito_pct.txt"
    if mito_file.exists():
        with open(mito_file) as f:
            metrics["mitochondrial_pct"] = float(f.read().strip())
    
    return metrics

def evaluate_qc_status_atacseq(metrics: dict) -> tuple:
    """ATAC-seq QC 기준"""
    status = "PASS"
    issues = []
    
    if metrics.get("tss_enrichment", 0) < 5:
        status = "WARN"
        issues.append("Low TSS enrichment")
    
    if metrics.get("mitochondrial_pct", 0) > 30:
        status = "FAIL"
        issues.append("High mitochondrial contamination")
    
    return status, issues
```

### WGS/WES 파이프라인 표준화

```python
# Variant calling 메트릭
def extract_qc_metrics_wgs(sample_id: str, output_dir: Path) -> dict:
    metrics = {}
    
    # Coverage statistics
    cov_file = output_dir / "metadata" / f"{sample_id}_coverage.txt"
    if cov_file.exists():
        with open(cov_file) as f:
            for line in f:
                if "mean coverage" in line:
                    metrics["mean_coverage"] = float(line.split()[-1])
                if "pct_20x" in line:  # % bases at 20x
                    metrics["pct_20x"] = float(line.split()[-1])
    
    # Variant statistics
    vcf_stats = output_dir / "metadata" / f"{sample_id}_vcf_stats.txt"
    if vcf_stats.exists():
        with open(vcf_stats) as f:
            data = json.load(f)
            metrics["snp_count"] = data["snp_count"]
            metrics["indel_count"] = data["indel_count"]
            metrics["ti_tv_ratio"] = data["ti_tv_ratio"]
    
    return metrics

def evaluate_qc_status_wgs(metrics: dict) -> tuple:
    """WGS QC 기준"""
    status = "PASS"
    issues = []
    
    if metrics.get("mean_coverage", 0) < 30:
        status = "FAIL"
        issues.append("Low mean coverage")
    
    if metrics.get("pct_20x", 0) < 90:
        status = "WARN"
        issues.append("Low coverage uniformity")
    
    # Ti/Tv ratio should be ~2.0 for WGS
    ti_tv = metrics.get("ti_tv_ratio", 0)
    if ti_tv < 1.8 or ti_tv > 2.2:
        status = "WARN"
        issues.append("Unusual Ti/Tv ratio")
    
    return status, issues
```

### 공통 적용 패턴

모든 파이프라인에 동일하게 적용:

1. **Phase 1**: 3-tier 디렉토리 구조 (final_outputs, intermediate, metadata)
2. **Phase 2**: Manifest.json 생성 (파이프라인별 QC 메트릭만 다름)
3. **Phase 3**: Project summary 생성 (동일 스크립트, QC 기준만 조정)
4. **Phase 4**: Query agent (동일 스크립트)
5. **Phase 5**: Bridge to downstream (파이프라인별 커스터마이징)
6. **Phase 6**: LLM agent (동일 프레임워크, 도구만 추가)

---

## 체크리스트

### Phase 1: 데이터 구조화
- [ ] 표준 디렉토리 구조 정의
- [ ] Snakemake 규칙 수정 (final_outputs, intermediate, metadata)
- [ ] 디렉토리 자동 생성 스크립트 작성 (RNA-seq 참조)
- [ ] 기존 파일 마이그레이션 (필요시)
- [ ] 구조 검증 (tree 명령)

### Phase 2: 메타데이터 통합
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/generate_manifest.py`
- [ ] **[수정]** `extract_qc_metrics()` 함수 - 파이프라인별 로그 파싱
- [ ] **[수정]** `evaluate_qc_status()` 함수 - 파이프라인별 QC 기준
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/convert_sample_sheet.py` (수정 없음)
- [ ] Snakemake 규칙 추가
- [ ] 테스트 샘플로 검증

### Phase 3: 프로젝트 레벨 요약
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/generate_project_summary.py` (수정 없음)
- [ ] QC 상태 요약 함수
- [ ] 조건별 그룹화 함수
- [ ] 통계 지표 계산 함수
- [ ] 문제 샘플 식별 함수
- [ ] Project summary 생성 스크립트
- [ ] 전체 프로젝트로 테스트

### Phase 4: 쿼리 인터페이스
- [ ] 쿼리 타입 정의 (최소 8개)
- [ ] Query agent 클래스 작성
- [ ] 각 쿼리 핸들러 구현
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/generate_project_summary.py` (수정 없음)
- [ ] 전체 프로젝트로 테스트

### Phase 4: 쿼리 인터페이스
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/agent_query.py` (수정 없음)
- [ ] 모든 쿼리 타입 테스트 (status, failed, condition, compare, sample, stats, issues, conditions)
- [ ] 문서화 (사용 예시)

### Phase 5: Pipeline 연결
- [ ] **[복사]** `rna-seq-pipeline/scripts/bridge_to_de_pipeline.py` (템플릿으로 사용)
- [ ] **[수정]** Downstream 파이프라인 입력 형식에 맞춰 조정
- [ ] **[수정]** `prepare_de_input()` - 파일 복사 로직
- [ ] **[수정]** `generate_metadata_template()` - 메타데이터 형식
- [ ] **[수정]** `generate_de_config()` - Config 생성 로직
- [ ] Dry-run 테스트
- [ ] End-to-end 테스트

### Phase 6: LLM 통합 (선택)
- [ ] **[복사]** `rna-seq-pipeline/scripts/standardization/llm_agent.py` (수정 없음)
- [ ] Ollama 설치 및 설정
- [ ] 모델 다운로드 (llama3.1:8b)
- [ ] 테스트 쿼리 실행
- [ ] **[선택]** 파이프라인 특화 도구 함수 추가

### 재사용 코드 확인
- [ ] **그대로 복사 (5개)**: convert_sample_sheet.py, update_manifests_metadata.py, generate_project_summary.py, agent_query.py, llm_agent.py
- [ ] **함수 수정 (1개)**: generate_manifest.py (extract_qc_metrics, evaluate_qc_status)
- [ ] **템플릿 사용 (1개)**: bridge_to_de_pipeline.py (Downstream에 맞춰 수정)

### 문서화
- [ ] 파이프라인별 QC 메트릭 문서화
- [ ] Downstream 연결 가이드 작성
- [ ] 사용자 가이드 (각 Phase별)
- [ ] README 업데이트

### 배포 및 유지보수
- [ ] Git 저장소 정리
- [ ] 버전 태깅
- [ ] 팀 교육 자료 준비
- [ ] 파일럿 프로젝트 적용
- [ ] 피드백 수집
- [ ] 개선 사항 반영
- [ ] 다른 파이프라인으로 확장

---

## 예상 시간표

| Phase | 작업 시간 | 검증 시간 | 문서화 | 총 |
|-------|----------|----------|--------|-----|
| Phase 1 | 4시간 | 1시간 | 2시간 | 7시간 |
| Phase 2 | 8시간 | 2시간 | 3시간 | 13시간 |
| Phase 3 | 6시간 | 1시간 | 2시간 | 9시간 |
| Phase 4 | 6시간 | 2시간 | 2시간 | 10시간 |
| Phase 5 | 8시간 | 2시간 | 3시간 | 13시간 |
| Phase 6 | 12시간 | 3시간 | 4시간 | 19시간 |
| **총계** | **44시간** | **11시간** | **16시간** | **71시간** |

*Phase 6은 선택사항이므로 제외 시 52시간*

---

## 성공 지표

### 기술적 지표
- ✅ 모든 샘플에 manifest.json 생성
- ✅ Project summary 자동 생성
- ✅ 쿼리 에이전트 정상 작동 (8개 쿼리 타입)
- ✅ Pipeline 브리지 성공률 100%
- ✅ QC 자동 평가 정확도 > 95%

### 효율성 지표
- ⏱️ QC 상태 확인 시간: 수동 30분 → 자동 10초
- ⏱️ Pipeline 연결 시간: 수동 2시간 → 자동 5분
- ⏱️ 프로젝트 상태 파악: 수동 1시간 → 즉시

### 사용자 만족도
- 👍 팀원이 표준화된 구조를 쉽게 이해
- 👍 새 프로젝트 시작이 빨라짐
- 👍 수동 작업 최소화로 오류 감소
- 👍 재현성 향상

---

## 빠른 시작: 새 파이프라인에 표준화 적용하기

### 단계별 간단 가이드

#### Step 1: 스크립트 복사 (5분)

```bash
# 새 파이프라인 디렉토리로 이동
cd /path/to/new-pipeline

# standardization 디렉토리 생성
mkdir -p scripts/standardization

# RNA-seq 파이프라인에서 재사용 가능한 스크립트 복사
cp /path/to/rna-seq-pipeline/scripts/standardization/convert_sample_sheet.py \
   scripts/standardization/

cp /path/to/rna-seq-pipeline/scripts/standardization/generate_manifest.py \
   scripts/standardization/

cp /path/to/rna-seq-pipeline/scripts/standardization/update_manifests_metadata.py \
   scripts/standardization/

cp /path/to/rna-seq-pipeline/scripts/standardization/generate_project_summary.py \
   scripts/standardization/

cp /path/to/rna-seq-pipeline/scripts/standardization/agent_query.py \
   scripts/standardization/

cp /path/to/rna-seq-pipeline/scripts/standardization/llm_agent.py \
   scripts/standardization/  # Optional
```

#### Step 2: generate_manifest.py 수정 (30분)

```bash
# 텍스트 에디터로 열기
vim scripts/standardization/generate_manifest.py

# 수정할 부분:
# 1. extract_qc_metrics() 함수
#    - 함수 이름: extract_qc_metrics_YOURPIPELINE()
#    - 로그 파일 경로를 파이프라인에 맞게 수정
#    - QC 메트릭 파싱 로직 추가

# 2. evaluate_qc_status() 함수
#    - 함수 이름: evaluate_qc_status_YOURPIPELINE()
#    - QC threshold를 파이프라인에 맞게 조정

# 3. main() 함수
#    - 함수 호출 부분만 변경
#    manifest["qc_metrics"] = extract_qc_metrics_YOURPIPELINE(...)
#    manifest["qc_status"], manifest["issues"] = evaluate_qc_status_YOURPIPELINE(...)
```

#### Step 3: 샘플 시트 준비 (10분)

```bash
# 템플릿 생성
python scripts/standardization/convert_sample_sheet.py --create-template samples.csv

# 샘플 정보 입력 (Excel/텍스트 에디터)
# 필수 컬럼: sample_id, condition, replicate, tissue, age, sex

# 파이프라인 형식으로 변환
python scripts/standardization/convert_sample_sheet.py samples.csv -o config/
```

#### Step 4: 테스트 (10분)

```bash
# 1. Manifest 생성 테스트
python scripts/standardization/generate_manifest.py \
    --sample-id TEST_SAMPLE \
    --output-dir /path/to/output \
    --sample-sheet config/samples.tsv

# 2. Project summary 생성 테스트
python scripts/standardization/generate_project_summary.py \
    --project-dir /path/to/output \
    --project-id test-project

# 3. Query 테스트
python scripts/standardization/agent_query.py \
    --project-summary /path/to/output/project_summary.json \
    --query status
```

#### Step 5: 완료! 🎉

이제 파이프라인이 표준화되었습니다. 

**재사용률**: ~85% (7개 파일 중 6개 그대로 사용, 1개만 함수 수정)

---

## 재사용 가능한 파일 요약

| 파일 | 재사용 방법 | 수정 필요 | 시간 |
|------|------------|----------|------|
| `convert_sample_sheet.py` | ✅ 그대로 복사 | 없음 | 1분 |
| `update_manifests_metadata.py` | ✅ 그대로 복사 | 없음 | 1분 |
| `generate_manifest.py` | ⚠️ 복사 후 수정 | 함수 2개 | 30분 |
| `generate_project_summary.py` | ✅ 그대로 복사 | 없음 | 1분 |
| `agent_query.py` | ✅ 그대로 복사 | 없음 | 1분 |
| `llm_agent.py` | ✅ 그대로 복사 | 없음 (선택) | 1분 |
| `bridge_to_de_pipeline.py` | ⚠️ 템플릿 사용 | 전체 조정 | 2시간 |

**총 예상 시간**: 약 3시간 (Bridge 포함), 1시간 (Bridge 제외)

---

## 결론

이 가이드를 따르면 모든 NGS 파이프라인에 일관된 표준화를 적용할 수 있습니다.

**핵심 원칙**:
1. **재사용 우선**: RNA-seq 코드를 최대한 활용
2. **최소 수정**: 파이프라인 특화 부분만 수정
3. **일관성**: 모든 파이프라인이 동일한 구조
4. **자동화**: 수동 작업 최소화
5. **보안**: 민감한 데이터 로컬 관리

**다음 단계**:
1. ✅ **RNA-seq 파이프라인 코드 복사** - 위의 "빠른 시작" 참조
2. ⚠️ **generate_manifest.py 함수 2개만 수정** - 30분 작업
3. ✅ **나머지는 그대로 사용** - 추가 작업 없음
4. 📝 **테스트 및 검증**
5. 📚 **문서화 및 팀 공유**

**참고 자료**:
- **원본 파일**: `rna-seq-pipeline/scripts/standardization/` 디렉토리
- **기술 문서**: `docs/developer/PHASE6_LLM_INTEGRATION.md`
- **사용자 가이드**: `docs/user/OLLAMA_SETUP.md`
- **진행 상황**: `docs/developer/STANDARDIZATION_PROGRESS.md`

**핵심 메시지**: 
> 💡 **새로 작성하지 마세요!** RNA-seq 파이프라인의 코드를 복사하고 필요한 부분만 수정하세요. 
> 85%의 코드는 그대로 재사용 가능합니다.

표준화를 통해 연구 효율성을 극대화하고 재현 가능한 분석 환경을 구축하세요! 🚀
