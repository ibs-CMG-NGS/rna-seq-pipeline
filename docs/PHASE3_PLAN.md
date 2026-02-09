# Phase 3: Snakefile Integration Plan

## 🎯 목표
표준 디렉토리 구조를 Snakefile에 완전히 통합하여 각 샘플이 독립적인 final_outputs, intermediate, metadata 디렉토리를 갖도록 함.

---

## 📂 현재 구조 vs 표준 구조

### 현재 Legacy 구조
```
results/
├── trimmed/
│   ├── sample1_1.fastq.gz
│   └── sample1_2.fastq.gz
├── aligned/
│   └── sample1/
│       ├── Aligned.sortedByCoord.out.bam
│       └── Log.final.out
├── counts/
│   └── counts_matrix.txt
└── qc/
    ├── sample1_1_fastqc.html
    └── sample1_2_fastqc.html
```

### 표준 구조 (USE_STANDARD=true)
```
/home/ngs/data/results/
└── H2O2_human_2025/                    # PROJECT_ID
    ├── h_RNA_200_1/                    # SAMPLE_ID
    │   └── rna-seq/                    # PIPELINE_TYPE
    │       ├── final_outputs/
    │       │   ├── bam/
    │       │   │   ├── aligned.sorted.bam
    │       │   │   └── aligned.sorted.bam.bai
    │       │   ├── counts/
    │       │   │   └── gene_counts.csv
    │       │   ├── qc/
    │       │   │   └── qc_summary.json
    │       │   └── manifest.json
    │       ├── intermediate/
    │       │   ├── trimmed/
    │       │   │   ├── sample_1.fastq.gz
    │       │   │   └── sample_2.fastq.gz
    │       │   ├── fastqc/
    │       │   │   ├── raw_1_fastqc.html
    │       │   │   └── raw_2_fastqc.html
    │       │   └── logs/
    │       │       ├── cutadapt.log
    │       │       └── star.log
    │       └── metadata/
    │           └── sample_info.yaml
    ├── h_RNA_200_2/
    │   └── rna-seq/
    │       └── ...
    ├── project_summary/                # 프로젝트 전체 요약
    │   ├── counts/
    │   │   └── counts_matrix_all.txt   # 전체 샘플 counts
    │   ├── qc/
    │   │   └── multiqc_report.html
    │   └── differential_expression/
    ├── metadata/
    │   └── samples_master.csv
    └── logs/
        └── snakemake.log
```

---

## 🔧 구현 단계

### Step 1: 규칙별 output 경로 수정 ✅

#### 1.1 fastqc_raw (intermediate)
```python
# Before (Legacy)
output:
    html=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.html",
    zip=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.zip"

# After (Standard)
output:
    html=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/fastqc/{{sample}}_{{read}}_fastqc.html",
    zip=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/fastqc/{{sample}}_{{read}}_fastqc.zip"
```

#### 1.2 cutadapt (intermediate)
```python
# Before
output:
    r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
    r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz"

# After
output:
    r1=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/trimmed/{{sample}}_1.fastq.gz",
    r2=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/trimmed/{{sample}}_2.fastq.gz"
```

#### 1.3 star_align (final_outputs)
```python
# Before
output:
    bam=f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam",
    log_final=f"{ALIGNED_DIR}/{{sample}}/Log.final.out"

# After
output:
    bam=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam",
    bai=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam.bai",
    log_final=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/logs/star_final.log"
```

#### 1.4 featurecounts_quant (샘플별 counts + 프로젝트 전체)
```python
# 샘플별 counts (final_outputs)
rule featurecounts_sample:
    input:
        bam=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam"
    output:
        counts=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/counts/gene_counts.txt"

# 프로젝트 전체 counts matrix (project_summary)
rule featurecounts_all:
    input:
        bams=expand(lambda wildcards: f"{get_final_outputs_dir('{{sample}}')}/bam/aligned.sorted.bam", sample=SAMPLES)
    output:
        counts=f"{PROJECT_SUMMARY_DIR}/counts/counts_matrix_all.txt"
```

---

### Step 2: 새 규칙 추가 ✅

#### 2.1 BAM index 생성
```python
rule index_bam:
    input:
        bam=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam"
    output:
        bai=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam.bai"
    shell:
        "samtools index {input.bam}"
```

#### 2.2 QC summary 생성
```python
rule generate_qc_summary:
    input:
        star_log=lambda wildcards: f"{get_intermediate_dir(wildcards.sample)}/logs/star_final.log",
        fc_summary=f"{PROJECT_SUMMARY_DIR}/counts/counts_matrix_all.txt.summary"
    output:
        qc_json=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/qc/qc_summary.json"
    params:
        sample_id="{sample}"
    shell:
        """
        python3 scripts/generate_qc_summary.py \
            --sample-id {params.sample_id} \
            --star-log {input.star_log} \
            --featurecounts {input.fc_summary} \
            -o {output.qc_json}
        """
```

#### 2.3 Manifest 생성
```python
rule generate_manifest:
    input:
        bam=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam",
        bai=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam.bai",
        counts=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/counts/gene_counts.txt",
        qc_summary=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/qc/qc_summary.json"
    output:
        manifest=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/manifest.json"
    params:
        sample_id="{sample}",
        sample_dir=lambda wildcards: get_sample_dir(wildcards.sample)
    shell:
        """
        python3 scripts/generate_manifest.py \
            --sample-dir {params.sample_dir} \
            --sample-id {params.sample_id} \
            --project-id {PROJECT_ID} \
            --pipeline-type {PIPELINE_TYPE}
        """
```

---

### Step 3: rule all 업데이트 ✅

```python
rule all:
    input:
        # 표준 구조 사용 시
        (
            # 샘플별 manifest
            expand(lambda wildcards: f"{get_final_outputs_dir('{{sample}}')}/manifest.json", 
                   sample=SAMPLES),
            # 프로젝트 전체 요약
            f"{PROJECT_SUMMARY_DIR}/counts/counts_matrix_all.txt",
            f"{PROJECT_SUMMARY_DIR}/qc/multiqc_report.html"
        ) if USE_STANDARD else (
            # Legacy 구조
            expand(f"{QC_DIR}/{{sample}}_{{read}}_fastqc.html", sample=SAMPLES, read=[1, 2]),
            expand(f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
            f"{COUNTS_DIR}/counts_matrix.txt",
            f"{RESULTS_DIR}/multiqc_report.html" if config.get("generate_multiqc", True) else []
        )
```

---

### Step 4: 로그 디렉토리 통합 ✅

**현재**: `logs/fastqc/`, `logs/cutadapt/`, `logs/star/`

**표준 구조**:
- 샘플별 로그: `{SAMPLE_DIR}/intermediate/logs/`
- 프로젝트 전체 로그: `{PROJECT_DIR}/logs/`

```python
def get_sample_log(wildcards, tool):
    """샘플별 로그 파일 경로"""
    if USE_STANDARD:
        return f"{get_intermediate_dir(wildcards.sample)}/logs/{tool}.log"
    else:
        return f"{LOGS_DIR}/{tool}/{wildcards.sample}.log"
```

---

## 🧪 테스트 전략

### 1. 단일 샘플 테스트
```bash
# 표준 구조로 단일 샘플 실행
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --forceall \
    --cores 4 \
    -- h_RNA_200_1/rna-seq/final_outputs/manifest.json
```

### 2. 표준 구조 검증
```bash
# 디렉토리 구조 확인
tree -L 5 /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/

# Manifest 검증
python3 scripts/generate_manifest.py \
    --validate /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/manifest.json
```

### 3. 전체 파이프라인 테스트
```bash
# 모든 샘플 실행
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --cores 12
```

---

## 📋 체크리스트

Phase 3 구현 완료 확인:

- [ ] Step 1: 규칙별 output 경로 수정
  - [ ] fastqc_raw → intermediate/fastqc/
  - [ ] cutadapt → intermediate/trimmed/
  - [ ] star_align → final_outputs/bam/
  - [ ] featurecounts → final_outputs/counts/ + project_summary/

- [ ] Step 2: 새 규칙 추가
  - [ ] index_bam (samtools index)
  - [ ] generate_qc_summary
  - [ ] generate_manifest

- [ ] Step 3: rule all 업데이트
  - [ ] 표준 구조 타겟
  - [ ] Legacy 구조 호환

- [ ] Step 4: 로그 디렉토리 통합
  - [ ] 샘플별 로그 → intermediate/logs/
  - [ ] 프로젝트 로그 → PROJECT_DIR/logs/

- [ ] Step 5: 테스트
  - [ ] 단일 샘플 테스트
  - [ ] Manifest 검증
  - [ ] 전체 파이프라인 실행

---

## 🚀 구현 우선순위

### High Priority (즉시 구현)
1. **BAM index 생성** - samtools index rule 추가
2. **QC summary 생성** - generate_qc_summary rule 추가  
3. **Manifest 생성** - generate_manifest rule 추가
4. **rule all 업데이트** - 표준 구조 타겟 추가

### Medium Priority (다음 단계)
5. **featureCounts 샘플별 분리** - 샘플별 counts + 프로젝트 전체 counts
6. **로그 디렉토리 통합** - intermediate/logs/ 이동
7. **FastQC 경로 수정** - intermediate/fastqc/ 이동

### Low Priority (선택 사항)
8. **Metadata 파일 생성** - 샘플 정보 YAML/JSON
9. **프로젝트 요약 리포트** - 전체 샘플 QC 통합 리포트

---

## 💡 구현 시 주의사항

1. **Lambda 함수 사용**: output 경로에서 wildcards를 사용하려면 lambda 필수
   ```python
   output:
       bam=lambda wildcards: f"{get_final_outputs_dir(wildcards.sample)}/bam/aligned.sorted.bam"
   ```

2. **디렉토리 자동 생성**: Snakemake는 output 디렉토리를 자동 생성하므로 `os.makedirs()` 불필요

3. **Legacy 호환성**: `USE_STANDARD=False`일 때도 기존 규칙이 작동하도록 조건부 처리

4. **featureCounts 처리**: 전체 샘플 counts는 한 번만 계산 (project_summary/)

5. **Dependency chain**: 
   ```
   fastq → cutadapt → star_align → index_bam → featurecounts → generate_qc_summary → generate_manifest
   ```

