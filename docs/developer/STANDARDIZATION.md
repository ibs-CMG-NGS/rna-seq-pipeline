# Pipeline Standardization Guide

## 목적
에이전트 기반 파이프라인 관리를 위한 입출력 표준화 가이드

## 1. 샘플 시트 표준화

### 1.1 마스터 샘플 시트 (`samples.csv`)

모든 파이프라인(RNA-seq, WGS, ATAC-seq)에서 공통으로 사용하는 최상위 샘플 정의 파일

#### 필수 컬럼:
- `project_id`: 프로젝트 식별자
- `sample_id`: 고유한 샘플 ID
- `sample_name`: 샘플 이름 (분석 결과 표시용)
- `condition`: 실험 조건/그룹
- `replicate`: 생물학적 반복 번호
- `sequencing_platform`: 시퀀싱 플랫폼 (Illumina, PacBio 등)
- `library_type`: 라이브러리 타입 (RNA-seq, WGS, ATAC-seq)
- `read_type`: single-end 또는 paired-end
- `fastq_1`: Read 1 FASTQ 파일 경로 (절대경로)
- `fastq_2`: Read 2 FASTQ 파일 경로 (paired-end인 경우)
- `species`: 생물종 (human, mouse 등)
- `genome_build`: 참조 유전체 버전 (GRCh38, mm10 등)

#### 선택 컬럼:
- `batch`: 배치 효과 보정용
- `sex`: 성별
- `age`: 나이
- `tissue`: 조직 타입
- `treatment`: 처리 조건
- `time_point`: 시간대
- `notes`: 추가 메모

#### 예시 (`samples.csv`):
```csv
project_id,sample_id,sample_name,condition,replicate,sequencing_platform,library_type,read_type,fastq_1,fastq_2,species,genome_build,tissue,treatment,time_point,notes
H2O2_2025,H2O2_Cont_1,Control_Rep1,Control,1,Illumina,RNA-seq,paired-end,/data/fastq/Cont_1_R1.fq.gz,/data/fastq/Cont_1_R2.fq.gz,human,GRCh38,fibroblast,none,0h,Control replicate 1
H2O2_2025,H2O2_100_1,H2O2_100uM_Rep1,H2O2_100uM,1,Illumina,RNA-seq,paired-end,/data/fastq/100_1_R1.fq.gz,/data/fastq/100_1_R2.fq.gz,human,GRCh38,fibroblast,100uM_H2O2,24h,100uM H2O2 treatment
H2O2_2025,H2O2_200_1,H2O2_200uM_Rep1,H2O2_200uM,1,Illumina,RNA-seq,paired-end,/data/fastq/200_1_R1.fq.gz,/data/fastq/200_1_R2.fq.gz,human,GRCh38,fibroblast,200uM_H2O2,24h,200uM H2O2 treatment
```

### 1.2 파이프라인별 변환

에이전트는 `samples.csv`를 읽어서 각 파이프라인 형식으로 변환:

#### RNA-seq Snakemake → `samples.tsv`
```tsv
sample_id	condition	replicate	fastq_r1	fastq_r2	notes
H2O2_Cont_1	Control	1	/data/fastq/Cont_1_R1.fq.gz	/data/fastq/Cont_1_R2.fq.gz	Control replicate 1
```

#### WGS WDL → `inputs.json`
```json
{
  "samples": [
    {
      "sample_id": "H2O2_Cont_1",
      "fastq_1": "/data/fastq/Cont_1_R1.fq.gz",
      "fastq_2": "/data/fastq/Cont_1_R2.fq.gz"
    }
  ]
}
```

#### ATAC-seq Nextflow → `samplesheet.csv`
```csv
sample,fastq_1,fastq_2,replicate
H2O2_Cont_1,/data/fastq/Cont_1_R1.fq.gz,/data/fastq/Cont_1_R2.fq.gz,1
```

---

## 2. 출력 디렉토리 표준화

### 2.1 디렉토리 구조

```
{base_results_dir}/
├── {project_id}/                    # 프로젝트별 최상위 디렉토리
│   ├── metadata/                    # 메타데이터
│   │   ├── samples.csv              # 마스터 샘플 시트
│   │   ├── analysis_log.json        # 분석 이력
│   │   └── pipeline_config.yaml     # 파이프라인 설정
│   │
│   ├── {sample_id}/                 # 샘플별 디렉토리
│   │   ├── {pipeline_type}/         # 파이프라인 타입 (rna-seq, wgs, atac-seq)
│   │   │   │
│   │   │   ├── final_outputs/       # ⭐ 최종 결과물 (에이전트가 읽을 파일)
│   │   │   │   ├── counts/
│   │   │   │   │   └── gene_counts.csv
│   │   │   │   ├── bam/
│   │   │   │   │   └── aligned.sorted.bam
│   │   │   │   │   └── aligned.sorted.bam.bai
│   │   │   │   ├── qc/
│   │   │   │   │   ├── multiqc_report.html
│   │   │   │   │   └── qc_summary.json
│   │   │   │   └── manifest.json    # 최종 결과물 목록
│   │   │   │
│   │   │   ├── intermediate/        # 중간 파일 (임시 파일)
│   │   │   │   ├── trimmed/
│   │   │   │   ├── fastqc/
│   │   │   │   └── logs/
│   │   │   │
│   │   │   └── metadata/
│   │   │       ├── pipeline_version.txt
│   │   │       ├── execution_time.json
│   │   │       └── resources_used.json
│   │   │
│   │   └── {another_pipeline_type}/
│   │
│   └── project_summary/              # 프로젝트 전체 요약
│       ├── combined_counts.csv       # 모든 샘플 통합 counts
│       ├── multiqc_report.html       # 전체 샘플 QC
│       ├── pca_plot.pdf              # PCA 분석
│       └── de_analysis/              # Differential expression 결과
│
└── {another_project_id}/
```

### 2.2 `manifest.json` 표준

각 샘플의 `final_outputs/manifest.json`은 에이전트가 읽을 파일 목록을 정의:

```json
{
  "sample_id": "H2O2_Cont_1",
  "project_id": "H2O2_2025",
  "pipeline_type": "rna-seq",
  "pipeline_version": "1.0.0",
  "execution_date": "2026-02-09T10:30:00",
  "status": "completed",
  "final_outputs": {
    "aligned_bam": {
      "path": "bam/aligned.sorted.bam",
      "md5": "a1b2c3d4e5f6...",
      "size_bytes": 1234567890,
      "description": "Sorted and indexed BAM file"
    },
    "gene_counts": {
      "path": "counts/gene_counts.csv",
      "md5": "f6e5d4c3b2a1...",
      "size_bytes": 123456,
      "description": "Gene-level counts matrix"
    },
    "qc_report": {
      "path": "qc/multiqc_report.html",
      "md5": "1a2b3c4d5e6f...",
      "size_bytes": 234567,
      "description": "Comprehensive QC report"
    },
    "qc_summary": {
      "path": "qc/qc_summary.json",
      "md5": "6f5e4d3c2b1a...",
      "size_bytes": 12345,
      "description": "Machine-readable QC metrics"
    }
  },
  "qc_metrics": {
    "total_reads": 50000000,
    "mapped_reads": 47000000,
    "mapping_rate": 0.94,
    "assigned_reads": 40000000,
    "assignment_rate": 0.85,
    "median_quality": 35,
    "status": "PASS"
  },
  "next_steps": [
    "differential_expression",
    "pathway_analysis"
  ]
}
```

---

## 3. 파이프라인별 출력 매핑

### 3.1 RNA-seq

**Final Outputs:**
- `bam/aligned.sorted.bam` - STAR alignment 결과
- `bam/aligned.sorted.bam.bai` - BAM index
- `counts/gene_counts.csv` - Gene-level counts
- `counts/transcript_counts.csv` - Transcript-level counts (optional)
- `qc/multiqc_report.html` - QC 리포트
- `qc/qc_summary.json` - QC 메트릭

**Intermediate:**
- `trimmed/*.fastq.gz` - Trimmed reads
- `fastqc/*_fastqc.html` - FastQC 결과
- `logs/*.log` - 실행 로그

### 3.2 WGS

**Final Outputs:**
- `bam/aligned.sorted.dedup.bam` - Aligned & deduplicated BAM
- `vcf/variants.filtered.vcf.gz` - Filtered variants
- `vcf/variants.annotated.vcf.gz` - Annotated variants
- `qc/multiqc_report.html` - QC 리포트
- `qc/coverage_summary.txt` - Coverage statistics

**Intermediate:**
- `trimmed/*.fastq.gz`
- `bam/aligned.sorted.bam` - Before deduplication
- `logs/*.log`

### 3.3 ATAC-seq

**Final Outputs:**
- `bam/aligned.filtered.bam` - Filtered BAM
- `peaks/peaks.narrowPeak` - Called peaks
- `bigwig/coverage.bw` - Coverage track
- `qc/multiqc_report.html` - QC 리포트
- `qc/fragment_length.pdf` - Fragment size distribution

**Intermediate:**
- `trimmed/*.fastq.gz`
- `bam/aligned.sorted.bam` - Before filtering
- `logs/*.log`

---

## 4. QC 메트릭 표준화

모든 파이프라인의 `qc_summary.json`은 공통 필드를 포함:

```json
{
  "sample_id": "H2O2_Cont_1",
  "pipeline_type": "rna-seq",
  "qc_version": "1.0",
  "overall_status": "PASS",  // PASS, WARN, FAIL
  
  "sequencing": {
    "total_reads": 50000000,
    "read_length": 150,
    "median_quality": 35,
    "q30_percentage": 92.5
  },
  
  "alignment": {
    "mapped_reads": 47000000,
    "mapping_rate": 0.94,
    "uniquely_mapped": 45000000,
    "multi_mapped": 2000000
  },
  
  "pipeline_specific": {
    // RNA-seq: assignment_rate, detected_genes
    // WGS: coverage, variant_count
    // ATAC-seq: peak_count, frip_score
  },
  
  "issues": [
    {
      "severity": "warning",
      "category": "duplication",
      "message": "High duplication rate detected (15%)"
    }
  ],
  
  "recommendations": [
    "Sample quality is sufficient for downstream analysis"
  ]
}
```

---

## 5. 에이전트 인터페이스

### 5.1 샘플 시트 생성기

```python
# scripts/generate_pipeline_inputs.py
"""
마스터 samples.csv를 읽어 각 파이프라인 형식으로 변환
"""

def convert_to_snakemake_tsv(samples_csv, output_tsv):
    """Snakemake용 TSV 생성"""
    pass

def convert_to_wdl_json(samples_csv, output_json):
    """WDL용 JSON 생성"""
    pass

def convert_to_nextflow_csv(samples_csv, output_csv):
    """Nextflow용 CSV 생성"""
    pass
```

### 5.2 결과 수집기

```python
# scripts/collect_results.py
"""
표준화된 디렉토리에서 최종 결과물만 수집
"""

def collect_final_outputs(project_id, sample_id, pipeline_type):
    """manifest.json을 읽어 최종 결과물 경로 반환"""
    pass

def validate_outputs(manifest_path):
    """MD5 체크섬으로 파일 무결성 검증"""
    pass

def prepare_for_downstream(project_id):
    """3차 분석용 입력 준비"""
    pass
```

### 5.3 QC 검증기

```python
# scripts/validate_qc.py
"""
모든 샘플의 QC 상태 검증
"""

def check_sample_qc(qc_summary_json):
    """QC 통과 여부 확인"""
    pass

def find_failed_samples(project_id):
    """FAIL 샘플 리스트 반환"""
    pass
```

---

## 6. 구현 로드맵

### Phase 1: RNA-seq 파이프라인 표준화 (현재)
- [ ] 마스터 `samples.csv` 형식 정의
- [ ] 출력 디렉토리 구조 변경
- [ ] `manifest.json` 생성 규칙 추가
- [ ] `qc_summary.json` 생성
- [ ] 샘플 시트 변환 스크립트 작성

### Phase 2: WGS/ATAC-seq 확장
- [ ] WGS 파이프라인 적용
- [ ] ATAC-seq 파이프라인 적용
- [ ] 공통 라이브러리 추출

### Phase 3: 에이전트 통합
- [ ] 에이전트용 API 설계
- [ ] 자동 QC 검증 워크플로우
- [ ] 3차 분석 자동 트리거

---

## 7. 마이그레이션 가이드

### 기존 파이프라인 → 표준 구조

```bash
# 1. 기존 결과를 표준 구조로 이동
python scripts/migrate_to_standard.py \
    --old-results /old/results \
    --new-base /standard/results \
    --project-id H2O2_2025

# 2. Manifest 자동 생성
python scripts/generate_manifests.py \
    --project-id H2O2_2025

# 3. 검증
python scripts/validate_structure.py \
    --project-id H2O2_2025
```

---

## 8. 예제: Human H2O2 프로젝트

현재 프로젝트를 표준 구조로 변환:

**기존:**
```
/home/ngs/data/h-rna-seq-pipeline-results/
├── trimmed/
├── aligned/
├── counts/
├── qc/
└── logs/
```

**표준:**
```
/home/ngs/data/results/
└── H2O2_2025/
    ├── metadata/
    │   └── samples.csv
    ├── h_RNA_Cont_1/
    │   └── rna-seq/
    │       ├── final_outputs/
    │       ├── intermediate/
    │       └── metadata/
    ├── h_RNA_100_1/
    │   └── rna-seq/
    │       └── ...
    └── project_summary/
        └── multiqc_report.html
```
