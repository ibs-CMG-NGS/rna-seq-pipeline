# Pipeline Standardization Implementation

## ✅ 완료 사항

### 1. 마스터 샘플 시트 표준 정의
- ✅ `STANDARDIZATION.md` - 전체 표준화 가이드 작성
- ✅ `samples_master.csv` - H2O2 프로젝트용 마스터 시트 생성
- ✅ 필수 컬럼 17개 정의
- ✅ 프로젝트/샘플/실험 메타데이터 통합

### 2. 샘플 시트 변환 도구
- ✅ `scripts/convert_sample_sheet.py` 작성
- ✅ Snakemake TSV 변환 기능
- ✅ WDL JSON 변환 기능 (스켈레톤)
- ✅ Nextflow CSV 변환 기능 (스켈레톤)
- ✅ 템플릿 생성 기능

### 3. 사용 예제
```bash
# 템플릿 생성
python scripts/convert_sample_sheet.py --create-template my_samples.csv

# 모든 파이프라인 형식으로 변환
python scripts/convert_sample_sheet.py samples_master.csv -o config/

# Snakemake만 변환
python scripts/convert_sample_sheet.py samples_master.csv --snakemake samples.tsv

# 특정 프로젝트만 변환
python scripts/convert_sample_sheet.py samples_master.csv -o config/ --project H2O2_human_2025
```

## 📋 다음 단계 (Phase 1 완료를 위해)

### 출력 디렉토리 구조 표준화
- [ ] Snakefile 수정: 표준 디렉토리 구조로 출력
- [ ] `final_outputs/` vs `intermediate/` 분리
- [ ] `manifest.json` 생성 규칙 추가
- [ ] `qc_summary.json` 생성 규칙 추가

### 결과 수집 도구
- [ ] `scripts/generate_manifest.py` - manifest.json 생성
- [ ] `scripts/collect_results.py` - 최종 결과물 수집
- [ ] `scripts/validate_outputs.py` - 파일 무결성 검증

### 마이그레이션 도구
- [ ] `scripts/migrate_to_standard.py` - 기존 결과 → 표준 구조 변환
- [ ] H2O2 프로젝트 마이그레이션 테스트

## 🎯 표준 디렉토리 구조 (목표)

```
/home/ngs/data/results/
└── H2O2_human_2025/
    ├── metadata/
    │   ├── samples_master.csv
    │   ├── analysis_log.json
    │   └── pipeline_config.yaml
    │
    ├── h_RNA_Cont_1/
    │   └── rna-seq/
    │       ├── final_outputs/
    │       │   ├── bam/
    │       │   │   ├── aligned.sorted.bam
    │       │   │   └── aligned.sorted.bam.bai
    │       │   ├── counts/
    │       │   │   └── gene_counts.csv
    │       │   ├── qc/
    │       │   │   ├── multiqc_report.html
    │       │   │   └── qc_summary.json
    │       │   └── manifest.json
    │       │
    │       ├── intermediate/
    │       │   ├── trimmed/
    │       │   ├── fastqc/
    │       │   └── logs/
    │       │
    │       └── metadata/
    │           └── execution_time.json
    │
    └── project_summary/
        ├── multiqc_report.html
        ├── combined_counts.csv
        └── de_analysis/
```

## 📊 manifest.json 스키마 (목표)

```json
{
  "sample_id": "h_RNA_Cont_1",
  "project_id": "H2O2_human_2025",
  "pipeline_type": "rna-seq",
  "pipeline_version": "1.0.0",
  "execution_date": "2026-02-09",
  "status": "completed",
  
  "final_outputs": {
    "aligned_bam": {
      "path": "bam/aligned.sorted.bam",
      "md5": "...",
      "size_bytes": 1234567890
    },
    "gene_counts": {
      "path": "counts/gene_counts.csv",
      "md5": "...",
      "size_bytes": 123456
    }
  },
  
  "qc_metrics": {
    "overall_status": "PASS",
    "total_reads": 50000000,
    "mapping_rate": 0.94,
    "assignment_rate": 0.85
  }
}
```

## 🚀 적용 순서

1. **현재 완료**: 샘플 시트 표준화 ✅
2. **다음**: 출력 구조 표준화
3. **그 다음**: 결과 수집 도구
4. **마지막**: WGS, ATAC-seq 확장
