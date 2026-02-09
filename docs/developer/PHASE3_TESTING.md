# Phase 3: Snakefile Integration Testing Guide

## 🎯 구현 완료 항목

### ✅ High Priority (완료)
1. **BAM index 생성** - `rule index_bam` 추가
2. **QC summary 생성** - `rule generate_qc_summary` 추가
3. **Manifest 생성** - `rule generate_manifest` 추가
4. **표준 구조 BAM 복사** - `rule copy_bam_to_standard` 추가
5. **rule all 업데이트** - `get_all_targets()` 함수로 동적 타겟 관리

---

## 🧪 테스트 시나리오

### 테스트 1: Dry-run으로 규칙 확인

```bash
cd /home/ngs/ngs-pipeline/rna-seq-pipeline

# 표준 구조 모드로 dry-run
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --dry-run \
    --cores 1
```

**예상 결과**: 
- 각 샘플에 대해 `copy_bam_to_standard`, `index_bam`, `generate_qc_summary`, `generate_manifest` 규칙 실행 계획 표시
- 에러 없이 dependency graph 생성

---

### 테스트 2: DAG 시각화

```bash
# Dependency graph 생성
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --dag | dot -Tpng > dag_standard.png

# Legacy 구조와 비교
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=false \
    --dag | dot -Tpng > dag_legacy.png
```

**확인 사항**:
- 표준 구조: `star_align` → `copy_bam_to_standard` → `index_bam` → `generate_qc_summary` → `generate_manifest`
- Legacy 구조: `star_align` → (기존 규칙)

---

### 테스트 3: 단일 샘플로 전체 파이프라인 테스트

```bash
# 새 샘플 데이터로 테스트 (h_RNA_200_1)
# 주의: 기존 결과가 있으면 --forcerun 사용

# 1. 특정 샘플의 manifest만 생성
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --cores 12 \
    --forceall \
    -- /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/manifest.json
```

**예상 출력 구조**:
```
/home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/
├── final_outputs/
│   ├── bam/
│   │   ├── aligned.sorted.bam
│   │   └── aligned.sorted.bam.bai
│   ├── qc/
│   │   └── qc_summary.json
│   └── manifest.json
├── intermediate/
│   └── logs/
│       ├── copy_bam.log
│       ├── samtools_index.log
│       ├── qc_summary.log
│       ├── manifest.log
│       └── star_final.log
└── metadata/
```

**검증 명령**:
```bash
# 1. 디렉토리 구조 확인
tree -L 4 /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/

# 2. Manifest 검증
python3 scripts/generate_manifest.py \
    --validate /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/manifest.json

# 3. QC summary 확인
cat /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/qc/qc_summary.json | jq '.overall_status'

# 4. BAM index 확인
samtools idxstats /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/bam/aligned.sorted.bam | head
```

---

### 테스트 4: 전체 샘플 파이프라인 실행

```bash
# 모든 샘플에 대해 표준 구조로 실행
snakemake --configfile config_human_H2O2.yaml \
    --config use_standard_structure=true \
    --cores 12 \
    --keep-going \
    2>&1 | tee pipeline_standard.log
```

**실행 순서**:
1. FastQC (raw)
2. Cutadapt (trimming)
3. STAR alignment
4. **copy_bam_to_standard** (표준 구조로 복사)
5. **index_bam** (BAM index 생성)
6. featureCounts (전체 샘플)
7. **generate_qc_summary** (각 샘플)
8. **generate_manifest** (각 샘플)
9. MultiQC (프로젝트 전체)

---

### 테스트 5: 프로젝트 전체 요약 확인

```bash
# 모든 샘플의 QC 상태 요약
echo "=== Project-wide QC Summary ==="
for sample in h_RNA_Cont_1 h_RNA_Cont_2 h_RNA_Cont_3 h_RNA_100_1 h_RNA_100_2 h_RNA_100_3 \
              h_RNA_200_1 h_RNA_200_2 h_RNA_200_3 h_RNA_D1_1 h_RNA_D1_2 h_RNA_D1_3 \
              h_RNA_D3_1 h_RNA_D3_2 h_RNA_D3_3; do
    
    MANIFEST="/home/ngs/data/results/H2O2_human_2025/${sample}/rna-seq/final_outputs/manifest.json"
    
    if [ -f "$MANIFEST" ]; then
        STATUS=$(jq -r '.qc_metrics.overall_status' $MANIFEST)
        UNIQUE=$(jq -r '.qc_metrics.alignment.uniquely_mapped_pct // 0' $MANIFEST)
        ASSIGN=$(jq -r '.qc_metrics.quantification.assignment_rate // 0' $MANIFEST)
        
        printf "%-15s: %-6s (Unique: %5.1f%%, Assign: %5.1f%%)\n" \
            "$sample" "$STATUS" "$UNIQUE" "$ASSIGN"
    else
        echo "${sample}: MISSING"
    fi
done

# 프로젝트 디렉토리 구조 확인
tree -L 3 /home/ngs/data/results/H2O2_human_2025/
```

---

## 🔍 문제 해결 (Troubleshooting)

### 문제 1: `copy_bam_to_standard` 규칙 실행 안됨
```bash
# 원인: USE_STANDARD=False로 설정되어 있음
# 해결: config 파일 확인
grep "use_standard_structure" config_human_H2O2.yaml
# 출력: use_standard_structure: true 확인
```

### 문제 2: `generate_qc_summary` 실패
```bash
# 로그 확인
tail -50 /home/ngs/data/results/H2O2_human_2025/{SAMPLE}/rna-seq/intermediate/logs/qc_summary.log

# 원인 1: STAR log 파일 경로 오류
# 원인 2: featureCounts summary 파일 없음
# 해결: featurecounts_quant 규칙이 먼저 실행되었는지 확인
```

### 문제 3: `generate_manifest` 실패
```bash
# 로그 확인
tail -50 /home/ngs/data/results/H2O2_human_2025/{SAMPLE}/rna-seq/intermediate/logs/manifest.log

# 원인: BAM 파일 또는 QC summary 없음
# 해결: 선행 규칙들이 성공했는지 확인
ls -lh /home/ngs/data/results/H2O2_human_2025/{SAMPLE}/rna-seq/final_outputs/bam/
ls -lh /home/ngs/data/results/H2O2_human_2025/{SAMPLE}/rna-seq/final_outputs/qc/
```

### 문제 4: Lambda 함수 오류
```bash
# 에러: NameError: name 'wildcards' is not defined
# 원인: lambda 함수에서 wildcards 사용 오류
# 해결: Snakefile의 lambda 함수 문법 확인
```

---

## ✅ 검증 체크리스트

### Phase 3 완료 확인

- [ ] **Dry-run 성공**
  - [ ] 표준 구조 모드
  - [ ] Legacy 구조 모드
  - [ ] 규칙 dependency 올바름

- [ ] **단일 샘플 테스트**
  - [ ] `copy_bam_to_standard` 실행
  - [ ] `index_bam` 실행 (BAM.bai 생성)
  - [ ] `generate_qc_summary` 실행 (qc_summary.json 생성)
  - [ ] `generate_manifest` 실행 (manifest.json 생성)
  - [ ] Manifest 검증 통과

- [ ] **디렉토리 구조**
  - [ ] `final_outputs/bam/` 생성
  - [ ] `final_outputs/qc/` 생성
  - [ ] `intermediate/logs/` 생성
  - [ ] `manifest.json` 최상위 위치

- [ ] **전체 파이프라인 테스트**
  - [ ] 모든 샘플 성공적으로 처리
  - [ ] QC summary 모든 샘플 생성
  - [ ] Manifest 모든 샘플 생성
  - [ ] 프로젝트 요약 (counts, MultiQC) 생성

- [ ] **QC 메트릭 검증**
  - [ ] PASS 샘플 올바르게 판정
  - [ ] WARN 샘플 올바르게 판정
  - [ ] FAIL 샘플 올바르게 판정

---

## 📊 예상 실행 시간

새 샘플 15개 기준 (paired-end, ~60M reads/sample):

| 단계 | 시간 (cores=12) | 비고 |
|------|----------------|------|
| FastQC (raw) | ~30분 | 병렬 실행 |
| Cutadapt | ~45분 | 병렬 실행 |
| STAR alignment | ~2-3시간 | 메모리 집약적 |
| copy_bam_to_standard | ~30분 | I/O 병목 |
| index_bam | ~15분 | 병렬 실행 |
| featureCounts | ~20분 | 전체 샘플 한번에 |
| generate_qc_summary | ~5분 | 가벼운 작업 |
| generate_manifest | ~30분 | MD5 계산 (4GB BAM) |
| MultiQC | ~5분 | 전체 요약 |
| **총계** | **~5-6시간** | |

---

## 🚀 다음 단계

Phase 3 테스트 완료 후:

1. **Phase 4: FastQC/Cutadapt 경로 표준화**
   - intermediate/fastqc/
   - intermediate/trimmed/
   
2. **Phase 5: featureCounts 샘플별 분리**
   - final_outputs/counts/gene_counts.csv (샘플별)
   - project_summary/counts/counts_matrix_all.txt (전체)

3. **Phase 6: Metadata 자동 생성**
   - metadata/sample_info.yaml
   - metadata/pipeline_config.yaml

4. **WGS/ATAC-seq 파이프라인 확장**
   - 동일한 표준 구조 적용
   - 교차 검증

---

## 💡 팁

### Snakemake 명령어 단축키

```bash
# Dry-run (빠른 확인)
alias sn-dry='snakemake --configfile config_human_H2O2.yaml --dry-run --cores 1'

# 단일 샘플 테스트
alias sn-one='snakemake --configfile config_human_H2O2.yaml --cores 12 --config use_standard_structure=true'

# 전체 실행
alias sn-all='snakemake --configfile config_human_H2O2.yaml --cores 12 --config use_standard_structure=true --keep-going'

# 특정 규칙만 재실행
sn-one --forcerun generate_manifest -- /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/final_outputs/manifest.json
```

### 로그 모니터링

```bash
# 실시간 로그 추적
tail -f /home/ngs/data/results/H2O2_human_2025/h_RNA_200_1/rna-seq/intermediate/logs/*.log

# 에러만 필터링
grep -i "error\|fail\|critical" /home/ngs/data/results/H2O2_human_2025/*/rna-seq/intermediate/logs/*.log
```

### Manifest 일괄 검증

```bash
# 모든 샘플의 manifest 한번에 검증
for manifest in /home/ngs/data/results/H2O2_human_2025/*/rna-seq/final_outputs/manifest.json; do
    echo "Validating: $manifest"
    python3 scripts/generate_manifest.py --validate $manifest
done
```
