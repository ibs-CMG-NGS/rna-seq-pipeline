# RNA-seq Analysis Pipeline

Snakemake 기반 RNA-seq 데이터 분석 파이프라인

## 📋 파이프라인 단계

1. **Quality Control (FastQC)** - Raw 데이터 품질 확인
2. **Adapter Trimming (Cutadapt)** - 어댑터 제거 및 품질 필터링
3. **Alignment (STAR)** - Reference genome에 정렬
4. **Quantification (featureCounts)** - 유전자별 read count 계산
5. **QC Report** - 전체 분석 결과 HTML 리포트 생성

## 🔧 설정 파일 (config.yaml)

모든 파라미터는 `config.yaml`에서 관리됩니다:

```yaml
# Reference 파일
star_index: "genome/star_index/"
annotation_gtf: "genome/genes.gtf"

# Adapter 시퀀스
adapter_r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Quality control
quality_cutoff: 20
min_read_length: 20

# Threading
star_threads: 8
featurecounts_threads: 4

# QC Report
generate_qc_report: true
qc_report_output: "results/qc_report.html"
qc_top_genes: 10
```

## 🚀 실행 방법

### 1. 환경 설정

```bash
# Conda 환경 생성
conda env create -f environment.yaml
conda activate rna-seq-pipeline
```

### 2. 데이터 준비

```bash
# Raw FASTQ 파일을 data/raw/에 배치
# 파일명 형식: {sample}_1.fastq.gz, {sample}_2.fastq.gz
# 예: Ctrl_3_1.fastq.gz, Ctrl_3_2.fastq.gz
```

### 3. Dry-run (테스트)

```bash
# 실제 실행 없이 계획만 확인
snakemake -n
```

### 4. 파이프라인 실행

```bash
# 모든 코어 사용
snakemake --cores all

# 특정 코어 수 지정
snakemake --cores 8

# 특정 단계까지만 실행
snakemake --cores 8 results/counts/counts_matrix.txt
```

### 5. 성능 최적화 (서버 사양별 설정)

#### 5.1 리소스 분석

파이프라인의 각 단계는 서로 다른 컴퓨팅 리소스를 요구합니다:

| 단계 | CPU 사용량 | 메모리 사용량 | 특징 |
|------|-----------|---------------|------|
| FastQC | 낮음 (1 thread) | ~500MB | I/O 집약적 |
| cutadapt | 중간 (4 threads) | ~1-2GB | CPU 집약적 |
| **STAR align** | **높음 (8-16 threads)** | **~30GB/샘플** | **메모리 집약적** |
| featureCounts | 중간 (4-8 threads) | ~2-3GB | I/O 집약적 |

#### 5.2 서버 사양별 추천 설정

**예시: Intel Xeon E5-2630 v2 (12 cores/24 threads, 62GB RAM)**

`config.yaml` 설정:
```yaml
star_threads: 12              # CPU 코어 수에 맞춤
featurecounts_threads: 8      # 비교적 가벼운 작업이므로 여유있게
```

`Snakefile`의 STAR 규칙 (이미 적용됨):
```python
rule star_align:
    resources:
        mem_gb=35             # 샘플당 메모리 제한 (동시 실행 고려)
    shell:
        """
        STAR --limitBAMsortRAM 30000000000  # 30GB (메모리 활용 최적화)
        """
```

#### 5.3 실행 명령어 (리소스 제어)

```bash
# 추천: 안전한 설정 (메모리 부족 방지)
snakemake --cores 16 --jobs 2 --use-conda

# 샘플이 많을 때: 3개 샘플 동시 처리
snakemake --cores 18 --jobs 3 --resources mem_gb=60 --use-conda

# 보수적인 설정: 한 번에 하나씩 (메모리 부족 시)
snakemake --cores 12 --jobs 1 --use-conda
```

**파라미터 설명:**
- `--cores N`: 전체 사용할 최대 CPU 코어 수
- `--jobs N`: 동시에 실행할 최대 작업(샘플) 수
- `--resources mem_gb=N`: 전체 메모리 제한 (GB)
- `--use-conda`: Conda 환경 자동 활성화

#### 5.4 리소스 모니터링

실행 중 시스템 리소스를 모니터링하세요:

```bash
# 실시간 CPU/메모리 모니터링
htop

# 메모리 사용량 확인
watch -n 1 free -h

# STAR 프로세스 모니터링
watch -n 1 'ps aux | grep STAR | head -5'

# 디스크 I/O 모니터링
iostat -x 2
```

#### 5.5 일반적인 서버 사양별 가이드

| 서버 사양 | star_threads | featurecounts_threads | snakemake --jobs | 비고 |
|-----------|--------------|----------------------|------------------|------|
| 8 cores, 32GB RAM | 6 | 4 | 1 | 메모리 제약 큼 |
| 12 cores, 62GB RAM | 12 | 8 | 2 | **현재 서버** |
| 24 cores, 128GB RAM | 16 | 12 | 4 | 여유로운 처리 |
| 48+ cores, 256GB+ RAM | 24 | 16 | 8 | 대규모 분석 |

**주의사항:**
- STAR 정렬은 메모리를 가장 많이 사용하므로 `--jobs` 값을 신중히 설정
- 메모리 부족(OOM) 에러 발생 시 `--jobs` 값을 줄이거나 `star_threads` 감소
- 디스크 I/O가 병목이 될 수 있으므로 SSD 사용 권장

### 6. QC 리포트만 재생성

```bash
snakemake --cores 1 results/qc_report.html --force
```

## 📊 결과 파일

```
results/
├── trimmed/              # Adapter 제거된 FASTQ
│   ├── {sample}_1.fastq.gz
│   └── {sample}_2.fastq.gz
├── aligned/              # STAR 정렬 결과
│   └── {sample}/
│       ├── Aligned.sortedByCoord.out.bam
│       └── Log.final.out
├── counts/               # featureCounts 결과
│   ├── counts_matrix.txt
│   └── counts_matrix.txt.summary
└── qc_report.html        # QC 리포트 (자동 생성)

src/                      # Python 스크립트
├── generate_qc_report.py # QC 리포트 생성기 (Snakemake에서 사용)
├── check_results.py      # 결과 점검 유틸리티
├── check_fastq.py        # FASTQ 검증 유틸리티
├── find_read.py          # Read 검색 유틸리티
└── fix_fastq.py          # FASTQ 수정 유틸리티
```

## 📈 QC 리포트 내용

HTML 리포트에 포함된 정보:
- ✂️ Adapter trimming 통계 (샘플별)
- 🎯 Alignment 품질 (매핑률, 진행 바)
- 🧮 Gene quantification 통계
- 📊 유전자 발현 통계 (검출된 유전자 수)
- 🔝 고발현 유전자 Top N (설정 가능)
- 💾 파일 크기 정보

## ⚙️ 설정 커스터마이징

### QC 리포트 설정

```yaml
# config.yaml에서 수정
generate_qc_report: true    # false로 설정하면 리포트 생성 안 함
qc_report_output: "results/custom_report.html"  # 출력 경로 변경
qc_top_genes: 20           # 표시할 고발현 유전자 수
```

### Threading 조정

```yaml
star_threads: 4            # 서버 사양에 맞게 조정
featurecounts_threads: 2
```

### Quality 파라미터

```yaml
quality_cutoff: 30         # 더 엄격한 품질 필터링
min_read_length: 50        # 더 긴 최소 길이
```

## 🔍 문제 해결

### 파이프라인 실패 시

```bash
# 로그 확인
cat logs/cutadapt/{sample}.log
cat logs/star/{sample}.log
cat logs/featurecounts.log
cat logs/qc_report.log

# 특정 단계 다시 실행
snakemake --cores 8 --forcerun star_align
```

### QC 리포트가 생성되지 않을 때

```bash
# 수동으로 리포트 생성
python3 src/generate_qc_report.py

# 또는 강제 재생성
snakemake --cores 1 results/qc_report.html --force

# 결과 점검
python3 src/check_results.py
```

### FASTQ 파일 문제

```bash
# FASTQ 파일 무결성 검증
python3 src/check_fastq.py

# 비표준 형식 수정 (필요시)
python3 src/fix_fastq.py
```

## 📝 샘플 추가

1. 새 FASTQ 파일을 `data/raw/`에 추가
2. 파이프라인 재실행 (Snakemake가 자동으로 새 샘플 감지)
3. QC 리포트가 자동으로 모든 샘플 포함

## 🎯 Next Steps

분석 완료 후:
1. `results/qc_report.html` 브라우저에서 열기
2. `results/counts/counts_matrix.txt` → DESeq2/edgeR로 differential expression 분석
3. `samples.tsv` 파일에 샘플 메타데이터 정리

## 📚 참고 자료

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [STAR Manual](https://github.com/alexdobin/STAR)
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [featureCounts](http://subread.sourceforge.net/)
