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

### 5. QC 리포트만 재생성

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
