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

## 🤖 LLM Agent로 파이프라인 실행 (권장)

자연어 대화로 파이프라인의 모든 단계를 처리할 수 있습니다.  
명령어를 직접 입력하지 않아도 됩니다.

```bash
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline

python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen2.5:32b
```

**대화 예시 (처음부터 끝까지):**
```
You: /data_3tb/shared/chd8-rna-seq-raw-data/fastq/ 에서 FASTQ 파일 찾아줘
Agent: 38개 샘플 감지 (149GB, paired-end) ✅

You: 새 프로젝트 만들어줘. ID는 my-project, 결과는 /data_3tb/shared/output/my-project/
Agent: config/projects/my-project.yaml 생성 완료 ✅

You: 리소스 얼마나 필요해?
Agent: 출력 약 746GB, 소요 5-10시간, 디스크/RAM 충분 ✅

You: 샘플 시트 만들어줘
Agent: wildtype 20개, heterozygous 18개 자동 분류 → TSV 저장 ✅

You: 데이터 검증해줘
Agent: FASTQ ✅, 디스크 ✅, genome 경로 ⚠️ 확인 필요

You: dry-run 해줘
Agent: 270 jobs, 11 rules — 문제 없음 ✅

You: 16 cores로 실행해줘
Agent: 파이프라인 시작 🚀

You: 상태 보여줘
Agent: 진행률 45%, star_align 17/38 진행 중...
```

> **빠른 시작 가이드**: `docs/user/LLM_AGENT_QUICKSTART.md`

---

## 🔧 수동 실행 방법

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

# 네트워크 파일 시스템 사용 시 (NFS, 원격 스토리지)
snakemake --cores 18 --jobs 2 --resources mem_gb=60 --latency-wait 120

# 파일 시스템 latency가 심한 경우
snakemake --cores 16 --jobs 2 --resources mem_gb=60 --latency-wait 180
```

**파라미터 설명:**
- `--cores N`: 전체 사용할 최대 CPU 코어 수
- `--jobs N`: 동시에 실행할 최대 작업(샘플) 수
- `--resources mem_gb=N`: 전체 메모리 제한 (GB)
- `--use-conda`: Conda 환경 자동 활성화
- `--latency-wait N`: 출력 파일 대기 시간 (초) - 네트워크 파일 시스템용

**네트워크 파일 시스템 사용 시 주의사항:**
- 원격 스토리지(NFS, SMB 등)에 결과를 저장하는 경우
- `MissingOutputException` 에러 발생 시 `--latency-wait` 값을 증가
- 파일 시스템 동기화 지연으로 인해 파일이 즉시 보이지 않을 수 있음
- 권장값: 60-180초 (네트워크 상태에 따라 조정)

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
│   ├── counts_matrix.txt              # featureCounts 원본 출력 (메타데이터 포함)
│   ├── counts_matrix.txt.summary      # 매핑 통계
│   └── counts_matrix_clean.csv        # DE 분석용 클린 매트릭스 ⭐
└── qc_report.html        # QC 리포트 (자동 생성)

src/                      # Python 스크립트
├── generate_qc_report.py    # QC 리포트 생성기 (Snakemake에서 사용)
├── convert_counts_matrix.py # Count matrix 변환기 (DE 분석용) ⭐
├── check_results.py         # 결과 점검 유틸리티
├── check_fastq.py           # FASTQ 검증 유틸리티
├── find_read.py             # Read 검색 유틸리티
└── fix_fastq.py             # FASTQ 수정 유틸리티
```

### 🔄 Downstream 분석 연결

**DE (Differential Expression) 분석을 위한 파일:**
- `results/counts/counts_matrix_clean.csv` - DESeq2/edgeR/limma-voom용 정리된 count matrix
  - 유전자 ID를 행 이름(index)으로
  - 샘플 이름을 열 이름으로
  - 메타데이터 컬럼(Chr, Start, End 등) 제거됨
  - 샘플 경로가 깔끔한 이름으로 변환됨 (예: `GABA_8`, `Ctrl_3`)

**파일 형식 예시:**
```
Geneid,GABA_8,Ctrl_4,GABA_5,Ctrl_3,...
ENSMUSG00000104478,0,0,0,0,...
ENSMUSG00000104385,0,0,0,0,...
ENSMUSG00000086053,0,0,2,0,...
```

이 파일을 바로 `RNA-Seq_DE_GO_analysis` 파이프라인의 `config.yml`에서 `count_data_path`로 지정할 수 있습니다.

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

### MissingOutputException (파일 시스템 Latency 문제)

**증상:**
```
MissingOutputException: Job completed successfully, but some output files are missing.
Missing files after 5 seconds. This might be due to filesystem latency.
```

**원인:**
- 네트워크 파일 시스템(NFS, SMB 등) 사용 시 파일 동기화 지연
- 파일이 실제로 생성되었지만 메타데이터가 아직 전파되지 않음

**해결 방법:**
```bash
# 대기 시간을 60초로 증가
snakemake --cores 18 --jobs 2 --latency-wait 60

# 심한 경우 120-180초로 증가
snakemake --cores 16 --jobs 2 --latency-wait 180

# 동시 작업 수를 줄이는 것도 도움이 됨
snakemake --cores 12 --jobs 1 --latency-wait 120
```

### Command not found (cutadapt, fastqc, STAR 등)

**증상:**
```
/usr/bin/bash: line 2: cutadapt: command not found
exit status 127
```

**해결 방법:**
```bash
# Conda 환경이 활성화되어 있는지 확인
conda activate rna-seq-pipeline

# 또는 snakemake_env 등 사용 중인 환경 활성화
conda activate snakemake_env

# 필요한 도구가 설치되어 있는지 확인
conda list | grep -E "cutadapt|fastqc|star|subread"

# 도구가 없다면 설치
conda install -c bioconda cutadapt fastqc star subread samtools
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
