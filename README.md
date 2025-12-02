# RNA-seq Analysis Pipeline

Paired-end RNA-seq 데이터 분석을 위한 Snakemake 기반 자동화 파이프라인입니다.

## 📋 파이프라인 개요

이 파이프라인은 다음 단계를 자동으로 수행합니다:

1. **Quality Control (FastQC)** - 원본 데이터 품질 검사
2. **Adapter Trimming (cutadapt)** - 어댑터 제거 및 품질 기반 트리밍
3. **Quality Control (FastQC)** - 트리밍 후 품질 검사
4. **Alignment (STAR)** - 레퍼런스 게놈에 리드 정렬
5. **Quantification (featureCounts)** - 유전자 발현량 정량화
6. **Summary Report (MultiQC)** - 전체 분석 품질 리포트 생성

## 🔧 요구사항

### 소프트웨어
- Conda 또는 Mamba
- Linux 환경 (WSL 포함)

### 데이터
- Paired-end FASTQ 파일 (`.fastq.gz` 형식)
- STAR genome index
- Gene annotation GTF 파일

## 📁 디렉토리 구조

```
rna_seq_pipeline/
├── Snakefile              # 파이프라인 워크플로우 정의
├── config.yaml            # 설정 파일 (사용자 수정 필요)
├── environment.yaml       # Conda 환경 정의
├── data/
│   └── raw/              # 원본 FASTQ 파일 위치
│       ├── sample1_R1.fastq.gz
│       ├── sample1_R2.fastq.gz
│       ├── sample2_R1.fastq.gz
│       └── sample2_R2.fastq.gz
├── genome/
│   ├── star_index/       # STAR genome index
│   └── annotation.gtf    # Gene annotation 파일
├── results/
│   ├── qc/               # FastQC 및 MultiQC 결과
│   ├── trimmed/          # 트리밍된 FASTQ 파일
│   ├── aligned/          # STAR 정렬 결과 (BAM 파일)
│   └── counts/           # featureCounts 결과
└── logs/                 # 각 작업의 로그 파일
    ├── fastqc/           # FastQC 로그
    ├── cutadapt/         # cutadapt 로그
    └── star/             # STAR 정렬 로그
```

## 🚀 사용 방법

### 1. Conda 환경 설정

```bash
# Conda 환경 생성
conda env create -f environment.yaml

# 환경 활성화
conda activate rna-seq-pipeline
```

### 2. 데이터 준비

#### 방법 1: 자동 다운로드 (권장)

`data/raw/md5sum.txt` 파일에 다운로드 링크가 포함되어 있는 경우, 자동 다운로드 스크립트를 사용할 수 있습니다:

```bash
# Python 스크립트 사용 (권장)
python download_fastq.py

# 또는 Bash 스크립트 사용
chmod +x download_fastq.sh
./download_fastq.sh
```

**주요 기능:**
- 모든 FASTQ 파일 자동 다운로드
- MD5 checksum을 통한 파일 무결성 자동 검증
- 이미 다운로드된 파일은 자동으로 스킵
- 손상된 파일 자동 재다운로드

#### 방법 2: 수동 복사

```bash
# FASTQ 파일을 data/raw/ 디렉토리에 복사
# 파일명 규칙: {sample_name}_R1.fastq.gz, {sample_name}_R2.fastq.gz
cp /path/to/your/fastq/*_R*.fastq.gz data/raw/
```

### 3. 설정 파일 생성

#### 프로젝트별 설정 파일 만들기

`config.yaml`은 템플릿 파일이므로 직접 수정하지 마세요. 대신 프로젝트별로 복사하여 사용합니다:

```bash
# 템플릿을 복사하여 프로젝트별 설정 파일 생성
cp config.yaml config_my_project.yaml

# 설정 파일 수정
nano config_my_project.yaml  # 또는 vi, code 등 원하는 에디터 사용
```

#### 주요 설정 항목

```yaml
# === Directory Structure ===
# 대용량 스토리지를 사용하는 경우 절대 경로 지정
data_dir: "/home/ngs/data/rna-seq-pipeline/data/my_project"
results_dir: "/home/ngs/data/rna-seq-pipeline/results/my_project"
logs_dir: "/home/ngs/data/rna-seq-pipeline/logs/my_project"

# === Reference Files ===
star_index: "genome/star_index/"
annotation_gtf: "genome/genes.gtf"

# === Computational Resources ===
star_threads: 12           # 시스템 CPU 코어 수에 맞게 조정
star_memory_gb: 35         # 사용 가능한 RAM에 맞게 조정
featurecounts_threads: 8
cutadapt_threads: 4
```

**참고:** 
- `config_*.yaml` 파일은 Git에서 추적되지 않습니다 (`.gitignore`에 등록됨)
- 각 프로젝트/데이터셋마다 별도의 설정 파일을 만들어 관리하세요

### 4. 파이프라인 실행

```bash
# Dry-run (실제 실행하지 않고 작업 계획만 확인)
snakemake --configfile config_my_project.yaml -n

# 실제 실행 (설정 파일의 스레드 수만큼 자동 사용)
snakemake --configfile config_my_project.yaml -j 12

# 특정 결과물만 생성
snakemake --configfile config_my_project.yaml -j 12 results/qc_report.html
```

### 5. 워크플로우 시각화 (선택사항)

```bash
# DAG (Directed Acyclic Graph) 생성
snakemake --dag | dot -Tpdf > workflow.pdf

# Rulegraph 생성
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
```

## 📊 결과물

### 주요 출력 파일

- `results/qc/multiqc_report.html` - 전체 분석 품질 요약 리포트
- `results/counts/counts_matrix.txt` - 유전자별 read count 매트릭스
- `results/counts/counts_matrix.txt.summary` - featureCounts 통계

### 샘플별 출력 파일

- `results/qc/{sample}_{R1,R2}_*_fastqc.html` - FastQC 품질 리포트
- `results/trimmed/{sample}_{R1,R2}.fastq.gz` - 트리밍된 FASTQ 파일
- `results/aligned/{sample}/Aligned.sortedByCoord.out.bam` - 정렬된 BAM 파일
- `results/aligned/{sample}/Log.final.out` - STAR 정렬 통계

## ⚙️ 파라미터 설정

`config.yaml`에서 다음 파라미터를 조정할 수 있습니다:

### Quality Control & Trimming
- `quality_cutoff`: 최소 base quality (기본값: 20)
- `min_read_length`: 트리밍 후 최소 read 길이 (기본값: 20)
- `adapter_r1`, `adapter_r2`: Illumina adapter 서열

### Alignment
- `star_threads`: STAR 정렬에 사용할 스레드 수 (기본값: 8)

### Quantification
- `featurecounts_threads`: featureCounts에 사용할 스레드 수 (기본값: 4)
- `strandedness`: RNA-seq library의 strand 정보
  - `0`: unstranded (기본값)
  - `1`: stranded (forward)
  - `2`: stranded (reverse)
- `feature_type`: 정량화할 feature 타입 (기본값: "exon")
- `attribute_type`: GTF attribute 타입 (기본값: "gene_id")

## 🔍 문제 해결

### FASTQ 파일이 인식되지 않는 경우
- 파일명이 `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz` 형식인지 확인
- 파일이 `data/raw/` 디렉토리에 있는지 확인

### STAR index 오류
- `config.yaml`의 `star_index` 경로가 올바른지 확인
- STAR index가 사용하는 STAR 버전과 호환되는지 확인

### Annotation GTF 오류
- GTF 파일 경로가 올바른지 확인
- GTF 파일 형식이 유효한지 확인 (gene_id attribute 포함 여부)

## 📝 주의사항

1. **STAR genome index 생성**: 이 파이프라인은 STAR index가 이미 준비되어 있다고 가정합니다. Index 생성이 필요한 경우:
   ```bash
   STAR --runMode genomeGenerate \
        --genomeDir genome/star_index/ \
        --genomeFastaFiles genome/reference.fa \
        --sjdbGTFfile genome/annotation.gtf \
        --runThreadN 8
   ```

2. **메모리 요구사항**: STAR 정렬은 대용량 메모리가 필요합니다 (인간 게놈 기준 최소 32GB RAM 권장)

3. **Strandedness 확인**: RNA-seq library preparation 방법에 따라 `strandedness` 파라미터를 올바르게 설정해야 합니다.

## 📚 참고 문헌

- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- cutadapt: https://cutadapt.readthedocs.io/
- STAR: https://github.com/alexdobin/STAR
- featureCounts: http://subread.sourceforge.net/
- MultiQC: https://multiqc.info/
- Snakemake: https://snakemake.readthedocs.io/

## 📧 문의

문제가 발생하거나 질문이 있으시면 이슈를 등록해주세요.
