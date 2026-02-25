# RNA-seq Analysis Pipeline

Paired-end RNA-seq 데이터 분석을 위한 Snakemake 기반 자동화 파이프라인입니다.

## ✨ 주요 기능

- 📊 **자동화된 QC**: FastQC 결과 자동 평가 및 문제 샘플 식별
- 🔄 **프로젝트 수준 요약**: 모든 샘플의 QC 상태를 통합 관리
- 🤖 **LLM 기반 에이전트**: 자연어로 파이프라인 상태 조회 및 분석 (로컬 LLM 지원)
- 🔗 **파이프라인 연결**: RNA-seq → DE/GO 분석 자동 브리지
- 📁 **표준화된 출력**: 체계적인 메타데이터 및 결과 관리
- 🔒 **데이터 보안**: 로컬 LLM으로 민감한 genomic 데이터 보호

## 📋 파이프라인 개요

이 파이프라인은 다음 단계를 자동으로 수행합니다:

1. **Quality Control (FastQC)** - 원본 데이터 품질 검사
2. **🆕 Automatic QC Evaluation** - FastQC 결과 자동 평가 및 문제 샘플 식별
3. **Adapter Trimming (cutadapt)** - 어댑터 제거 및 품질 기반 트리밍
4. **Quality Control (FastQC)** - 트리밍 후 품질 검사
5. **Alignment (STAR)** - 레퍼런스 게놈에 리드 정렬
6. **Quantification (featureCounts)** - 유전자 발현량 정량화
7. **Summary Report (MultiQC)** - 전체 분석 품질 리포트 생성
8. **🆕 LLM Agent** - 자연어로 결과 조회 및 downstream 분석 시작

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
rna-seq-pipeline/
├── Snakefile              # 파이프라인 워크플로우 정의
├── environment.yaml       # Conda 환경 정의
├── config/                # 설정 파일
│   ├── default.yaml       # 기본 설정 템플릿
│   ├── projects/          # 프로젝트별 설정
│   │   └── H2O2_human_2025.yaml
│   └── samples/           # 샘플 시트
│       ├── master.csv     # 마스터 샘플 시트
│       └── template.tsv   # 샘플 시트 템플릿
├── src/                   # 파이프라인 소스 코드
│   ├── qc/                # QC 스크립트
│   ├── preprocessing/     # 전처리 스크립트
│   ├── quantification/    # 정량화 스크립트
│   └── utils/             # 유틸리티
├── scripts/               # 독립 실행 도구
│   ├── standardization/   # 표준화 도구
│   └── data/              # 데이터 관리 도구
├── docs/                  # 문서
│   ├── user/              # 사용자 가이드
│   └── developer/         # 개발자 문서
├── tests/                 # 테스트
├── data/
│   └── raw/              # 원본 FASTQ 파일 위치
│       ├── sample1_1.fastq.gz
│       ├── sample1_2.fastq.gz
│       ├── sample2_1.fastq.gz
│       └── sample2_2.fastq.gz
├── genome/
│   ├── star_index/       # STAR genome index
│   └── genes.gtf         # Gene annotation 파일
├── results/              # 분석 결과
│   ├── qc/               # FastQC 및 MultiQC 결과
│   ├── trimmed/          # 트리밍된 FASTQ 파일
│   ├── aligned/          # STAR 정렬 결과 (BAM 파일)
│   └── counts/           # featureCounts 결과
└── logs/                 # 각 작업의 로그 파일
    ├── fastqc/           # FastQC 로그
    ├── cutadapt/         # Cutadapt 로그
    └── star/             # STAR 로그
```
    ├── cutadapt/         # cutadapt 로그
    └── star/             # STAR 정렬 로그
```

## 🆕 새 프로젝트 시작하기

새로운 RNA-seq 프로젝트를 시작하는 경우 다음 단계를 따라주세요:

### Step 1: 프로젝트 설정 파일 생성

```bash
# 템플릿을 복사하여 새 프로젝트 설정 파일 생성
cp config/config.yaml config/projects/config_MY_PROJECT.yaml
```

> **참고:** `config_*.yaml` 파일은 Git에서 추적되지 않으므로 프로젝트별로 자유롭게 생성 가능합니다.

### Step 2: 프로젝트 설정 수정

`config/projects/config_MY_PROJECT.yaml` 파일을 열어 다음 항목들을 수정하세요:

#### 필수 수정 항목

```yaml
# 1. 데이터 경로
data_dir: "/home/ngs/data/MY_PROJECT/fastq"  # FASTQ 파일 위치 (절대경로 권장)
raw_data_subdir: ""  # FASTQ가 data_dir 바로 아래 있으면 비워둠

# 2. 프로젝트 정보
project_id: "MY_PROJECT_2025"  # 고유한 프로젝트 ID
pipeline_type: "rna-seq"  # 파이프라인 타입

# 3. 결과 저장 경로
base_results_dir: "/home/ngs/data/results"  # 결과 저장 베이스 디렉토리
use_standard_structure: true  # 표준 구조 사용 (권장)

# 4. Reference 파일 경로
genome_dir: "/home/ngs/data/genome/SPECIES_BUILD"  # 예: human_GRCh38, mouse_GRCm39
star_index: "/home/ngs/data/genome/SPECIES_BUILD/star_index/"
annotation_gtf: "/home/ngs/data/genome/SPECIES_BUILD/genes.gtf"
```

#### Species별 조정 항목

```yaml
# 5. FastQC 평가 기준 (species에 따라 GC content 조정 필요)
fastqc_evaluation:
  # Human (GRCh38): 35-65% (GC ~40-45%)
  # Mouse (GRCm39): 40-55% (GC ~42%)
  # Rat (Rnor_6.0): 40-55% (GC ~42%)
  min_gc_content: 35  # Species에 맞게 조정
  max_gc_content: 65  # Species에 맞게 조정
  min_total_sequences: 5000000  # 최소 5M reads
```

#### 선택 수정 항목

```yaml
# 6. 계산 리소스 (서버 사양에 맞게 조정)
star_threads: 12  # CPU 코어 수
star_memory_gb: 35  # 사용 가능 RAM
featurecounts_threads: 8
cutadapt_threads: 4

# 7. Strandedness (시퀀싱 프로토콜 확인 필요)
strandedness: 0  # 0=unstranded, 1=forward, 2=reverse
```

### Step 3: 샘플 정보 기록 (선택, 문서화용)

```bash
# 샘플 정보를 TSV 파일로 기록 (파이프라인 실행에는 불필요하지만 문서화에 유용)
cp config/samples/template.tsv config/samples/MY_PROJECT.tsv
# 에디터로 샘플 정보 작성
```

> **중요:** 파이프라인은 FASTQ 파일명에서 자동으로 샘플 리스트를 추출하므로, 
> 샘플 시트는 문서화 목적으로만 사용됩니다.

### Step 4: FASTQ 파일 준비

파이프라인은 다양한 FASTQ 파일명 패턴을 자동 인식합니다:

#### 지원되는 파일명 패턴

```bash
# Pattern 1 (기본)
{sample}_1.fastq.gz / {sample}_2.fastq.gz

# Pattern 2 (일반적)
{sample}_R1.fastq.gz / {sample}_R2.fastq.gz

# Pattern 3 (Illumina)
{sample}_R1_001.fastq.gz / {sample}_R2_001.fastq.gz

# Pattern 4
{sample}.1.fastq.gz / {sample}.2.fastq.gz

# Pattern 5 (단축형)
{sample}_1.fq.gz / {sample}_2.fq.gz
{sample}_R1.fq.gz / {sample}_R2.fq.gz
```

#### 파일명 예시

```bash
# ✅ 올바른 예시 (자동 인식됨)
Control_1_1.fastq.gz, Control_1_2.fastq.gz
Treatment_R1.fastq.gz, Treatment_R2.fastq.gz
Sample01_R1_001.fastq.gz, Sample01_R2_001.fastq.gz

# ❌ 잘못된 예시 (인식 안 됨)
sample.fa.gz  # .fastq.gz 또는 .fq.gz 사용 필요
sample_R1.fasta.gz  # .fastq.gz 또는 .fq.gz 사용 필요
sample_forward.fastq.gz  # _1, _2 또는 _R1, _R2 규칙 필요
```

> **중요:** 
> - 확장자는 반드시 `.fastq.gz` 또는 `.fq.gz`여야 합니다.
> - Read 1과 Read 2 파일은 숫자만 다르고 나머지는 동일해야 합니다.
> - 파이프라인은 첫 번째 매칭되는 패턴을 자동 감지하여 사용합니다.

#### 패턴 감지 확인

Dry-run 실행 시 감지된 패턴이 출력됩니다:

```bash
# 출력 예시:
#   Detected FASTQ Pattern: {sample}_R1.fastq.gz
#   Found 15 samples
```

### Step 5: 설정 검증 (Dry-run)

```bash
# 파이프라인이 올바르게 설정되었는지 확인
snakemake --configfile config/projects/config_MY_PROJECT.yaml \
  --config use_standard_structure=true \
  --dry-run --cores 1

# 출력 예시:
# ================================================================================
# PIPELINE CONFIGURATION:
#   Project ID: MY_PROJECT_2025
#   Found 15 samples
#   Sample list: ['sample1', 'sample2', ...]
# ================================================================================
# Job stats: ... (총 작업 수 표시)
```

### Step 6: 파이프라인 실행

```bash
# 실제 실행
snakemake --configfile config/projects/config_MY_PROJECT.yaml \
  --config use_standard_structure=true \
  --cores 8
```

---

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
- **🆕 `results/qc/fastqc_evaluation.txt`** - FastQC 자동 평가 리포트 (PASS/WARN/FAIL)
- **🆕 `results/qc/fastqc_evaluation.json`** - FastQC 평가 결과 (JSON 형식)
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

## 📖 추가 문서

### 사용자 가이드
- **[PIPELINE_GUIDE.md](docs/user/PIPELINE_GUIDE.md)** - 파이프라인 상세 사용법
- **[FASTQC_GUIDE.md](docs/user/FASTQC_GUIDE.md)** - FastQC 리포트 해석 상세 가이드
- **🆕 [FASTQC_AUTO_EVAL_GUIDE.md](docs/user/FASTQC_AUTO_EVAL_GUIDE.md)** - FastQC 자동 평가 기능 사용법
- **[QC_REPORT_GUIDE.md](docs/user/QC_REPORT_GUIDE.md)** - QC 리포트 해석 가이드
- **🆕 [OLLAMA_SETUP.md](docs/user/OLLAMA_SETUP.md)** - LLM 에이전트 설정 가이드 (로컬 LLM)
- **🆕 [LLM_AGENT_QUICKSTART.md](docs/user/LLM_AGENT_QUICKSTART.md)** - LLM 에이전트 빠른 시작 가이드

### 개발자 문서
- **[PROJECT_STRUCTURE.md](docs/developer/PROJECT_STRUCTURE.md)** - 프로젝트 구조 설명
- **🆕 [PHASE6_LLM_INTEGRATION.md](docs/developer/PHASE6_LLM_INTEGRATION.md)** - LLM 통합 아키텍처 및 구현
- **[STANDARDIZATION.md](docs/developer/STANDARDIZATION.md)** - 파이프라인 표준화 가이드
- **[PROJECT_REORGANIZATION.md](docs/developer/PROJECT_REORGANIZATION.md)** - 프로젝트 재구성 계획

## 📧 문의

문제가 발생하거나 질문이 있으시면 이슈를 등록해주세요.
