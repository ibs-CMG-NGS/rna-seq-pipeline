# RNA-seq Pipeline - Project Structure

```
rna-seq-pipeline/
├── 📄 Snakefile                    # 메인 파이프라인 워크플로우
├── ⚙️ config.yaml                  # 파이프라인 설정 파일
├── 📦 environment.yaml              # Conda 환경 설정
├── 📖 README.md                     # 프로젝트 소개
├── 📚 PIPELINE_GUIDE.md             # 상세 사용 가이드
├── 📋 samples.tsv                   # 샘플 메타데이터 (선택)
│
├── 🐍 src/                          # Python 스크립트
│   ├── README.md                    # 스크립트 설명서
│   ├── generate_qc_report.py        # ⭐ QC 리포트 생성 (Snakemake)
│   ├── check_results.py             # 결과 점검 유틸리티
│   ├── check_fastq.py               # FASTQ 검증
│   ├── find_read.py                 # Read 검색
│   └── fix_fastq.py                 # FASTQ 수정
│
├── 📂 data/                         # 원본 데이터
│   └── raw/                         # Raw FASTQ 파일
│       ├── {sample}_1.fastq.gz
│       ├── {sample}_2.fastq.gz
│       └── backup/                  # 원본 백업 (선택)
│           └── md5sum.txt           # MD5 체크섬
│
├── 🧬 genome/                       # Reference 파일
│   ├── genome.fa                    # Reference genome
│   ├── genes.gtf                    # Gene annotation
│   └── star_index/                  # STAR index
│       ├── Genome
│       ├── SA
│       ├── SAindex
│       └── ...
│
├── 📊 results/                      # 분석 결과
│   ├── qc_report.html              # ⭐ HTML QC 리포트
│   │
│   ├── trimmed/                    # Cutadapt 결과
│   │   ├── {sample}_1.fastq.gz
│   │   └── {sample}_2.fastq.gz
│   │
│   ├── aligned/                    # STAR 정렬 결과
│   │   └── {sample}/
│   │       ├── Aligned.sortedByCoord.out.bam
│   │       ├── Log.final.out
│   │       ├── Log.out
│   │       └── SJ.out.tab
│   │
│   └── counts/                     # featureCounts 결과
│       ├── counts_matrix.txt       # ⭐ Count matrix
│       └── counts_matrix.txt.summary
│
└── 📋 logs/                        # 실행 로그
    ├── qc_report.log               # QC 리포트 로그
    ├── featurecounts.log           # featureCounts 로그
    ├── cutadapt/                   # Cutadapt 로그
    │   └── {sample}.log
    ├── star/                       # STAR 로그
    │   └── {sample}.log
    └── fastqc/                     # FastQC 로그
        └── {sample}_{read}_raw.log
```

## 📌 주요 파일 설명

### 설정 파일
| 파일 | 설명 |
|------|------|
| `Snakefile` | 파이프라인 워크플로우 정의 (7개 규칙) |
| `config.yaml` | 모든 파라미터 중앙 관리 |
| `environment.yaml` | Conda 환경 및 의존성 |

### Python 스크립트 (`src/`)
| 파일 | 용도 | 실행 방식 |
|------|------|----------|
| `generate_qc_report.py` | QC 리포트 생성 | Snakemake 자동 / 수동 |
| `check_results.py` | 결과 검증 | 수동 |
| `check_fastq.py` | FASTQ 무결성 검사 | 수동 |
| `find_read.py` | 특정 read 검색 | 수동 |
| `fix_fastq.py` | FASTQ 형식 수정 | 수동 |

### 결과 파일
| 파일 | 설명 | 크기 (예시) |
|------|------|-------------|
| `qc_report.html` | 종합 QC 리포트 | 15 KB |
| `counts_matrix.txt` | 유전자별 count matrix | 25 MB |
| `*.bam` | 정렬된 reads | 2 GB/샘플 |
| `trimmed/*.fastq.gz` | Trimmed reads | 1.9 GB/샘플 |

## 🔄 워크플로우

```
Raw FASTQ
    ↓
[cutadapt] → trimmed FASTQ
    ↓
[STAR] → BAM files
    ↓
[featureCounts] → count matrix
    ↓
[generate_qc_report] → qc_report.html
```

## 📝 파일 명명 규칙

### 입력 파일
- `{sample}_1.fastq.gz` - Read 1
- `{sample}_2.fastq.gz` - Read 2

### 출력 파일
- `results/trimmed/{sample}_1.fastq.gz`
- `results/aligned/{sample}/Aligned.sortedByCoord.out.bam`
- `results/counts/counts_matrix.txt` (모든 샘플 통합)
- `results/qc_report.html` (모든 샘플 통합)

## 🎯 디렉토리 용도

| 디렉토리 | 용도 | 백업 필요? |
|---------|------|-----------|
| `data/raw/` | 원본 데이터 | ✅ 필수 |
| `genome/` | Reference 파일 | ✅ 권장 |
| `src/` | 스크립트 | ✅ Git |
| `results/` | 분석 결과 | ⚠️ 재생성 가능 |
| `logs/` | 실행 로그 | ❌ 불필요 |
| `.snakemake/` | Snakemake 메타데이터 | ❌ 불필요 |

## 💾 디스크 사용량 예상

**단일 샘플 (Ctrl_3 기준):**
- Raw FASTQ: 4.1 GB
- Trimmed FASTQ: 3.8 GB
- Aligned BAM: 2.0 GB
- **소계: ~10 GB/샘플**

**전체 프로젝트 (18개 샘플):**
- Raw data: ~74 GB
- Results: ~180 GB
- **총합: ~260 GB 권장**
