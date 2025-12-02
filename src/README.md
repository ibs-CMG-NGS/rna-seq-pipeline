# Source Scripts

이 폴더는 RNA-seq 파이프라인에서 사용하는 Python 스크립트들을 포함합니다.

## 📁 파일 설명

### 🔬 분석 및 QC

**`generate_qc_report.py`** - QC 리포트 생성기
- Snakemake 파이프라인에서 자동 실행
- 독립 실행도 가능: `python3 src/generate_qc_report.py`
- 모든 샘플의 QC 통계를 HTML 리포트로 생성
- 입력: cutadapt 로그, STAR 로그, featureCounts 결과
- 출력: `results/qc_report.html`

**`check_results.py`** - 결과 점검 스크립트
- 파이프라인 실행 후 결과 검증
- 터미널에 요약 통계 출력
- 사용법: `python3 src/check_results.py`

### 🔧 유틸리티

**`check_fastq.py`** - FASTQ 파일 형식 검증
- FASTQ 파일의 무결성 확인
- 시퀀스/품질 점수 길이 일치 검사
- 사용법: `python3 src/check_fastq.py`

**`find_read.py`** - 특정 read 검색
- 문제가 있는 read를 찾아서 검사
- FASTQ 파일에서 특정 read ID 검색
- 사용법: `python3 src/find_read.py`

**`fix_fastq.py`** - FASTQ 파일 수정
- 비표준 형식의 FASTQ 파일 수정
- 2줄로 나뉜 시퀀스/품질 점수를 1줄로 병합
- 사용법: `python3 src/fix_fastq.py`

## 🚀 사용 방법

### Snakemake 파이프라인에서 자동 실행
```bash
# generate_qc_report.py는 자동으로 실행됨
snakemake --cores 8
```

### 수동 실행

```bash
# QC 리포트 생성
python3 src/generate_qc_report.py

# 결과 점검
python3 src/check_results.py

# FASTQ 파일 검증
python3 src/check_fastq.py

# 특정 read 찾기
python3 src/find_read.py

# FASTQ 파일 수정 (필요시)
python3 src/fix_fastq.py
```

## 📝 스크립트 수정

스크립트를 수정한 후:
1. Snakemake가 자동으로 변경 감지
2. 필요한 규칙만 재실행
3. QC 리포트 강제 재생성: `snakemake --cores 1 results/qc_report.html --force`

## 🔍 문제 해결

### ImportError 발생 시
```bash
# Python 경로 확인
which python3

# 필요한 패키지 설치
pip install -r requirements.txt
```

### 스크립트 실행 권한
```bash
chmod +x src/*.py
```

## 📚 의존성

- Python 3.7+
- Standard library only (외부 패키지 불필요)
- 모든 스크립트는 독립적으로 실행 가능
