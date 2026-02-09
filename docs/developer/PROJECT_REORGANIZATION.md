# Project Reorganization Plan

## 🎯 목표
파이프라인 코드베이스를 체계적으로 정리하여 유지보수성과 확장성 향상

---

## 📂 현재 구조 분석

### 현재 디렉토리 구조
```
rna-seq-pipeline/
├── src/                          # Snakemake 규칙에서 사용하는 스크립트
│   ├── check_fastq.py
│   ├── check_results.py
│   ├── check_samples.py
│   ├── convert_counts_matrix.py
│   ├── evaluate_fastqc.py
│   ├── find_read.py
│   ├── fix_fastq.py
│   ├── generate_qc_report.py
│   └── README.md
├── scripts/                      # Phase 2/3에서 추가된 표준화 스크립트
│   ├── convert_sample_sheet.py
│   ├── generate_manifest.py
│   └── generate_qc_summary.py
├── docs/                         # Phase 2/3 문서만 있음
│   ├── PHASE2_TESTING.md
│   ├── PHASE3_PLAN.md
│   ├── PHASE3_TESTING.md
│   └── STANDARDIZATION_PROGRESS.md
├── config.yaml                   # 기본 config (루트)
├── config_human_H2O2.yaml        # 프로젝트별 config (루트)
├── config_H2O2_astrocyte.yaml    # 프로젝트별 config (루트)
├── samples.tsv                   # 샘플 시트 (루트)
├── samples_human_H2O2.tsv        # 샘플 시트 (루트)
├── samples_master.csv            # 마스터 샘플 시트 (루트)
├── download_fastq.py             # 유틸리티 (루트)
├── download_fastq.sh             # 유틸리티 (루트)
├── test_evaluate_fastqc.py       # 테스트 (루트)
├── README.md                     # 메인 문서 (루트)
├── PIPELINE_GUIDE.md             # 사용자 가이드 (루트)
├── FASTQC_GUIDE.md               # 사용자 가이드 (루트)
├── FASTQC_AUTO_EVAL_GUIDE.md     # 사용자 가이드 (루트)
├── QC_REPORT_GUIDE.md            # 사용자 가이드 (루트)
├── PROJECT_STRUCTURE.md          # 개발 문서 (루트)
├── STANDARDIZATION.md            # 개발 문서 (루트)
└── Snakefile                     # 파이프라인 정의
```

---

## 🆚 src/ vs scripts/ 차이점

### src/ (Source Code)
- **목적**: Snakemake 파이프라인의 **내부 로직**
- **사용**: Snakefile에서 `script:` directive로 직접 호출
- **특징**: 
  - 파이프라인 실행에 필수적
  - Snakemake context 접근 가능 (snakemake.input, snakemake.output 등)
  - 주로 데이터 처리/변환 로직
  
**현재 파일**:
- `evaluate_fastqc.py` - FastQC 결과 자동 평가
- `generate_qc_report.py` - HTML QC 리포트 생성
- `convert_counts_matrix.py` - Count matrix 형식 변환
- `check_*.py`, `fix_*.py` - 데이터 검증/수정 유틸리티

### scripts/ (Standalone Scripts)
- **목적**: **독립적으로 실행 가능한** 도구
- **사용**: 커맨드라인에서 직접 실행 (`python3 scripts/xxx.py`)
- **특징**:
  - CLI 인터페이스 (argparse 등)
  - 파이프라인 외부에서도 사용 가능
  - 표준화 프레임워크 도구
  
**현재 파일**:
- `generate_manifest.py` - Manifest.json 생성 (Phase 2)
- `generate_qc_summary.py` - QC summary JSON 생성 (Phase 2)
- `convert_sample_sheet.py` - 샘플 시트 형식 변환 (Phase 1)

---

## ✅ 개선된 구조 제안

```
rna-seq-pipeline/
├── config/                       # ✨ 새 폴더: 모든 설정 파일
│   ├── default.yaml              # (구 config.yaml)
│   ├── projects/
│   │   ├── H2O2_human_2025.yaml  # (구 config_human_H2O2.yaml)
│   │   └── H2O2_astrocyte.yaml   # (구 config_H2O2_astrocyte.yaml)
│   └── samples/
│       ├── master.csv            # (구 samples_master.csv)
│       ├── H2O2_human.tsv        # (구 samples_human_H2O2.tsv)
│       └── template.tsv          # (구 samples.template.tsv)
│
├── src/                          # Snakemake 내부 스크립트
│   ├── qc/
│   │   ├── evaluate_fastqc.py
│   │   └── generate_qc_report.py
│   ├── preprocessing/
│   │   ├── check_fastq.py
│   │   └── fix_fastq.py
│   ├── quantification/
│   │   └── convert_counts_matrix.py
│   └── utils/
│       ├── check_results.py
│       ├── check_samples.py
│       └── find_read.py
│
├── scripts/                      # 독립 실행 도구
│   ├── standardization/
│   │   ├── generate_manifest.py
│   │   ├── generate_qc_summary.py
│   │   └── convert_sample_sheet.py
│   └── data/
│       ├── download_fastq.py     # (루트에서 이동)
│       └── download_fastq.sh     # (루트에서 이동)
│
├── docs/                         # ✨ 모든 문서 통합
│   ├── user/                     # 사용자 가이드
│   │   ├── PIPELINE_GUIDE.md
│   │   ├── FASTQC_GUIDE.md
│   │   ├── FASTQC_AUTO_EVAL_GUIDE.md
│   │   └── QC_REPORT_GUIDE.md
│   ├── developer/                # 개발자 문서
│   │   ├── PROJECT_STRUCTURE.md
│   │   ├── STANDARDIZATION.md
│   │   ├── STANDARDIZATION_PROGRESS.md
│   │   ├── PHASE2_TESTING.md
│   │   ├── PHASE3_PLAN.md
│   │   └── PHASE3_TESTING.md
│   └── README.md                 # Docs 디렉토리 인덱스
│
├── tests/                        # ✨ 새 폴더: 모든 테스트
│   ├── test_evaluate_fastqc.py   # (루트에서 이동)
│   ├── test_manifest.py          # 추가 예정
│   └── test_qc_summary.py        # 추가 예정
│
├── Snakefile                     # 파이프라인 메인 정의
├── environment.yaml              # Conda 환경
├── .gitignore
└── README.md                     # 프로젝트 메인 README
```

---

## 📋 이동 계획

### Phase 1: Config 파일 통합
```bash
mkdir -p config/projects config/samples

# Config 파일 이동
mv config.yaml config/default.yaml
mv config_human_H2O2.yaml config/projects/H2O2_human_2025.yaml
mv config_H2O2_astrocyte.yaml config/projects/H2O2_astrocyte.yaml

# 샘플 시트 이동
mv samples_master.csv config/samples/master.csv
mv samples_human_H2O2.tsv config/samples/H2O2_human.tsv
mv samples.template.tsv config/samples/template.tsv
rm samples.tsv samples_converted.tsv  # 임시 파일 제거
```

### Phase 2: src/ 재구성
```bash
mkdir -p src/qc src/preprocessing src/quantification src/utils

# QC 관련
mv src/evaluate_fastqc.py src/qc/
mv src/generate_qc_report.py src/qc/

# 전처리
mv src/check_fastq.py src/preprocessing/
mv src/fix_fastq.py src/preprocessing/

# 정량화
mv src/convert_counts_matrix.py src/quantification/

# 유틸리티
mv src/check_results.py src/utils/
mv src/check_samples.py src/utils/
mv src/find_read.py src/utils/
```

### Phase 3: scripts/ 재구성
```bash
mkdir -p scripts/standardization scripts/data

# 표준화 스크립트
mv scripts/generate_manifest.py scripts/standardization/
mv scripts/generate_qc_summary.py scripts/standardization/
mv scripts/convert_sample_sheet.py scripts/standardization/

# 데이터 다운로드 스크립트
mv download_fastq.py scripts/data/
mv download_fastq.sh scripts/data/
```

### Phase 4: docs/ 통합
```bash
mkdir -p docs/user docs/developer

# 사용자 가이드
mv PIPELINE_GUIDE.md docs/user/
mv FASTQC_GUIDE.md docs/user/
mv FASTQC_AUTO_EVAL_GUIDE.md docs/user/
mv QC_REPORT_GUIDE.md docs/user/

# 개발자 문서
mv PROJECT_STRUCTURE.md docs/developer/
mv STANDARDIZATION.md docs/developer/
mv docs/STANDARDIZATION_PROGRESS.md docs/developer/
mv docs/PHASE2_TESTING.md docs/developer/
mv docs/PHASE3_PLAN.md docs/developer/
mv docs/PHASE3_TESTING.md docs/developer/
```

### Phase 5: tests/ 생성
```bash
mkdir -p tests

mv test_evaluate_fastqc.py tests/
```

---

## 🔧 파일 경로 업데이트 필요 목록

### 1. Snakefile
```python
# Before
configfile: "config.yaml"
script: "src/generate_qc_report.py"
shell: "python3 scripts/generate_qc_summary.py ..."

# After
configfile: "config/default.yaml"
script: "src/qc/generate_qc_report.py"
shell: "python3 scripts/standardization/generate_qc_summary.py ..."
```

### 2. Config 파일들
```yaml
# Before (config_human_H2O2.yaml)
star_index: "/home/ngs/data/genome/..."

# After (config/projects/H2O2_human_2025.yaml)
# 경로는 동일하지만, 파일 위치만 변경
```

### 3. Scripts (import 경로)
```python
# scripts/standardization/generate_manifest.py
# 상대 경로 import가 있다면 업데이트 필요
```

### 4. 문서 내 경로 참조
```markdown
# Before
See [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)

# After
See [PIPELINE_GUIDE.md](docs/user/PIPELINE_GUIDE.md)
```

### 5. README.md
```markdown
# Before
- [Pipeline Guide](PIPELINE_GUIDE.md)

# After
- [Pipeline Guide](docs/user/PIPELINE_GUIDE.md)
```

---

## 📊 이점

### 1. **명확한 책임 분리**
- `config/` - 설정만
- `src/` - 파이프라인 로직만
- `scripts/` - 독립 도구만
- `docs/` - 문서만
- `tests/` - 테스트만

### 2. **확장성 향상**
- 새 프로젝트: `config/projects/새프로젝트.yaml` 추가
- 새 도구: `scripts/카테고리/새도구.py` 추가
- 새 문서: `docs/user/` 또는 `docs/developer/` 추가

### 3. **찾기 쉬움**
- "Config는?" → `config/`
- "사용법은?" → `docs/user/`
- "개발 가이드는?" → `docs/developer/`
- "표준화 도구는?" → `scripts/standardization/`

### 4. **유지보수 용이**
- 관련 파일이 같은 폴더에 모여있음
- 명확한 파일 명명 규칙
- 중복 파일 제거 (samples.tsv, samples_converted.tsv)

---

## ⚠️ 주의사항

### 1. Git 이력 보존
```bash
# 파일 이동 시 git mv 사용 (이력 보존)
git mv config.yaml config/default.yaml

# 일반 mv 사용 시 이력 손실
mv config.yaml config/default.yaml  # ❌ 피하기
```

### 2. 하위 호환성
```bash
# 기존 실행 명령이 깨지지 않도록
# Old: snakemake --configfile config_human_H2O2.yaml
# New: snakemake --configfile config/projects/H2O2_human_2025.yaml

# 또는 루트에 symlink 생성
ln -s config/projects/H2O2_human_2025.yaml config_human_H2O2.yaml
```

### 3. CI/CD 업데이트
- GitHub Actions workflow 경로 수정 필요
- 테스트 스크립트 경로 수정

---

## 🚀 실행 순서

1. **백업 생성**
   ```bash
   git checkout -b refactor/project-reorganization
   ```

2. **Phase 1-5 순차적 실행**
   - 각 Phase마다 커밋
   - 각 Phase마다 테스트

3. **경로 업데이트**
   - Snakefile 수정
   - README.md 수정
   - 문서 내 링크 수정

4. **테스트**
   ```bash
   # Dry-run 확인
   snakemake --configfile config/projects/H2O2_human_2025.yaml --dry-run
   
   # 단위 테스트
   python3 -m pytest tests/
   ```

5. **문서 업데이트**
   - README.md에 새 구조 설명
   - docs/README.md 인덱스 생성

6. **Pull Request**
   - 변경사항 리뷰
   - 팀원 피드백 반영

---

## 💡 추가 제안

### 1. src/ 내부에 __init__.py 추가
```python
# src/qc/__init__.py
from .evaluate_fastqc import evaluate_fastqc
from .generate_qc_report import generate_report
```

### 2. scripts/ 도구에 --help 표준화
```python
# 모든 스크립트에 argparse 사용
parser = argparse.ArgumentParser(
    description="Generate manifest.json for pipeline outputs",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
```

### 3. VERSION 파일 추가
```
# VERSION
1.0.0-alpha
```

### 4. CHANGELOG.md 추가
```markdown
# Changelog

## [Unreleased]
### Added
- Project reorganization
- Standardization framework (Phase 1-3)

## [0.9.0] - 2025-01-27
### Added
- FastQC auto-evaluation
```

---

## 📝 체크리스트

- [ ] Phase 1: Config 파일 이동
- [ ] Phase 2: src/ 재구성
- [ ] Phase 3: scripts/ 재구성
- [ ] Phase 4: docs/ 통합
- [ ] Phase 5: tests/ 생성
- [ ] Snakefile 경로 업데이트
- [ ] Config 파일 경로 업데이트
- [ ] README.md 업데이트
- [ ] 문서 링크 업데이트
- [ ] Dry-run 테스트 성공
- [ ] 단위 테스트 성공
- [ ] docs/README.md 인덱스 생성
- [ ] CHANGELOG.md 업데이트
