# FastQC 자동 평가 가이드

## 📚 목차
- [개요](#개요)
- [자동 평가 기능](#자동-평가-기능)
- [평가 기준](#평가-기준)
- [사용 방법](#사용-방법)
- [결과 해석](#결과-해석)
- [설정 커스터마이징](#설정-커스터마이징)

---

## 개요

FastQC 자동 평가 기능은 모든 샘플의 FastQC 결과를 자동으로 분석하고, `FASTQC_GUIDE.md`에 명시된 기준에 따라 품질을 판단합니다. 이를 통해:

✅ **시간 절약**: 수십~수백 개의 FastQC HTML 리포트를 일일이 확인할 필요 없음  
✅ **객관적 판단**: 문서화된 기준에 따른 일관된 평가  
✅ **우선순위 지정**: 문제가 있는 샘플만 선별하여 집중 검토  
✅ **RNA-seq 특성 고려**: 정상적인 FAIL 패턴 자동 허용  

---

## 자동 평가 기능

### 🎯 주요 기능

1. **FastQC 데이터 파싱**
   - `fastqc_data.txt` 파일에서 수치 데이터 추출
   - 11개 모듈 상태 (PASS/WARN/FAIL) 분석

2. **RNA-seq 특화 평가**
   - 정상적인 FAIL 패턴 자동 허용:
     - Per base sequence content (랜덤 프라이밍 편향)
     - Sequence Duplication Levels (고발현 유전자)
     - Overrepresented sequences (rRNA, 하우스키핑 유전자)
     - Adapter Content (trimming 전)

3. **3단계 판정 시스템**
   - **PASS**: 분석 사용 가능
   - **WARN**: 경미한 문제, 리포트 확인 권장
   - **FAIL**: 심각한 문제, 재시퀀싱 또는 제외 고려

4. **요약 리포트 생성**
   - 텍스트 리포트: 사람이 읽기 쉬운 형식
   - JSON 파일: 프로그래밍 방식 처리 가능

---

## 평가 기준

### ✅ PASS 기준 (모두 충족 필요)

#### 1. Basic Statistics
- **Total Sequences** ≥ 1,000,000 (1M reads)
- **GC Content**: 30-70% (종-특이적)

#### 2. Per base sequence quality
- **Median Quality** ≥ Q28 (대부분의 위치)
- **Lower Quartile** ≥ Q20
- Q30+ 리드 비율 ≥ 75%

#### 3. Per base N content
- **N content** < 5% (모든 위치)
- Critical threshold: < 10%

#### 4. Adapter Content (trimming 후)
- **Adapter %** < 1%

### ⚠️ WARN 기준

- Median quality Q20-Q28
- N content 5-10%
- Q30+ 리드 비율 < 75%
- RNA-seq 정상 FAIL 외 모듈에서 FAIL

### ❌ FAIL 기준 (심각한 문제)

- Median quality < Q20
- Lower quartile < Q10
- N content > 10%
- Total sequences < 1M
- GC content < 30% 또는 > 70%
- Adapter content > 1% (trimming 후)

---

## 사용 방법

### 1️⃣ 파이프라인 통합 실행 (권장)

파이프라인을 실행하면 자동으로 FastQC 평가가 수행됩니다:

```bash
# 파이프라인 실행
snakemake --cores 12

# FastQC와 평가만 실행
snakemake --cores 4 results/qc/fastqc_evaluation.txt
```

**자동으로 생성되는 파일:**
```
results/qc/
├── *_fastqc.html                    # 개별 샘플 FastQC 리포트
├── *_fastqc.zip                     # FastQC 원본 데이터
├── fastqc_evaluation.txt            # ⭐ 자동 평가 리포트 (텍스트)
└── fastqc_evaluation.json           # 자동 평가 결과 (JSON)
```

### 2️⃣ 독립 실행

FastQC가 이미 실행된 경우, 평가만 따로 수행할 수 있습니다:

```bash
# 기본 실행
python3 src/evaluate_fastqc.py results/qc/

# 출력 경로 지정
python3 src/evaluate_fastqc.py results/qc/ \
    -o results/qc/my_evaluation.txt \
    --json results/qc/my_evaluation.json

# Trimming 후 데이터 평가 (어댑터 엄격 체크)
python3 src/evaluate_fastqc.py results/qc/ \
    --trimmed \
    -o results/qc/trimmed_evaluation.txt
```

### 3️⃣ 커스텀 설정 사용

```bash
# JSON 설정 파일 생성
cat > my_qc_config.json << EOF
{
  "min_total_sequences": 5000000,
  "min_median_quality": 30,
  "min_q30_percentage": 85
}
EOF

# 커스텀 설정으로 평가
python3 src/evaluate_fastqc.py results/qc/ \
    --config my_qc_config.json \
    -o results/qc/strict_evaluation.txt
```

---

## 결과 해석

### 📊 리포트 구조

자동 생성되는 `fastqc_evaluation.txt`는 다음 섹션으로 구성됩니다:

```
# FastQC 자동 평가 리포트

================================================================================

## 📊 요약

- 총 샘플 수: 24
- ✅ PASS: 18 샘플
- ⚠️ WARN: 4 샘플
- ❌ FAIL: 2 샘플
- 🔴 ERROR: 0 샘플

================================================================================

## ❌ FAIL - 즉시 확인 필요
[심각한 문제가 있는 샘플들 상세 정보]

## ⚠️ WARN - 리포트 확인 권장
[경미한 문제가 있는 샘플들 상세 정보]

## ✅ PASS - 분석 사용 가능
[정상 샘플 목록]

================================================================================

## 📝 참고사항
[RNA-seq 정상 FAIL 패턴 설명]
```

### 🔍 샘플별 상세 정보 예시

```markdown
### Sample_ABC_1

**권장사항**: REVIEW REQUIRED - 심각한 문제 발견. 재시퀀싱 또는 제외 고려

**문제점**:
- 🔴 Critical: Low quality at positions: 80-101 (Q17.5)
- 🔴 Total Sequences too low: 856,432 (minimum: 1,000,000)

**경고**:
- ⚠️ Warning: N content: 6.2% (recommended: <5.0%)

**기본 통계**:
- Total Sequences: 856,432
- GC Content: 47%
- Sequence Length: 101
```

### 📈 다음 단계 결정

#### ✅ PASS 샘플
```bash
# 다음 단계로 진행
# Trimming, Alignment, Quantification
```

#### ⚠️ WARN 샘플
```bash
# 1. FastQC HTML 리포트 확인
firefox results/qc/Sample_XYZ_1_fastqc.html

# 2. 경고 내용 검토
# - RNA-seq 정상 패턴인가?
# - Trimming으로 개선 가능한가?

# 3. 진행 여부 결정
# - 대부분의 경우: 진행 가능
# - 의심스러우면: PI와 상의
```

#### ❌ FAIL 샘플
```bash
# 1. FastQC HTML 리포트 상세 확인
firefox results/qc/Sample_FAIL_1_fastqc.html

# 2. 원인 분석
# - 시퀀싱 품질 문제? → 재시퀀싱
# - 샘플 오염? → 원본 샘플 확인
# - 라이브러리 문제? → 라이브러리 재제작

# 3. 조치
# - 재시퀀싱 요청
# - 또는 분석에서 제외
```

---

## 설정 커스터마이징

### 📝 Config 파일 편집

`config_H2O2_astrocyte.yaml` (또는 사용 중인 config 파일):

```yaml
# === FastQC Auto-Evaluation Parameters ===
fastqc_evaluation:
  enabled: true  # 자동 평가 활성화/비활성화
  
  # Basic Statistics thresholds
  min_total_sequences: 1000000  # 최소 리드 수
  min_gc_content: 30            # 최소 GC%
  max_gc_content: 70            # 최대 GC%
  
  # Quality Score thresholds
  min_median_quality: 28        # 최소 median Phred score
  min_lower_quartile: 20        # 최소 lower quartile
  min_q30_percentage: 75        # 최소 Q30+ 비율 (%)
  
  # N content thresholds
  max_n_content: 5.0           # 허용 N% (경고)
  critical_n_content: 10.0     # Critical N% (FAIL)
  
  # Adapter content thresholds
  max_adapter_trimmed: 1.0     # Trimming 후 최대 adapter%
  warn_adapter_raw: 10.0       # Raw data adapter% 경고
  
  # Output files
  evaluation_report: "fastqc_evaluation.txt"
  evaluation_json: "fastqc_evaluation.json"
```

### 🔧 프로젝트별 맞춤 설정

#### 예제 1: 고품질 요구사항

```yaml
fastqc_evaluation:
  min_total_sequences: 5000000  # 5M reads 이상
  min_median_quality: 30        # Q30 이상
  min_q30_percentage: 85        # 85% Q30+
  max_n_content: 1.0           # N < 1%
```

#### 예제 2: 사람 게놈 (GC 범위 조정)

```yaml
fastqc_evaluation:
  min_gc_content: 38  # 사람: ~41%
  max_gc_content: 44
```

#### 예제 3: 낮은 Input 샘플

```yaml
fastqc_evaluation:
  min_total_sequences: 500000  # 0.5M reads 허용
  min_q30_percentage: 70       # 70% Q30+ 허용
```

### 🎛️ 평가 비활성화

자동 평가를 사용하지 않으려면:

```yaml
fastqc_evaluation:
  enabled: false
```

또는 Snakefile에서 해당 규칙만 제외:

```bash
# FastQC만 실행, 평가는 스킵
snakemake --cores 4 \
    --until fastqc_raw
```

---

## 실전 예제

### 시나리오 1: 정상 프로젝트

```
## 📊 요약
- 총 샘플 수: 12
- ✅ PASS: 12 샘플
- ⚠️ WARN: 0 샘플
- ❌ FAIL: 0 샘플
```

**조치**: 모든 샘플로 분석 진행 ✅

---

### 시나리오 2: 일부 문제 샘플

```
## 📊 요약
- 총 샘플 수: 24
- ✅ PASS: 20 샘플
- ⚠️ WARN: 3 샘플
- ❌ FAIL: 1 샘플
```

**조치**:
1. FAIL 샘플 (Sample_18): FastQC HTML 확인 → 재시퀀싱 요청
2. WARN 샘플 (Sample_03, 07, 15): 리포트 확인
   - RNA-seq 정상 패턴 → 진행
   - 실제 문제 → PI와 상의

---

### 시나리오 3: 전반적 품질 문제

```
## 📊 요약
- 총 샘플 수: 12
- ✅ PASS: 2 샘플
- ⚠️ WARN: 5 샘플
- ❌ FAIL: 5 샘플
```

**조치**:
1. 시퀀싱 런(run) 전체 문제 의심
2. 시퀀싱 센터에 문의
3. 전체 재시퀀싱 고려

---

## JSON 출력 활용

### Python에서 결과 읽기

```python
import json

# JSON 결과 로드
with open('results/qc/fastqc_evaluation.json', 'r') as f:
    results = json.load(f)

# FAIL 샘플만 필터링
fail_samples = [
    r['sample'] for r in results 
    if r['status'] == 'FAIL'
]

print(f"FAIL 샘플: {fail_samples}")

# 샘플별 문제점 출력
for result in results:
    if result['status'] == 'FAIL':
        print(f"\n{result['sample']}:")
        for issue in result['issues']:
            print(f"  - {issue}")
```

### R에서 결과 읽기

```r
library(jsonlite)

# JSON 결과 로드
results <- fromJSON('results/qc/fastqc_evaluation.json')

# 상태별 샘플 수
table(sapply(results, function(x) x$status))

# FAIL 샘플 추출
fail_samples <- results[sapply(results, function(x) x$status == 'FAIL')]

# 데이터프레임으로 변환
df <- data.frame(
  sample = sapply(results, function(x) x$sample),
  status = sapply(results, function(x) x$status),
  total_seq = sapply(results, function(x) x$basic_stats$`Total Sequences`),
  gc_content = sapply(results, function(x) x$basic_stats$`%GC`)
)
```

---

## 문제 해결

### Q1: 모든 샘플이 FAIL로 나옵니다

**원인**:
- 평가 기준이 너무 엄격함
- 프로젝트 특성이 기본 설정과 다름

**해결**:
```yaml
# config 파일에서 기준 완화
fastqc_evaluation:
  min_total_sequences: 500000  # 1M → 0.5M
  min_median_quality: 25       # 28 → 25
```

---

### Q2: RNA-seq 정상 FAIL도 문제로 표시됩니다

**확인**:
```bash
# 리포트에서 "참고사항" 섹션 확인
# RNA-seq 정상 FAIL 패턴은 자동으로 제외됨
```

스크립트가 다음 모듈은 자동으로 허용합니다:
- Per base sequence content
- Sequence Duplication Levels
- Overrepresented sequences
- Adapter Content (trimming 전)

---

### Q3: 평가 리포트가 생성되지 않습니다

**확인 사항**:
1. FastQC가 완료되었는지 확인
```bash
ls results/qc/*_fastqc.zip
```

2. Config에서 활성화 여부 확인
```yaml
fastqc_evaluation:
  enabled: true  # false가 아닌지 확인
```

3. 로그 확인
```bash
cat logs/fastqc/evaluate_raw.log
```

---

## 참고 자료

- **FastQC 해석 가이드**: `FASTQC_GUIDE.md`
- **파이프라인 가이드**: `PIPELINE_GUIDE.md`
- **QC 리포트 가이드**: `QC_REPORT_GUIDE.md`

---

**문서 버전**: 1.0  
**최종 수정**: 2026-01-26  
**작성자**: RNA-seq Pipeline Team
