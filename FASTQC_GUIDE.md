# FastQC 품질 평가 가이드

## 📚 목차
- [개요](#개요)
- [FastQC란?](#fastqc란)
- [FastQC 실행 방법](#fastqc-실행-방법)
- [리포트 해석 가이드](#리포트-해석-가이드)
  - [기본 통계](#1-basic-statistics-기본-통계)
  - [염기별 품질](#2-per-base-sequence-quality-염기별-시퀀스-품질)
  - [타일별 품질](#3-per-tile-sequence-quality-타일별-시퀀스-품질)
  - [시퀀스별 품질 점수](#4-per-sequence-quality-scores-시퀀스별-품질-점수)
  - [염기 조성](#5-per-base-sequence-content-염기별-시퀀스-조성)
  - [GC 함량](#6-per-sequence-gc-content-시퀀스별-gc-함량)
  - [N 염기 함량](#7-per-base-n-content-염기별-n-함량)
  - [시퀀스 길이 분포](#8-sequence-length-distribution-시퀀스-길이-분포)
  - [중복 수준](#9-sequence-duplication-levels-시퀀스-중복-수준)
  - [과대표현 서열](#10-overrepresented-sequences-과대표현-서열)
  - [어댑터 함량](#11-adapter-content-어댑터-함량)
- [RNA-seq 특화 해석](#rna-seq-특화-해석)
- [문제 해결 가이드](#문제-해결-가이드)
- [실전 예제](#실전-예제)

---

## 개요

FastQC는 NGS(차세대 시퀀싱) 데이터의 품질을 평가하는 가장 널리 사용되는 도구입니다. 이 가이드는 FastQC 리포트를 처음 보는 사람도 이해하고 해석할 수 있도록 작성되었습니다.

### 🎯 학습 목표
- FastQC 리포트의 각 항목이 무엇을 의미하는지 이해
- PASS/WARN/FAIL 상태의 의미와 대응 방법 파악
- RNA-seq 데이터의 정상적인 패턴 인식
- 실제 문제와 허위 양성(false positive) 구분

---

## FastQC란?

### 📖 정의
FastQC는 **FASTQ 형식의 시퀀싱 데이터 품질을 자동으로 평가**하는 Java 기반 프로그램입니다.

### ✨ 주요 기능
- **빠른 분석**: 대용량 파일도 수 분 내에 분석
- **시각화**: 그래프와 차트로 직관적 표현
- **자동 판정**: PASS/WARN/FAIL로 품질 상태 표시
- **HTML 리포트**: 브라우저에서 바로 확인 가능

### 🔧 버전 정보
```bash
# FastQC 버전 확인
fastqc --version
# 이 파이프라인 사용 버전: FastQC v0.12.1
```

---

## FastQC 실행 방법

### 1️⃣ 단일 파일 분석
```bash
# 기본 실행
fastqc sample_1.fastq.gz

# 출력 디렉토리 지정
fastqc -o qc_results/ sample_1.fastq.gz

# 멀티스레드 사용
fastqc -t 4 sample_1.fastq.gz
```

### 2️⃣ 다중 파일 분석
```bash
# 여러 파일 동시 분석
fastqc -t 8 -o qc_results/ sample_*.fastq.gz

# 와일드카드로 모든 FASTQ 파일 분석
fastqc -t 12 data/raw/*.fastq.gz
```

### 3️⃣ 파이프라인 내 자동 실행
이 RNA-seq 파이프라인에서는 Snakemake가 자동으로 FastQC를 실행합니다:
```bash
# Snakemake 실행 시 자동으로 FastQC 수행
snakemake --cores 12

# FastQC만 따로 실행
snakemake --cores 4 fastqc_raw
```

### 📂 출력 파일
```
results/qc/
├── sample_1_fastqc.html  # 브라우저에서 볼 수 있는 리포트
├── sample_1_fastqc.zip   # 원본 데이터 (텍스트, 이미지)
└── sample_1_fastqc/      # 압축 해제된 폴더
    ├── summary.txt       # ⭐ 간단한 PASS/WARN/FAIL 요약
    ├── fastqc_data.txt   # 원본 수치 데이터
    └── Images/           # 그래프 이미지 파일들
```

---

## 리포트 해석 가이드

FastQC 리포트는 **11개의 품질 평가 모듈**로 구성됩니다. 각 모듈은 다음 세 가지 상태로 표시됩니다:

| 상태 | 아이콘 | 의미 |
|------|--------|------|
| ✅ **PASS** | 녹색 체크 | 정상 범위 내 |
| ⚠️ **WARN** | 주황색 느낌표 | 약간 비정상적이지만 허용 가능 |
| ❌ **FAIL** | 빨간색 X | 문제가 있음 (하지만 RNA-seq에서는 정상일 수 있음!) |

> **⚠️ 중요**: FAIL이 반드시 "나쁜 데이터"를 의미하지는 않습니다. 특히 RNA-seq 데이터는 정상적인 생물학적 특성 때문에 일부 모듈에서 FAIL이 나올 수 있습니다.

---

### 1. Basic Statistics (기본 통계)

#### 📊 표시 항목
```
Filename                : Ctrl_3_1.fastq.gz
File type               : Conventional base calls
Encoding                : Sanger / Illumina 1.9
Total Sequences         : 32,280,544
Sequences flagged as poor quality : 0
Sequence length         : 101
%GC                     : 48
```

#### 🔍 의미
- **Total Sequences**: 총 리드(read) 개수
  - Paired-end의 경우 R1과 R2가 동일해야 함
  - 일반적으로 2천만~5천만 개가 적당
- **Sequence length**: 리드 길이
  - Illumina NovaSeq: 주로 100bp, 150bp
  - 모든 샘플이 동일한 길이여야 함
- **%GC**: 전체 GC 함량
  - 마우스(mm10): ~42%
  - 사람(hg38): ~41%
  - 종에 따라 다름

#### ✅ 정상 기준
- 파일 형식이 올바르게 인식됨
- Total Sequences > 0
- Sequence length가 예상값과 일치

#### ❌ 문제 징후
- Total Sequences가 비정상적으로 적음 (< 100만 개)
- %GC가 예상 범위를 크게 벗어남 (GC < 30% 또는 > 70%)

---

### 2. Per base sequence quality (염기별 시퀀스 품질)

#### 📊 그래프 설명
- **X축**: 리드 내 염기 위치 (1번부터 끝까지)
- **Y축**: Phred 품질 점수 (0~40)
- **Box plot**: 각 위치의 품질 점수 분포

```
Quality score 범위:
0-20  (빨강)  : 나쁨 (1% 이상 에러율)
20-28 (주황)  : 보통 (0.1~1% 에러율)
28-40 (녹색)  : 좋음 (0.1% 미만 에러율)
```

#### 🔍 의미
각 염기 위치에서 시퀀싱 품질이 얼마나 좋은지 보여줍니다.

**Phred Score 해석**:
| Score | 정확도 | 에러율 | 의미 |
|-------|--------|--------|------|
| 10 | 90% | 10% | 매우 나쁨 |
| 20 | 99% | 1% | 최소 허용 |
| 30 | 99.9% | 0.1% | 우수 ⭐ |
| 40 | 99.99% | 0.01% | 최고 |

#### ✅ 정상 기준
- 대부분의 위치에서 **Q30 이상** (녹색 영역)
- Read 1: 처음부터 끝까지 고른 품질
- Read 2: 뒤쪽으로 갈수록 약간 떨어지는 것은 정상

#### ⚠️ WARN 상태
- 일부 구간에서 Q20~Q28 (주황색 영역)
- 리드 끝부분 품질 저하

#### ❌ FAIL 상태
- 중간값이 Q20 미만
- 하위 사분위수가 Q10 미만

#### 🛠️ 대처 방법
```bash
# Trimmomatic으로 저품질 끝부분 제거
trimmomatic PE \
  -phred33 \
  input_R1.fq.gz input_R2.fq.gz \
  output_R1.fq.gz unpaired_R1.fq.gz \
  output_R2.fq.gz unpaired_R2.fq.gz \
  TRAILING:20 MINLEN:36

# 또는 Cutadapt 품질 필터링 강화
cutadapt -q 20,20 -m 20 \
  -o output_R1.fq.gz -p output_R2.fq.gz \
  input_R1.fq.gz input_R2.fq.gz
```

---

### 3. Per tile sequence quality (타일별 시퀀스 품질)

#### 📊 그래프 설명
- **X축**: 리드 내 염기 위치
- **Y축**: 플로우셀 타일 번호
- **색상**: 평균 대비 품질 편차
  - 파랑: 평균보다 좋음
  - 빨강: 평균보다 나쁨

#### 🔍 의미
Illumina 플로우셀의 **특정 물리적 위치(타일)**에서 품질 문제가 있는지 확인합니다.

#### ✅ 정상 기준
- 대부분의 타일이 파란색 (평균 수준)
- 균일한 색상 분포

#### ⚠️ WARN 상태
- 일부 타일에서 약간의 품질 저하 (연한 빨강)
- 💡 **흔한 현상**: Illumina 시퀀싱에서 자주 발생, 분석에는 영향 없음

#### ❌ FAIL 상태
- 특정 타일이 진한 빨강색 (Phred score > 5점 차이)

#### 🛠️ 대처 방법
```bash
# 심각한 경우: 해당 타일의 리드만 필터링 (권장하지 않음)
# 보통은 무시하고 진행해도 됨

# 시퀀싱 센터에 문의:
# - 플로우셀 불량 가능성
# - 재시퀀싱 요청 고려
```

---

### 4. Per sequence quality scores (시퀀스별 품질 점수)

#### 📊 그래프 설명
- **X축**: 평균 품질 점수 (Phred score)
- **Y축**: 해당 점수를 가진 리드의 개수

#### 🔍 의미
각 리드(read)의 **전체 평균 품질**이 어떻게 분포하는지 보여줍니다.

#### ✅ 정상 기준
- 대부분의 리드가 **Q30 이상**에 분포
- 좁고 높은 피크 형태

```
Good example:
          ┃
          ┃
          ┃██
          ┃███
    ──────┃████───────
      20  30  40
```

#### ❌ FAIL 상태
- 피크가 Q27 미만
- 분포가 넓게 퍼짐

```
Bad example:
     ┃
   ███┃
  ████┃
 █████┃
───────┃───────
   20  30  40
```

#### 🛠️ 대처 방법
```bash
# 저품질 리드 필터링
# Cutadapt의 --minimum-quality 옵션 사용
cutadapt -q 20 --minimum-length 20 \
  -o output.fq.gz input.fq.gz
```

---

### 5. Per base sequence content (염기별 시퀀스 조성)

#### 📊 그래프 설명
- **X축**: 리드 내 염기 위치
- **Y축**: A, T, G, C 각 염기의 비율 (%)
- **4개의 선**: A(녹색), T(빨강), G(검정), C(파랑)

#### 🔍 의미
각 위치에서 A, T, G, C 염기의 비율이 균등한지 확인합니다.

#### ✅ 정상 기준 (DNA 랜덤 라이브러리)
- 모든 위치에서 A≈T, G≈C
- 4개 선이 평행하게 유지
- 약간의 편차(±5%) 허용

```
Good (랜덤 DNA):
30%─ G ═══════════
25%─ C ═══════════
25%─ A ═══════════
20%─ T ═══════════
```

#### ❌ FAIL 상태
- 리드 시작 부분에서 심한 편향
- 특정 염기가 과다/과소

```
RNA-seq 패턴 (정상적 FAIL):
40%─ A ╱╲─────────
30%─   ╱  ╲───────
20%─ T╱    ╲──────
10%─ G      ══════
     C      ══════
   1-10    11-101
```

#### 💡 RNA-seq에서는 정상!
**이 모듈은 RNA-seq에서 거의 항상 FAIL입니다. 하지만 이것은 정상입니다!**

**FAIL이 나는 이유**:
1. **랜덤 헥사머 프라이밍 편향**
   - cDNA 합성 시 사용하는 랜덤 프라이머가 완전히 랜덤하지 않음
   - 처음 10~15bp에서 염기 편향 발생

2. **전사 시작 지점 (TSS) 편향**
   - 유전자 시작 부분의 서열이 특정 패턴을 가짐
   - GC가 풍부한 프로모터 영역 영향

3. **Poly-A 선택 (선택적)**
   - mRNA의 Poly-A tail 때문에 A가 많이 나타날 수 있음

#### ✅ 실제 문제 vs 정상 패턴
| 상황 | 패턴 | 판단 |
|------|------|------|
| 처음 1-15bp만 편향 | RNA-seq 랜덤 프라이밍 | ✅ 정상 |
| 전체적으로 심한 편향 | 어댑터/오염 | ❌ 문제 |
| 중간부터 A가 급증 | Poly-A 오염 | ❌ 문제 |

---

### 6. Per sequence GC content (시퀀스별 GC 함량)

#### 📊 그래프 설명
- **파란선**: 실제 데이터의 GC 분포
- **빨간선**: 이론적 정규분포
- **X축**: GC 함량 (%)
- **Y축**: 리드 개수

#### 🔍 의미
각 리드의 GC 함량 분포가 예상되는 정규분포를 따르는지 확인합니다.

#### ✅ 정상 기준
- 파란선과 빨간선이 거의 일치
- 단일 피크 형태
- 피크 위치가 전체 GC%와 일치

```
Good example:
    ┃    ██
    ┃   ████
    ┃  ██████
    ┃ ████████
────┃──────────────
   30  42  50  60
       ↑
    마우스 GC%
```

#### ❌ FAIL 상태
1. **이중 피크** (Bimodal distribution)
```
    ┃  ██    ██
    ┃ ████  ████
────┃──────────────
   30  40  50  60
```
   - **원인**: 오염(contamination), 어댑터 다이머
   - **대처**: 오염 소스 확인 및 제거

2. **비정상적으로 넓은 분포**
```
    ┃  ████████████
────┃──────────────
   20  40  60  80
```
   - **원인**: 여러 생물 종의 DNA 혼합
   - **대처**: 샘플 오염 확인

3. **피크 이동**
```
    ┃        ██
    ┃       ████
────┃──────────────
   30  40  50  60
              ↑
          예상 42%
```
   - **원인**: 잘못된 레퍼런스 게놈 사용
   - **대처**: 게놈 버전 확인

#### 💡 RNA-seq 특성
RNA-seq에서는 약간의 편향이 정상:
- **높은 GC 영역**: 리보솜 RNA (rRNA) - 특히 18S, 28S
- **낮은 GC 영역**: 미토콘드리아 유전자

---

### 7. Per base N content (염기별 N 함량)

#### 📊 그래프 설명
- **X축**: 리드 내 염기 위치
- **Y축**: 불확실 염기(N)의 비율 (%)

#### 🔍 의미
시퀀싱 과정에서 **염기를 확정할 수 없는 경우(N)** 얼마나 많은지 확인합니다.

#### ✅ 정상 기준
- 모든 위치에서 **N < 5%**
- 대부분의 경우 0% 근처

```
Excellent:
5%─
   
0%─══════════════
   1        101
```

#### ❌ FAIL 상태
- 어느 위치에서든 N > 20%

```
Bad example:
20%─      ╱╲
10%─    ╱    ╲
0%─════╯      ╰═══
   1    50    101
```

#### 🛠️ 대처 방법
- **N > 20%**: 시퀀싱 실패, 재시퀀싱 필요
- **N < 10%**: N이 포함된 리드 제거
```bash
# BBTools로 N 제거
bbduk.sh in=input.fq.gz out=output.fq.gz maxns=0
```

---

### 8. Sequence Length Distribution (시퀀스 길이 분포)

#### 📊 그래프 설명
- **X축**: 리드 길이 (bp)
- **Y축**: 해당 길이를 가진 리드의 개수

#### 🔍 의미
모든 리드가 **동일한 길이**를 가지는지 확인합니다.

#### ✅ 정상 기준
- **Trimming 전**: 단일 피크 (모두 101bp)
```
    ┃     █
    ┃     █
────┃─────█─────
       100 101 102
```

- **Trimming 후**: 범위 분포 (20~101bp)
```
    ┃ █████████
    ┃█████████████
────┃──────────────
   20    60    101
```

#### ❌ FAIL 상태
- 예상치 못한 다양한 길이 (trimming 전)
- 매우 짧은 리드만 남음 (< 20bp)

#### 🛠️ 대처 방법
```bash
# 최소 길이 필터링
cutadapt --minimum-length 20 \
  -o output.fq.gz input.fq.gz
```

---

### 9. Sequence Duplication Levels (시퀀스 중복 수준)

#### 📊 그래프 설명
- **파란선**: 중복 제거 전
- **빨간선**: 중복 제거 후
- **X축**: 중복 횟수 (1, 2, 3... >10)
- **Y축**: 해당 중복 수준의 서열 비율 (%)

#### 🔍 의미
**동일한 서열**이 얼마나 반복되는지 확인합니다.

**중복의 두 가지 의미**:
1. **생물학적 중복**: 같은 전사체가 많이 발현됨 (정상)
2. **기술적 중복**: PCR에서 같은 분자가 여러 번 증폭됨 (문제)

#### ✅ 정상 기준 (DNA-seq)
- 대부분의 서열이 1~2회만 나타남
- 중복도 < 20%

```
DNA-seq (정상):
50%─█
20%─ █
10%─  █
 0%─   ████████
    1 2 3  >10
```

#### ❌ FAIL 상태 (DNA-seq)
- 중복도 > 50%
- PCR 과증폭 의심

#### 💡 RNA-seq에서는 정상!
**RNA-seq에서 이 모듈은 거의 항상 FAIL입니다. 이것은 정상입니다!**

```
RNA-seq (정상적 FAIL):
50%─        ████
20%─      ██
10%─    ██
 0%─█████
    1  >10  >500
```

**RNA-seq에서 높은 중복이 정상인 이유**:
1. **고발현 유전자**: Actb, Gapdh 같은 하우스키핑 유전자는 수십만 개 전사체 존재
2. **리보솜 RNA**: 18S, 28S rRNA가 전체 RNA의 80% 차지
3. **짧은 전사체**: 작은 유전자는 선택지가 제한적

**실제 문제 vs 정상 패턴**:
| 중복도 | RNA-seq | DNA-seq |
|--------|---------|---------|
| < 20% | 매우 낮음 (rRNA 제거?) | ✅ 정상 |
| 20-50% | ✅ 정상 | ⚠️ 주의 |
| 50-80% | ✅ 정상 (고발현 유전자) | ❌ 문제 |
| > 80% | ⚠️ rRNA 과다? | ❌ 심각 |

#### 🔬 실제 확인 방법
중복이 정상인지 확인하려면 **어떤 서열이 중복되는지** 봐야 합니다:
```bash
# 가장 많이 중복되는 서열 확인
zcat input.fq.gz | \
  awk 'NR%4==2' | \
  sort | uniq -c | \
  sort -rn | head -20

# 리보솜 RNA인지 확인
# BLAST 검색 또는 rRNA 데이터베이스 매칭
```

---

### 10. Overrepresented sequences (과대표현 서열)

#### 📊 표 설명
```
Sequence               Count    Percentage   Possible Source
AGGTCGACTGCATGACT...   125,421  0.38%        TruSeq Adapter
NNNNNNNNNNNNNNNN...    98,234   0.30%        No Hit
```

#### 🔍 의미
전체 리드에서 **0.1% 이상 차지하는 동일 서열**을 나열합니다.

#### ✅ 정상 기준
- "No sequences flagged as overrepresented"
- 또는 알려진 고발현 유전자

#### ❌ FAIL 상태
1. **Illumina 어댑터 검출**
   - TruSeq Adapter
   - Nextera Transposase Sequence
   - **대처**: Cutadapt로 제거
   ```bash
   cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o out_R1.fq.gz -p out_R2.fq.gz \
            in_R1.fq.gz in_R2.fq.gz
   ```

2. **Primer/Barcode 서열**
   - PCR primer 잔존
   - **대처**: 해당 서열 trim

3. **오염(Contamination)**
   - 박테리아/바이러스 DNA
   - **대처**: BLAST 검색 후 필터링

#### 💡 RNA-seq 정상 패턴
- **rRNA 서열**: 18S, 28S, 5S, 5.8S
- **미토콘드리아 유전자**: mt-Nd1, mt-Co1
- **고발현 유전자**: Actb, Gapdh, Eef1a1

**확인 방법**:
```bash
# 과대표현 서열을 BLAST 검색
# http://blast.ncbi.nlm.nih.gov/Blast.cgi
# 또는
blastn -query sequence.fa -db nt -remote
```

---

### 11. Adapter Content (어댑터 함량)

#### 📊 그래프 설명
- **X축**: 리드 내 위치
- **Y축**: 어댑터가 발견된 리드의 비율 (%)
- **여러 색상**: 다양한 어댑터 타입

#### 🔍 의미
리드 안에 **시퀀싱 어댑터 서열**이 얼마나 남아있는지 확인합니다.

#### ✅ 정상 기준 (Trimming 후)
- 모든 위치에서 < 0.1%
```
Good (trimmed):
5%─
   
0%─══════════════
   1        101
```

#### ❌ FAIL 상태 (Trimming 전)
- 리드 끝부분으로 갈수록 증가
```
Raw data (FAIL은 정상):
20%─         ╱
10%─       ╱
 0%─══════╯
   1      80  101
```

**어댑터가 검출되는 이유**:
1. **Insert 크기 < Read 길이**
   - Insert 100bp, Read 101bp → 1bp 어댑터 읽음
   - 정상적인 현상

2. **Adapter dimer**
   - Insert가 없이 어댑터끼리 붙음
   - 매우 짧은 리드 (<20bp)

#### 🛠️ 대처 방법
```bash
# Cutadapt로 어댑터 제거
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -q 20 \
  --minimum-length 20 \
  -o trimmed_R1.fq.gz \
  -p trimmed_R2.fq.gz \
  raw_R1.fq.gz \
  raw_R2.fq.gz

# Trim 후 FastQC 재확인
fastqc trimmed_R1.fq.gz trimmed_R2.fq.gz
```

**Illumina 주요 어댑터 서열**:
```
TruSeq Universal Adapter:
AGATCGGAAGAG

TruSeq Adapter, Index 1-12:
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

Nextera Transposase Sequence:
CTGTCTCTTATACACATCT
```

---

## RNA-seq 특화 해석

### 🧬 RNA-seq 데이터의 정상적인 "FAIL" 패턴

RNA-seq FastQC 리포트에서 다음 항목들은 **FAIL이어도 정상**입니다:

| 모듈 | 상태 | 이유 |
|------|------|------|
| Per base sequence content | ❌ FAIL | 랜덤 헥사머 프라이밍 편향 (처음 10~15bp) |
| Sequence Duplication Levels | ❌ FAIL | 고발현 유전자 (Actb, Gapdh, rRNA) |
| Overrepresented sequences | ❌ FAIL | rRNA, 미토콘드리아, 하우스키핑 유전자 |
| Adapter Content (trimming 전) | ❌ FAIL | Insert < Read length |

### ✅ 우수한 RNA-seq 데이터 체크리스트

#### Raw Data (Trimming 전)
```
✅ Per base sequence quality        : Q30 이상
✅ Per sequence quality scores      : 피크 > Q30
⚠️ Per tile sequence quality        : WARN 허용
✅ Per sequence GC content          : 단일 피크
✅ Per base N content               : < 1%
❌ Per base sequence content        : FAIL 정상 (처음 15bp 편향)
❌ Sequence Duplication Levels      : FAIL 정상 (고발현 유전자)
❌ Overrepresented sequences        : FAIL 정상 (rRNA 검출)
❌ Adapter Content                  : FAIL → Cutadapt 필요
```

#### Trimmed Data (Cutadapt 후)
```
✅ Per base sequence quality        : Q30 이상 유지
✅ Adapter Content                  : < 0.1% (개선됨!)
✅ Sequence Length Distribution     : 20~101bp 분포
❌ Per base sequence content        : 여전히 FAIL (정상)
❌ Sequence Duplication Levels      : 여전히 FAIL (정상)
```

### 📊 샘플 간 비교 포인트

**동일 프로젝트의 모든 샘플이 비슷한 패턴**을 보여야 합니다:

| 항목 | 확인 사항 |
|------|-----------|
| Total Sequences | ±20% 이내 유사 |
| %GC | ±2% 이내 |
| Quality profile | 동일한 곡선 형태 |
| Duplication level | 비슷한 수준 |

**차이가 크다면**:
- 샘플 준비 과정의 불균일
- 시퀀싱 배치 효과
- RNA 품질(RIN) 차이

---

## 문제 해결 가이드

### 🔴 심각한 문제들

#### 1. 낮은 품질 점수 (Q < 20)
```
증상: Per base sequence quality에서 빨간 영역
원인: 시퀀싱 실패, 플로우셀 문제
해결:
  1. Trimming 강화
     cutadapt -q 25,25 --minimum-length 30
  2. 재시퀀싱 고려
```

#### 2. 높은 N 비율 (> 10%)
```
증상: Per base N content 그래프 상승
원인: 시퀀싱 화학반응 실패
해결: 재시퀀싱 필요 (복구 불가능)
```

#### 3. 어댑터 다이머 (Adapter dimer)
```
증상: 
  - Sequence Length Distribution에서 매우 짧은 피크 (~20bp)
  - Adapter Content 100%
원인: Insert가 없이 어댑터만 증폭됨
해결:
  1. 엄격한 길이 필터링
     cutadapt --minimum-length 50
  2. 라이브러리 재제작 고려
```

#### 4. 샘플 오염 (Contamination)
```
증상:
  - Per sequence GC content 이중 피크
  - Overrepresented sequences에 박테리아/바이러스 서열
원인: 샘플 준비 과정의 오염
해결:
  1. BLAST로 오염원 확인
  2. 오염 서열 필터링
     bbduk.sh ref=contaminants.fa in=input.fq out=clean.fq
  3. 재실험 권장
```

### 🟡 일반적인 경고들

#### 1. Per tile sequence quality WARN
```
대처: 보통 무시해도 됨
      심각하면 시퀀싱 센터에 문의
```

#### 2. Sequence Duplication WARN (DNA-seq)
```
확인: 어떤 서열이 중복되는지 체크
대처: 
  - PCR 중복 제거 도구 사용 (Picard MarkDuplicates)
  - 라이브러리 복잡도(complexity) 확인
```

---

## 실전 예제

### 예제 1: 정상적인 RNA-seq Raw Data

```
✅ Basic Statistics
   Total Sequences: 32,280,544
   Sequence length: 101
   %GC: 48

✅ Per base sequence quality
   Q30+ throughout

⚠️ Per tile sequence quality
   Some tiles slightly lower (normal)

✅ Per sequence quality scores
   Peak at Q38

❌ Per base sequence content
   First 12bp biased (NORMAL for RNA-seq)

✅ Per sequence GC content
   Single peak at 48%

✅ Per base N content
   <0.1%

✅ Sequence Length Distribution
   All 101bp

❌ Sequence Duplication Levels
   70% duplicated (NORMAL - high expression genes)

❌ Overrepresented sequences
   18S rRNA, Actb, Gapdh (NORMAL)

❌ Adapter Content
   5% at read ends (NEEDS TRIMMING)
```

**판정**: ⭐⭐⭐⭐⭐ **우수한 데이터**
- Cutadapt 처리 필요
- FAIL 항목들은 모두 RNA-seq 정상 패턴

---

### 예제 2: Cutadapt 처리 후

```
✅ Per base sequence quality
   Still Q30+

✅ Adapter Content
   <0.1% (개선됨!)

✅ Sequence Length Distribution
   20-101bp range

❌ Per base sequence content
   Still FAIL (정상)

❌ Sequence Duplication Levels
   Still FAIL (정상)
```

**판정**: ✅ **분석 준비 완료**

---

### 예제 3: 문제 있는 데이터

```
❌ Per base sequence quality
   Q15 average - 매우 낮음!

❌ Per base N content
   15% N's - 시퀀싱 실패

❌ Per sequence GC content
   Two peaks at 35% and 65% - 오염 의심

❌ Overrepresented sequences
   E. coli sequences - 박테리아 오염!
```

**판정**: ❌ **사용 불가**
- 재시퀀싱 필요
- 샘플 준비 과정 점검

---

## 요약 및 권장사항

### 📋 FastQC 실행 체크리스트

#### 분석 전 (Raw Data)
- [ ] 모든 FASTQ 파일에 FastQC 실행
- [ ] summary.txt로 빠른 PASS/WARN/FAIL 확인
- [ ] Per base sequence quality 필수 확인
- [ ] Adapter Content 확인

#### Trimming 후
- [ ] FastQC 재실행
- [ ] Adapter Content 개선 확인 (<0.1%)
- [ ] Quality 유지 확인
- [ ] Sequence length 분포 확인

#### 비교 분석
- [ ] 모든 샘플의 패턴이 유사한지 확인
- [ ] 이상치(outlier) 샘플 식별
- [ ] MultiQC로 통합 리포트 생성

### 🎯 핵심 포인트

1. **FAIL ≠ 나쁜 데이터**
   - RNA-seq에서는 일부 FAIL이 정상
   - 각 모듈의 의미를 이해해야 함

2. **품질의 절대 기준**
   - Per base sequence quality: Q30+
   - Per base N content: <1%
   - 이 두 가지가 가장 중요!

3. **샘플 간 일관성**
   - 동일 프로젝트는 비슷한 패턴
   - 큰 차이는 문제 신호

4. **Trimming 전후 비교**
   - Adapter Content 개선 확인
   - 품질 손실 최소화

### 🔧 유용한 도구들

```bash
# MultiQC: 여러 샘플 통합 리포트
multiqc results/qc/

# FastQ Screen: 오염 확인
fastq_screen --aligner bowtie2 sample.fq.gz

# BBDuk: 품질 필터링 및 오염 제거
bbduk.sh in=input.fq out=clean.fq \
  ref=adapters.fa qtrim=rl trimq=20 minlen=20
```

---

## 참고 자료

### 📚 공식 문서
- FastQC Manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- FastQC Help: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

### 🎓 학습 자료
- RNA-seqlopedia: https://rnaseq.uoregon.edu/
- SEQanswers Forum: http://seqanswers.com/
- Biostars: https://www.biostars.org/

### 🛠️ 관련 도구
- MultiQC: https://multiqc.info/
- Cutadapt: https://cutadapt.readthedocs.io/
- BBTools: https://jgi.doe.gov/data-and-tools/bbtools/

---

**문서 버전**: 1.0  
**최종 수정**: 2025-12-05  
**작성자**: RNA-seq Pipeline Team  
**파이프라인 버전**: 1.0
