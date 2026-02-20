# 결과 검증 가이드

파이프라인 실행 후 결과가 올바르게 생성되었는지 확인하는 방법입니다.

## 빠른 시작

### 방법 1: 자동 검증 스크립트 (추천)

**Bash 스크립트:**
```bash
bash scripts/verify_results_structure.sh
```

**Python 스크립트:**
```bash
python scripts/verify_results.py
```

두 스크립트 모두:
- ✅ 필수 파일/디렉토리 존재 확인
- ✅ 파일 크기 및 내용 검증
- ✅ 샘플별 결과 확인
- ✅ 최종 성공/실패 판정

### 방법 2: 수동 명령어

`scripts/VERIFICATION_COMMANDS.md` 파일에 있는 명령어들을 복사해서 사용하세요.

---

## 검증 항목

### 1. 필수 검증 항목 (PASS 필요)

#### 디렉토리 구조
```
✓ results/
  ✓ trimmed/      - Cutadapt 결과
  ✓ aligned/      - STAR 정렬 결과
  ✓ counts/       - featureCounts 결과
  ✓ qc/           - QC 결과 (선택)
✓ logs/           - 로그 파일
```

#### 파일 존재
- `results/trimmed/*.fastq.gz` - 모든 샘플의 trimmed FASTQ
- `results/aligned/{sample}/Aligned.sortedByCoord.out.bam` - 각 샘플의 BAM 파일
- `results/counts/counts_matrix.txt` - Raw counts matrix
- `results/counts/counts_matrix_clean.csv` - Clean counts matrix (CSV)
- `results/counts/counts_matrix.txt.summary` - featureCounts summary

#### 내용 검증
- Counts matrix에 모든 샘플이 포함되어 있는가?
- BAM 파일 크기가 합리적인가? (최소 수백 MB)
- 유전자 개수가 적절한가? (수만 개 수준)

### 2. 선택 검증 항목 (WARN)

- `results/qc_report.html` - QC HTML 리포트
- `results/qc/qc_summary.json` - QC JSON summary
- `results/qc/fastqc_evaluation.json` - FastQC 자동 평가 (Phase3 feature)
- FastQC HTML 리포트들
- STAR Log.final.out 파일들

---

## 검증 스크립트 출력 예시

### 성공한 경우
```
======================================
RNA-seq Pipeline 결과 구조 검증
======================================

1. 기본 디렉토리 구조 확인
----------------------------------------
✓ Results 디렉토리: results
✓ Trimmed 디렉토리: results/trimmed
✓ Aligned 디렉토리: results/aligned
✓ Counts 디렉토리: results/counts
✓ Logs 디렉토리: logs

2. 샘플별 결과 파일 검증
----------------------------------------
📁 Trimmed FASTQ 파일:
✓ Trimmed FASTQ 파일: 12 개 발견
...

3. Counts Matrix 검증
----------------------------------------
✓ Raw counts matrix: results/counts/counts_matrix.txt (25.3MB)
✓ Clean counts matrix (CSV): results/counts/counts_matrix_clean.csv (26.1MB)

  📊 Counts matrix 정보:
  ✓ 유전자 수: 60,483
  ✓ 샘플 수: 6
...

======================================
검증 결과 요약
======================================
✓ PASS: 45
⚠ WARN: 3
✗ FAIL: 0

성공률: 93%

🎉 모든 필수 검증 통과!
⚠️  3 개의 선택사항 파일이 누락되었습니다.

✅ 브랜치 병합(merge) 준비 완료
```

### 실패한 경우
```
======================================
검증 결과 요약
======================================
✓ PASS: 32
⚠ WARN: 5
✗ FAIL: 3

성공률: 80%

❌ 3 개의 필수 검증 실패

⚠️  위 문제를 해결한 후 병합하세요.
```

---

## 수동 확인 명령어

### 기본 구조 확인
```bash
# 디렉토리 트리
tree -L 2 results

# 디렉토리별 용량
du -h --max-depth=1 results | sort -rh
```

### 파일 개수 확인
```bash
# BAM 파일
find results/aligned -name "*.bam" | wc -l

# Trimmed FASTQ
ls results/trimmed/*.fastq.gz | wc -l

# FastQC 리포트
find results/qc -name "*_fastqc.html" | wc -l
```

### Counts Matrix 확인
```bash
# 헤더 (샘플 이름)
head -n 1 results/counts/counts_matrix_clean.csv

# 유전자 개수
tail -n +2 results/counts/counts_matrix_clean.csv | wc -l

# 샘플 개수
head -n 1 results/counts/counts_matrix_clean.csv | tr ',' '\n' | tail -n +2 | wc -l

# 처음 5개 유전자
head -n 6 results/counts/counts_matrix_clean.csv | column -t -s,
```

### BAM 파일 확인
```bash
# 각 샘플별 BAM 크기
for dir in results/aligned/*/; do
    sample=$(basename "$dir")
    if [ -f "$dir/Aligned.sortedByCoord.out.bam" ]; then
        size=$(du -h "$dir/Aligned.sortedByCoord.out.bam" | cut -f1)
        echo "$sample: $size"
    fi
done
```

### STAR 정렬 통계
```bash
# Mapping rate 확인
for dir in results/aligned/*/; do
    sample=$(basename "$dir")
    echo "=== $sample ==="
    if [ -f "$dir/Log.final.out" ]; then
        grep "Uniquely mapped reads %" "$dir/Log.final.out"
        grep "Number of input reads" "$dir/Log.final.out"
    fi
done
```

### featureCounts Summary
```bash
# Summary 전체 내용
cat results/counts/counts_matrix.txt.summary

# 샘플별 할당된 reads 개수
python3 << 'EOF'
with open("results/counts/counts_matrix.txt.summary", "r") as f:
    lines = f.readlines()
    
header = lines[0].strip().split('\t')
samples = header[1:]

assigned = lines[1].strip().split('\t')[1:]

print("샘플별 할당된 reads:")
for sample, count in zip(samples, assigned):
    print(f"  {sample}: {int(count):,}")
EOF
```

### QC 결과 확인
```bash
# QC HTML 리포트
ls -lh results/qc_report.html

# FastQC 평가 결과
python3 -m json.tool results/qc/fastqc_evaluation.json | head -n 50
```

---

## 문제 해결

### BAM 파일이 너무 작은 경우
```bash
# 각 샘플의 STAR 로그 확인
cat results/aligned/{sample}/Log.final.out

# Mapping rate 확인
grep "Uniquely mapped reads %" results/aligned/{sample}/Log.final.out
```

**해결 방법:**
- Mapping rate가 낮은 경우: 참조 유전체 확인
- Input reads가 적은 경우: 원본 FASTQ 확인

### Counts matrix가 비어있는 경우
```bash
# featureCounts 로그 확인
cat logs/featurecounts.log

# Summary 파일 확인
cat results/counts/counts_matrix.txt.summary
```

**해결 방법:**
- "Assigned" 비율이 0인 경우: GTF 파일 확인
- Strandedness 문제: config.yaml의 featureCounts 파라미터 확인

### FastQC 평가가 FAIL인 경우
```bash
# 평가 결과 상세 확인
python3 -m json.tool results/qc/fastqc_evaluation.json

# 실패한 모듈 확인
python3 << 'EOF'
import json
with open("results/qc/fastqc_evaluation.json") as f:
    data = json.load(f)
    
for sample, info in data["samples"].items():
    if info["status"] != "PASS":
        print(f"{sample}: {info['status']}")
        print(f"  Failed: {info.get('failed_modules', [])}")
        print(f"  Warned: {info.get('warned_modules', [])}")
EOF
```

---

## 브랜치 병합 체크리스트

병합하기 전에 다음을 확인하세요:

- [ ] `bash scripts/verify_results_structure.sh` 실행 결과 PASS
- [ ] `python scripts/verify_results.py` 실행 결과 성공 (exit code 0)
- [ ] 모든 샘플의 BAM 파일이 존재하고 크기가 합리적
- [ ] Counts matrix에 모든 샘플이 포함됨
- [ ] 유전자 개수가 적절 (수만 개 수준)
- [ ] Mapping rate가 적절 (>70% 권장)
- [ ] 로그 파일에 치명적 에러 없음
- [ ] QC 리포트 생성됨 (선택사항)

**모든 체크리스트를 통과했다면:**
```bash
git add .
git commit -m "feat: FastQC auto-evaluation feature completed"
git push origin feature/fastqc-auto-evaluation

# GitHub에서 Pull Request 생성 후 merge
```

---

## 추가 정보

자세한 명령어는 `scripts/VERIFICATION_COMMANDS.md` 참조
