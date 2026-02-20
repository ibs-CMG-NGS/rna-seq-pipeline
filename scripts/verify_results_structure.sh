#!/bin/bash
# 결과 폴더 표준화 검증 스크립트
# 사용법: bash scripts/verify_results_structure.sh

set -e

echo "======================================"
echo "RNA-seq Pipeline 결과 구조 검증"
echo "======================================"
echo ""

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 성공/실패 카운터
PASS=0
FAIL=0
WARN=0

# 검증 함수
check_exists() {
    local path="$1"
    local name="$2"
    local is_critical="${3:-true}"
    
    if [ -e "$path" ]; then
        echo -e "${GREEN}✓${NC} $name: $path"
        ((PASS++))
        return 0
    else
        if [ "$is_critical" = "true" ]; then
            echo -e "${RED}✗${NC} $name: $path (NOT FOUND)"
            ((FAIL++))
        else
            echo -e "${YELLOW}⚠${NC} $name: $path (OPTIONAL - NOT FOUND)"
            ((WARN++))
        fi
        return 1
    fi
}

check_not_empty() {
    local path="$1"
    local name="$2"
    
    if [ -e "$path" ] && [ -s "$path" ]; then
        local size=$(du -h "$path" | cut -f1)
        echo -e "${GREEN}✓${NC} $name: $path (size: $size)"
        ((PASS++))
        return 0
    else
        echo -e "${RED}✗${NC} $name: $path (EMPTY or NOT FOUND)"
        ((FAIL++))
        return 1
    fi
}

# ====================
# 1. 기본 디렉토리 구조
# ====================
echo "1. 기본 디렉토리 구조 확인"
echo "----------------------------------------"

check_exists "results" "Results 디렉토리"
check_exists "results/trimmed" "Trimmed 디렉토리"
check_exists "results/aligned" "Aligned 디렉토리"
check_exists "results/counts" "Counts 디렉토리"
check_exists "results/qc" "QC 디렉토리" "false"
check_exists "logs" "Logs 디렉토리"

echo ""

# ====================
# 2. 샘플 파일 검증
# ====================
echo "2. 샘플별 결과 파일 검증"
echo "----------------------------------------"

# Trimmed FASTQ 파일 확인
echo "📁 Trimmed FASTQ 파일:"
TRIMMED_COUNT=$(find results/trimmed -name "*.fastq.gz" 2>/dev/null | wc -l)
if [ $TRIMMED_COUNT -gt 0 ]; then
    echo -e "${GREEN}✓${NC} Trimmed FASTQ 파일: $TRIMMED_COUNT 개 발견"
    ((PASS++))
    # 샘플 목록
    ls -lh results/trimmed/*.fastq.gz | head -n 5
    if [ $TRIMMED_COUNT -gt 5 ]; then
        echo "  ... (총 $TRIMMED_COUNT 개)"
    fi
else
    echo -e "${RED}✗${NC} Trimmed FASTQ 파일이 없습니다"
    ((FAIL++))
fi
echo ""

# Aligned BAM 파일 확인
echo "📁 Aligned BAM 파일:"
BAM_SAMPLES=$(find results/aligned -name "*.bam" 2>/dev/null | wc -l)
if [ $BAM_SAMPLES -gt 0 ]; then
    echo -e "${GREEN}✓${NC} BAM 파일: $BAM_SAMPLES 개 발견"
    ((PASS++))
    # 각 샘플 디렉토리 확인
    for sample_dir in results/aligned/*/; do
        if [ -d "$sample_dir" ]; then
            sample=$(basename "$sample_dir")
            bam="$sample_dir/Aligned.sortedByCoord.out.bam"
            log="$sample_dir/Log.final.out"
            
            if [ -f "$bam" ]; then
                size=$(du -h "$bam" | cut -f1)
                echo -e "  ${GREEN}✓${NC} $sample: $bam ($size)"
                ((PASS++))
            else
                echo -e "  ${RED}✗${NC} $sample: BAM 파일 없음"
                ((FAIL++))
            fi
            
            if [ -f "$log" ]; then
                echo -e "    ${GREEN}✓${NC} Log.final.out 존재"
                ((PASS++))
            else
                echo -e "    ${YELLOW}⚠${NC} Log.final.out 없음"
                ((WARN++))
            fi
        fi
    done
else
    echo -e "${RED}✗${NC} BAM 파일이 없습니다"
    ((FAIL++))
fi
echo ""

# ====================
# 3. Counts Matrix 검증
# ====================
echo "3. Counts Matrix 검증"
echo "----------------------------------------"

check_not_empty "results/counts/counts_matrix.txt" "Raw counts matrix"

if [ -f "results/counts/counts_matrix_clean.csv" ]; then
    check_not_empty "results/counts/counts_matrix_clean.csv" "Clean counts matrix (CSV)"
    
    # CSV 파일 내용 미리보기
    echo ""
    echo "  📊 Counts matrix 미리보기 (처음 5줄):"
    head -n 5 results/counts/counts_matrix_clean.csv
    echo ""
    
    # 행/열 개수 확인
    GENE_COUNT=$(tail -n +2 results/counts/counts_matrix_clean.csv | wc -l)
    SAMPLE_COUNT=$(head -n 1 results/counts/counts_matrix_clean.csv | tr ',' '\n' | tail -n +2 | wc -l)
    echo -e "  ${GREEN}✓${NC} 유전자 수: $GENE_COUNT"
    echo -e "  ${GREEN}✓${NC} 샘플 수: $SAMPLE_COUNT"
    ((PASS+=2))
fi

if [ -f "results/counts/counts_matrix.txt.summary" ]; then
    check_exists "results/counts/counts_matrix.txt.summary" "featureCounts summary"
    echo ""
    echo "  📊 featureCounts summary:"
    cat results/counts/counts_matrix.txt.summary
fi

echo ""

# ====================
# 4. QC 리포트 검증
# ====================
echo "4. QC 리포트 검증"
echo "----------------------------------------"

# QC HTML 리포트
if [ -f "results/qc_report.html" ]; then
    check_not_empty "results/qc_report.html" "QC HTML 리포트"
    
    # HTML 파일 크기
    SIZE=$(du -h results/qc_report.html | cut -f1)
    echo "  크기: $SIZE"
    
    # 내용 확인 (샘플 개수)
    if command -v grep &> /dev/null; then
        SAMPLE_IN_REPORT=$(grep -o "sample-row" results/qc_report.html 2>/dev/null | wc -l)
        echo -e "  ${GREEN}✓${NC} 리포트 내 샘플 수: $SAMPLE_IN_REPORT"
        ((PASS++))
    fi
else
    echo -e "${YELLOW}⚠${NC} QC HTML 리포트가 없습니다 (선택사항)"
    ((WARN++))
fi

# QC JSON summary
if [ -f "results/qc/qc_summary.json" ]; then
    check_not_empty "results/qc/qc_summary.json" "QC JSON summary"
    
    if command -v python3 &> /dev/null; then
        echo ""
        echo "  📊 QC Summary 내용:"
        python3 -m json.tool results/qc/qc_summary.json 2>/dev/null | head -n 20
    fi
fi

echo ""

# ====================
# 5. FastQC 결과 검증
# ====================
echo "5. FastQC 결과 검증"
echo "----------------------------------------"

# Raw FastQC
RAW_FASTQC_COUNT=$(find results/qc -name "*_raw_fastqc.html" 2>/dev/null | wc -l)
if [ $RAW_FASTQC_COUNT -gt 0 ]; then
    echo -e "${GREEN}✓${NC} Raw FastQC 리포트: $RAW_FASTQC_COUNT 개"
    ((PASS++))
else
    echo -e "${YELLOW}⚠${NC} Raw FastQC 리포트 없음"
    ((WARN++))
fi

# Trimmed FastQC
TRIMMED_FASTQC_COUNT=$(find results/qc -name "*_trimmed_fastqc.html" 2>/dev/null | wc -l)
if [ $TRIMMED_FASTQC_COUNT -gt 0 ]; then
    echo -e "${GREEN}✓${NC} Trimmed FastQC 리포트: $TRIMMED_FASTQC_COUNT 개"
    ((PASS++))
else
    echo -e "${YELLOW}⚠${NC} Trimmed FastQC 리포트 없음"
    ((WARN++))
fi

# FastQC 자동 평가 결과
if [ -f "results/qc/fastqc_evaluation.json" ]; then
    check_not_empty "results/qc/fastqc_evaluation.json" "FastQC 자동 평가 결과"
    
    if command -v python3 &> /dev/null; then
        echo ""
        echo "  📊 FastQC 평가 요약:"
        python3 << 'EOF'
import json
import sys

try:
    with open("results/qc/fastqc_evaluation.json", "r") as f:
        data = json.load(f)
    
    # 전체 상태
    if "overall_status" in data:
        status = data["overall_status"]
        emoji = "✓" if status == "PASS" else ("⚠" if status == "WARN" else "✗")
        print(f"    {emoji} 전체 상태: {status}")
    
    # 샘플별 요약
    if "samples" in data:
        total = len(data["samples"])
        passed = sum(1 for s in data["samples"].values() if s.get("status") == "PASS")
        warned = sum(1 for s in data["samples"].values() if s.get("status") == "WARN")
        failed = sum(1 for s in data["samples"].values() if s.get("status") == "FAIL")
        
        print(f"    총 샘플: {total}")
        print(f"    PASS: {passed}, WARN: {warned}, FAIL: {failed}")
        
except Exception as e:
    print(f"    오류: {e}", file=sys.stderr)
EOF
    fi
fi

echo ""

# ====================
# 6. 로그 파일 검증
# ====================
echo "6. 로그 파일 검증"
echo "----------------------------------------"

# Cutadapt 로그
CUTADAPT_LOGS=$(find logs/cutadapt -name "*.log" 2>/dev/null | wc -l)
echo -e "Cutadapt 로그: $CUTADAPT_LOGS 개"

# STAR 로그
STAR_LOGS=$(find logs/star -name "*.log" 2>/dev/null | wc -l)
echo -e "STAR 로그: $STAR_LOGS 개"

# FastQC 로그
FASTQC_LOGS=$(find logs/fastqc -name "*.log" 2>/dev/null | wc -l)
echo -e "FastQC 로그: $FASTQC_LOGS 개"

# featureCounts 로그
if [ -f "logs/featurecounts.log" ]; then
    echo -e "${GREEN}✓${NC} featureCounts 로그 존재"
    ((PASS++))
fi

# QC 리포트 로그
if [ -f "logs/qc_report.log" ]; then
    echo -e "${GREEN}✓${NC} QC 리포트 로그 존재"
    ((PASS++))
fi

echo ""

# ====================
# 7. 디스크 사용량
# ====================
echo "7. 디스크 사용량"
echo "----------------------------------------"

echo "📊 디렉토리별 용량:"
du -h --max-depth=1 results 2>/dev/null | sort -rh

echo ""
echo "📊 상위 10개 대용량 파일:"
find results -type f -exec du -h {} + 2>/dev/null | sort -rh | head -n 10

echo ""

# ====================
# 최종 요약
# ====================
echo "======================================"
echo "검증 결과 요약"
echo "======================================"
echo -e "${GREEN}✓ PASS:${NC} $PASS"
echo -e "${YELLOW}⚠ WARN:${NC} $WARN"
echo -e "${RED}✗ FAIL:${NC} $FAIL"
echo ""

TOTAL=$((PASS + WARN + FAIL))
if [ $TOTAL -gt 0 ]; then
    SUCCESS_RATE=$((PASS * 100 / TOTAL))
    echo "성공률: ${SUCCESS_RATE}%"
fi

echo ""

if [ $FAIL -eq 0 ]; then
    echo -e "${GREEN}🎉 모든 필수 검증 통과!${NC}"
    if [ $WARN -gt 0 ]; then
        echo -e "${YELLOW}⚠️  $WARN 개의 선택사항 파일이 누락되었습니다.${NC}"
    fi
    echo ""
    echo "✅ 브랜치 병합(merge) 준비 완료"
    exit 0
else
    echo -e "${RED}❌ $FAIL 개의 필수 검증 실패${NC}"
    echo ""
    echo "⚠️  위 문제를 해결한 후 병합하세요."
    exit 1
fi
