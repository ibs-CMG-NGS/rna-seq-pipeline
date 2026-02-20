# 결과 폴더 검증 명령어 모음
# 서버에서 복사해서 사용하세요

# ==========================================
# 1. 빠른 구조 확인 (트리 구조)
# ==========================================

# 기본 트리 (2단계 깊이)
tree -L 2 results

# 전체 트리 (파일 크기 포함)
tree -h results

# 디렉토리만 표시
tree -d results


# ==========================================
# 2. 주요 파일 존재 확인
# ==========================================

# 모든 BAM 파일 나열
find results/aligned -name "*.bam" -ls

# Counts matrix 확인
ls -lh results/counts/

# QC 리포트 확인
ls -lh results/*.html results/qc/*.json 2>/dev/null

# Trimmed FASTQ 개수
ls results/trimmed/*.fastq.gz | wc -l


# ==========================================
# 3. 파일 크기 및 개수 확인
# ==========================================

# 디렉토리별 용량
du -h --max-depth=1 results | sort -rh

# 파일 타입별 개수
echo "BAM files: $(find results/aligned -name '*.bam' | wc -l)"
echo "Trimmed FASTQ: $(find results/trimmed -name '*.fastq.gz' | wc -l)"
echo "FastQC reports: $(find results/qc -name '*_fastqc.html' | wc -l)"

# 총 용량
du -sh results


# ==========================================
# 4. Counts Matrix 검증
# ==========================================

# CSV 파일 헤더 확인 (샘플 이름)
head -n 1 results/counts/counts_matrix_clean.csv

# 유전자 개수
tail -n +2 results/counts/counts_matrix_clean.csv | wc -l

# 샘플 개수
head -n 1 results/counts/counts_matrix_clean.csv | tr ',' '\n' | tail -n +2 | wc -l

# 처음 5개 유전자 데이터
head -n 6 results/counts/counts_matrix_clean.csv | column -t -s,


# ==========================================
# 5. BAM 파일 상세 확인
# ==========================================

# 각 샘플별 BAM 파일 크기
for dir in results/aligned/*/; do
    sample=$(basename "$dir")
    if [ -f "$dir/Aligned.sortedByCoord.out.bam" ]; then
        size=$(du -h "$dir/Aligned.sortedByCoord.out.bam" | cut -f1)
        echo "$sample: $size"
    fi
done

# STAR 정렬 통계 확인 (Log.final.out)
for dir in results/aligned/*/; do
    sample=$(basename "$dir")
    echo "=== $sample ==="
    if [ -f "$dir/Log.final.out" ]; then
        grep "Uniquely mapped reads %" "$dir/Log.final.out"
        grep "Number of input reads" "$dir/Log.final.out"
    fi
done


# ==========================================
# 6. QC 결과 확인
# ==========================================

# QC HTML 리포트 크기
ls -lh results/qc_report.html

# QC JSON summary 내용 (Python 사용)
python3 -m json.tool results/qc/qc_summary.json | head -n 30

# FastQC 평가 결과 (있다면)
if [ -f "results/qc/fastqc_evaluation.json" ]; then
    python3 -m json.tool results/qc/fastqc_evaluation.json | head -n 50
fi


# ==========================================
# 7. featureCounts Summary 확인
# ==========================================

# Summary 파일 내용 전체
cat results/counts/counts_matrix.txt.summary

# 샘플별 할당률 계산
python3 << 'EOF'
with open("results/counts/counts_matrix.txt.summary", "r") as f:
    lines = f.readlines()
    
header = lines[0].strip().split('\t')
samples = header[1:]  # 첫 컬럼 제외

assigned = lines[1].strip().split('\t')[1:]  # Assigned 행

print("샘플별 할당된 reads:")
for sample, count in zip(samples, assigned):
    print(f"  {sample}: {int(count):,}")
EOF


# ==========================================
# 8. 로그 파일 확인
# ==========================================

# 최신 에러 확인
find logs -name "*.log" -exec grep -l "error\|Error\|ERROR\|failed\|Failed" {} \;

# Snakemake 로그 확인
ls -lt .snakemake/log/*.snakemake.log | head -n 1

# 최신 Snakemake 로그 내용
latest_log=$(ls -t .snakemake/log/*.snakemake.log | head -n 1)
echo "=== Latest Snakemake Log: $latest_log ==="
tail -n 50 "$latest_log"


# ==========================================
# 9. 표준화 체크리스트 (수동 확인용)
# ==========================================

cat << 'CHECKLIST'
□ results/trimmed/ 에 모든 샘플의 FASTQ 파일이 있는가?
□ results/aligned/ 에 각 샘플별 디렉토리가 있는가?
□ 각 샘플 디렉토리에 Aligned.sortedByCoord.out.bam 파일이 있는가?
□ results/counts/counts_matrix_clean.csv 파일이 존재하고 비어있지 않은가?
□ results/qc_report.html 파일이 생성되었는가? (선택사항)
□ results/qc/fastqc_evaluation.json 파일이 있는가? (Phase3 feature)
□ 모든 샘플이 counts matrix에 포함되어 있는가?
□ BAM 파일 크기가 합리적인가? (너무 작지 않은가?)
□ Log.final.out 에서 mapping rate가 적절한가? (>70% 권장)
□ 로그 파일에 치명적 에러가 없는가?
CHECKLIST


# ==========================================
# 10. 종합 검증 스크립트 실행
# ==========================================

# 자동 검증 스크립트 실행 (위에서 생성한 스크립트)
bash scripts/verify_results_structure.sh

# 또는 실행 권한 부여 후
chmod +x scripts/verify_results_structure.sh
./scripts/verify_results_structure.sh
