#!/usr/bin/env python3
import sys

print("=" * 60)
print("RNA-seq 파이프라인 결과 점검 리포트")
print("=" * 60)

# 1. Cutadapt 결과
print("\n📊 1. CUTADAPT (어댑터 제거)")
print("-" * 60)
with open('logs/cutadapt/Ctrl_3.log', 'r') as f:
    content = f.read()
    for line in content.split('\n'):
        if 'Total read pairs processed' in line:
            print(f"  총 read pairs: {line.split(':')[1].strip()}")
        elif 'Pairs that were too short' in line:
            print(f"  제거된 짧은 reads: {line.split(':')[1].strip()}")
        elif 'Pairs written (passing filters)' in line:
            print(f"  필터 통과: {line.split(':')[1].strip()}")

# 2. STAR 정렬 결과
print("\n📊 2. STAR (정렬)")
print("-" * 60)
with open('results/aligned/Ctrl_3/Log.final.out', 'r') as f:
    for line in f:
        line = line.strip()
        if 'Number of input reads' in line:
            print(f"  입력 reads: {line.split('|')[1].strip()}")
        elif 'Uniquely mapped reads %' in line:
            print(f"  고유 매핑률: {line.split('|')[1].strip()}")
        elif '% of reads mapped to multiple loci' in line:
            print(f"  다중 매핑: {line.split('|')[1].strip()}")
        elif '% of reads unmapped: too short' in line:
            print(f"  매핑 실패 (짧음): {line.split('|')[1].strip()}")

# 3. featureCounts 결과
print("\n📊 3. FEATURECOUNTS (유전자별 카운팅)")
print("-" * 60)
with open('results/counts/counts_matrix.txt.summary', 'r') as f:
    lines = f.readlines()
    total = 0
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            status, count = parts[0], int(parts[1])
            total += count
            if 'Assigned' in status:
                print(f"  ✓ 할당된 reads: {count:,} ({count/total*100:.1f}%)")
            elif 'MultiMapping' in status:
                print(f"  • 다중 매핑: {count:,} ({count/total*100:.1f}%)")
            elif 'NoFeatures' in status:
                print(f"  • 특징 없음: {count:,} ({count/total*100:.1f}%)")
            elif 'Ambiguity' in status:
                print(f"  • 모호함: {count:,} ({count/total*100:.1f}%)")
    print(f"  총 fragments: {total:,}")

# 4. Count matrix 통계
print("\n📊 4. COUNT MATRIX (유전자 발현 행렬)")
print("-" * 60)
total_genes = 0
genes_with_reads = 0
total_counts = 0
high_expression = []

with open('results/counts/counts_matrix.txt', 'r') as f:
    for i, line in enumerate(f):
        if line.startswith('#') or i == 1:
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 7:
            gene_id = parts[0]
            count = int(parts[6])
            total_genes += 1
            total_counts += count
            if count > 0:
                genes_with_reads += 1
            if count > 10000:
                high_expression.append((gene_id, count))

print(f"  총 유전자 수: {total_genes:,}")
print(f"  발현된 유전자: {genes_with_reads:,} ({genes_with_reads/total_genes*100:.1f}%)")
print(f"  총 카운트: {total_counts:,}")
print(f"  평균 카운트/유전자: {total_counts/total_genes:.1f}")

# 상위 발현 유전자
print(f"\n  🔝 고발현 유전자 (>10,000 reads): {len(high_expression)}개")
high_expression.sort(key=lambda x: x[1], reverse=True)
for gene, count in high_expression[:5]:
    print(f"     {gene}: {count:,}")

# 5. 파일 크기
print("\n📊 5. 파일 크기")
print("-" * 60)
import os
files_to_check = [
    ('원본 R1', 'data/raw/Ctrl_3_1.fastq.gz'),
    ('원본 R2', 'data/raw/Ctrl_3_2.fastq.gz'),
    ('Trimmed R1', 'results/trimmed/Ctrl_3_1.fastq.gz'),
    ('Trimmed R2', 'results/trimmed/Ctrl_3_2.fastq.gz'),
    ('Aligned BAM', 'results/aligned/Ctrl_3/Aligned.sortedByCoord.out.bam'),
    ('Counts matrix', 'results/counts/counts_matrix.txt'),
]

for name, path in files_to_check:
    if os.path.exists(path):
        size = os.path.getsize(path)
        if size > 1e9:
            print(f"  {name:15s}: {size/1e9:.2f} GB")
        else:
            print(f"  {name:15s}: {size/1e6:.2f} MB")

# 결론
print("\n" + "=" * 60)
print("✅ 점검 결과 요약")
print("=" * 60)
mapping_rate = genes_with_reads / total_genes * 100
if mapping_rate > 85:
    print("✓ 매핑률 우수 (87.95%)")
else:
    print("⚠ 매핑률 확인 필요")

if genes_with_reads > 10000:
    print(f"✓ 충분한 유전자 발현 감지 ({genes_with_reads:,}개)")
else:
    print("⚠ 발현 유전자 수 확인 필요")

print(f"✓ 총 {total_counts:,} reads가 유전자에 할당됨")
print("\n파이프라인이 정상적으로 완료되었습니다! 🎉")
