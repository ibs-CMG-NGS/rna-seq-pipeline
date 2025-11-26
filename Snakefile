# Snakefile

# --- 1. 전역 변수 설정 ---
# 샘플 이름 정의
SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")

# 경로 변수 (수정 필요!)
STAR_INDEX = "genome/star_index/"
ANNOTATION_GTF = "/path/to/your/annotation.gtf" # GTF 파일 경로 추가


# --- 2. 파이프라인의 최종 목표 정의 (Rule all) ---
# 최종 목표를 featureCounts 결과물인 counts matrix로 변경합니다.
rule all:
    input:
        "results/counts/counts_matrix.txt"


# --- 3. Raw 데이터 QC (FastQC) 규칙 ---
# (이전과 동일)
rule fastqc_raw:
    input:
        "data/raw/{sample}_{read}.fastq.gz"
    output:
        html="results/qc/{sample}_{read}_raw_fastqc.html",
        zip="results/qc/{sample}_{read}_raw_fastqc.zip"
    log:
        "logs/fastqc/{sample}_{read}_raw.log"
    threads: 1
    shell:
        "fastqc {input} -o results/qc/ > {log} 2>&1"


# --- 4. 어댑터 제거 (cutadapt) 규칙 ---
# (이전과 동일)
rule cutadapt:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4
    shell:
        """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                 --cores={threads} \
                 -o {output.r1} -p {output.r2} \
                 {input.r1} {input.r2} > {log} 2>&1
        """

# --- 5. STAR 정렬 규칙 ---
# (이전과 동일)
rule star_align:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz"
    output:
        # directory()는 삭제하고, 실제 BAM 파일 경로를 명시해주는 것이
        # 다음 단계인 featureCounts에서 사용하기 더 명확합니다.
        bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}.log"
    threads: 8
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {STAR_INDEX} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix results/aligned/{wildcards.sample}/ \
             --outSAMtype BAM SortedByCoordinate > {log} 2>&1
        """

# --- 6. 유전자 발현량 계산 (featureCounts) 규칙 (⭐ 새로 추가된 부분) ---
rule featurecounts_quant:
    input:
        # 모든 샘플의 BAM 파일을 입력으로 받습니다.
        bams=expand("results/aligned/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        "results/counts/counts_matrix.txt"
    log:
        "logs/featurecounts.log"
    threads: 8
    shell:
        """
        featureCounts -T {threads} -p -t exon -g gene_id \
                      -a {ANNOTATION_GTF} \
                      -o {output} \
                      {input.bams} > {log} 2>&1
        """