# Snakefile

# --- 설정 파일 로드 ---
configfile: "config.yaml"

# --- 1. 전역 변수 설정 ---
# 샘플 이름 정의
SAMPLES, = glob_wildcards("data/raw/{sample}_1.fastq.gz")

# 경로 변수 (config.yaml에서 로드)
STAR_INDEX = config["star_index"]
ANNOTATION_GTF = config["annotation_gtf"]


# --- 2. 파이프라인의 최종 목표 정의 (Rule all) ---
# QC 리포트 생성 여부에 따라 최종 목표 설정
rule all:
    input:
        "results/counts/counts_matrix.txt",
        config["qc_report_output"] if config.get("generate_qc_report", True) else []


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
rule cutadapt:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4
    shell:
        """
        cutadapt -a {config[adapter_r1]} \
                 -A {config[adapter_r2]} \
                 --minimum-length {config[min_read_length]} \
                 --quality-cutoff {config[quality_cutoff]} \
                 --cores={threads} \
                 -o {output.r1} -p {output.r2} \
                 {input.r1} {input.r2} > {log} 2>&1
        """

# --- 5. STAR 정렬 규칙 ---
rule star_align:
    input:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz"
    output:
        bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}.log"
    threads: config["star_threads"]
    resources:
        mem_gb=35  # Memory limit per STAR job (35GB for 62GB total RAM)
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {STAR_INDEX} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix results/aligned/{wildcards.sample}/ \
             --outSAMtype BAM SortedByCoordinate \
             --limitBAMsortRAM 30000000000 \
             --outBAMsortingThreadN {threads} > {log} 2>&1
        """

# --- 6. 유전자 발현량 계산 (featureCounts) 규칙 ---
rule featurecounts_quant:
    input:
        bams=expand("results/aligned/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        "results/counts/counts_matrix.txt"
    log:
        "logs/featurecounts.log"
    threads: config["featurecounts_threads"]
    shell:
        """
        featureCounts -T {threads} -p \
                      -t {config[feature_type]} \
                      -g {config[attribute_type]} \
                      -s {config[strandedness]} \
                      -a {ANNOTATION_GTF} \
                      -o {output} \
                      {input.bams} > {log} 2>&1
        """

# --- 7. QC 리포트 생성 규칙 ---
rule generate_qc_report:
    input:
        counts="results/counts/counts_matrix.txt",
        counts_summary="results/counts/counts_matrix.txt.summary",
        cutadapt_logs=expand("logs/cutadapt/{sample}.log", sample=SAMPLES),
        star_logs=expand("results/aligned/{sample}/Log.final.out", sample=SAMPLES)
    output:
        config["qc_report_output"]
    params:
        top_genes=config.get("qc_top_genes", 10)
    log:
        "logs/qc_report.log"
    script:
        "src/generate_qc_report.py"