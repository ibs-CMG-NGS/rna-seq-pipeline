# Snakefile

# --- 설정 파일 로드 ---
configfile: "config_H2O2_astrocyte.yaml"

# --- 1. 전역 변수 설정 ---
# 경로 변수 구성
DATA_DIR = config.get("data_dir", "data")
RAW_DATA_DIR = f"{DATA_DIR}/{config.get('raw_data_subdir', 'raw')}"

RESULTS_DIR = config.get("results_dir", "results")
TRIMMED_DIR = f"{RESULTS_DIR}/{config.get('trimmed_subdir', 'trimmed')}"
ALIGNED_DIR = f"{RESULTS_DIR}/{config.get('aligned_subdir', 'aligned')}"
COUNTS_DIR = f"{RESULTS_DIR}/{config.get('counts_subdir', 'counts')}"
QC_DIR = f"{RESULTS_DIR}/{config.get('qc_subdir', 'qc')}"

LOGS_DIR = config.get("logs_dir", "logs")

# 샘플 이름 정의
SAMPLES, = glob_wildcards(f"{RAW_DATA_DIR}/{{sample}}_1.fastq.gz")

# 경로 변수 (config.yaml에서 로드)
STAR_INDEX = config["star_index"]
ANNOTATION_GTF = config["annotation_gtf"]


# --- 2. 파이프라인의 최종 목표 정의 (Rule all) ---
# QC 리포트 생성 여부에 따라 최종 목표 설정
rule all:
    input:
        f"{COUNTS_DIR}/counts_matrix.txt",
        f"{COUNTS_DIR}/counts_matrix_clean.csv",  # DE analysis용 clean matrix 추가
        f"{RESULTS_DIR}/{config.get('qc_report_filename', 'qc_report.html')}" if config.get("generate_qc_report", True) else []


# --- 3. Raw 데이터 QC (FastQC) 규칙 ---
# (이전과 동일)
rule fastqc_raw:
    input:
        f"{RAW_DATA_DIR}/{{sample}}_{{read}}.fastq.gz"
    output:
        html=f"{QC_DIR}/{{sample}}_{{read}}_raw_fastqc.html",
        zip=f"{QC_DIR}/{{sample}}_{{read}}_raw_fastqc.zip"
    log:
        f"{LOGS_DIR}/fastqc/{{sample}}_{{read}}_raw.log"
    threads: config.get("fastqc_threads", 2)
    shell:
        f"fastqc {{input}} -o {QC_DIR}/ > {{log}} 2>&1"


# --- 4. 어댑터 제거 (cutadapt) 규칙 ---
rule cutadapt:
    input:
        r1=f"{RAW_DATA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DATA_DIR}/{{sample}}_2.fastq.gz"
    output:
        r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz"
    log:
        f"{LOGS_DIR}/cutadapt/{{sample}}.log"
    threads: config.get("cutadapt_threads", 4)
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
        r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz"
    output:
        bam=f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam"
    log:
        f"{LOGS_DIR}/star/{{sample}}.log"
    threads: config["star_threads"]
    resources:
        mem_gb=config.get("star_memory_gb", 35)  # Memory limit per STAR job (GB)
    shell:
        f"""
        STAR --runThreadN {{threads}} \
             --genomeDir {STAR_INDEX} \
             --readFilesIn {{input.r1}} {{input.r2}} \
             --readFilesCommand zcat \
             --outFileNamePrefix {ALIGNED_DIR}/{{wildcards.sample}}/ \
             --outSAMtype BAM SortedByCoordinate \
             --limitBAMsortRAM {{config[star_sort_memory_bytes]}} \
             --outBAMsortingThreadN {{threads}} > {{log}} 2>&1
        """

# --- 6. 유전자 발현량 계산 (featureCounts) 규칙 ---
rule featurecounts_quant:
    input:
        bams=expand(f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        counts=f"{COUNTS_DIR}/counts_matrix.txt",
        summary=f"{COUNTS_DIR}/counts_matrix.txt.summary"
    log:
        f"{LOGS_DIR}/featurecounts.log"
    threads: config["featurecounts_threads"]
    shell:
        f"""
        featureCounts -T {{threads}} -p \
                      -t {{config[feature_type]}} \
                      -g {{config[attribute_type]}} \
                      -s {{config[strandedness]}} \
                      -a {ANNOTATION_GTF} \
                      -o {{output.counts}} \
                      {{input.bams}} > {{log}} 2>&1
        """

# --- 6-1. Count matrix 변환 (DE 분석용) ---
rule convert_counts_matrix:
    input:
        f"{COUNTS_DIR}/counts_matrix.txt"
    output:
        f"{COUNTS_DIR}/counts_matrix_clean.csv"
    log:
        f"{LOGS_DIR}/convert_counts.log"
    shell:
        """
        python3 src/convert_counts_matrix.py {input} {output} > {log} 2>&1
        """

# --- 7. QC 리포트 생성 규칙 ---
rule generate_qc_report:
    input:
        counts=f"{COUNTS_DIR}/counts_matrix.txt",
        counts_summary=f"{COUNTS_DIR}/counts_matrix.txt.summary",
        cutadapt_logs=expand(f"{LOGS_DIR}/cutadapt/{{sample}}.log", sample=SAMPLES),
        star_logs=expand(f"{ALIGNED_DIR}/{{sample}}/Log.final.out", sample=SAMPLES)
    output:
        f"{RESULTS_DIR}/{config.get('qc_report_filename', 'qc_report.html')}"
    params:
        top_genes=config.get("qc_top_genes", 10)
    log:
        f"{LOGS_DIR}/qc_report.log"
    script:
        "src/generate_qc_report.py"