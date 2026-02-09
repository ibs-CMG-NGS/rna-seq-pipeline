# Snakefile

# --- 설정 파일 로드 ---
# 기본값: config/default.yaml, 실행 시 --configfile 옵션으로 오버라이드 가능
# 예: snakemake --configfile config/projects/H2O2_human_2025.yaml -j 12
configfile: "config/default.yaml"

# --- 1. 전역 변수 설정 ---
import os

# 입력 데이터 경로
DATA_DIR = config.get("data_dir", "data")
RAW_DATA_SUBDIR = config.get("raw_data_subdir", "raw")
RAW_DATA_DIR = f"{DATA_DIR}/{RAW_DATA_SUBDIR}" if RAW_DATA_SUBDIR else DATA_DIR

# 표준 디렉토리 구조 사용 여부
USE_STANDARD = config.get("use_standard_structure", False)
PROJECT_ID = config.get("project_id", "default_project")
PIPELINE_TYPE = config.get("pipeline_type", "rna-seq")

# 디렉토리 경로 함수
def get_sample_dir(sample_id):
    """샘플별 디렉토리 경로 반환"""
    if USE_STANDARD:
        base = config.get("base_results_dir", "/home/ngs/data/results")
        return f"{base}/{PROJECT_ID}/{sample_id}/{PIPELINE_TYPE}"
    else:
        return config.get("results_dir", "results")

def get_final_outputs_dir(sample_id):
    """Final outputs 디렉토리 경로"""
    if USE_STANDARD:
        return f"{get_sample_dir(sample_id)}/final_outputs"
    else:
        return config.get("results_dir", "results")

def get_intermediate_dir(sample_id):
    """Intermediate 디렉토리 경로"""
    if USE_STANDARD:
        return f"{get_sample_dir(sample_id)}/intermediate"
    else:
        return config.get("results_dir", "results")

def get_metadata_dir(sample_id):
    """Metadata 디렉토리 경로"""
    if USE_STANDARD:
        return f"{get_sample_dir(sample_id)}/metadata"
    else:
        return config.get("results_dir", "results")

# Legacy 구조 경로 (USE_STANDARD=False인 경우)
if not USE_STANDARD:
    RESULTS_DIR = config.get("results_dir", "results")
    TRIMMED_DIR = f"{RESULTS_DIR}/{config.get('trimmed_subdir', 'trimmed')}"
    ALIGNED_DIR = f"{RESULTS_DIR}/{config.get('aligned_subdir', 'aligned')}"
    COUNTS_DIR = f"{RESULTS_DIR}/{config.get('counts_subdir', 'counts')}"
    QC_DIR = f"{RESULTS_DIR}/{config.get('qc_subdir', 'qc')}"
    LOGS_DIR = config.get("logs_dir", "logs")
else:
    # 표준 구조 경로
    BASE_RESULTS = config.get("base_results_dir", "/home/ngs/data/results")
    PROJECT_DIR = f"{BASE_RESULTS}/{PROJECT_ID}"
    PROJECT_SUMMARY_DIR = f"{PROJECT_DIR}/project_summary"
    METADATA_DIR = f"{PROJECT_DIR}/metadata"
    
    # 프로젝트 레벨 디렉토리는 샘플과 무관하게 생성
    RESULTS_DIR = PROJECT_SUMMARY_DIR  # 프로젝트 전체 요약
    COUNTS_DIR = f"{PROJECT_SUMMARY_DIR}/counts"
    QC_DIR = f"{PROJECT_SUMMARY_DIR}/qc"
    LOGS_DIR = f"{PROJECT_DIR}/logs"
    
    # 표준 구조에서도 legacy 규칙 호환성을 위한 임시 경로
    # (실제로는 샘플별 디렉토리를 사용하지만, expand 등에서 필요)
    TRIMMED_DIR = f"{PROJECT_DIR}/intermediate/trimmed"  # 사용하지 않음
    ALIGNED_DIR = f"{PROJECT_DIR}/intermediate/aligned"  # 사용하지 않음

# 샘플 이름 정의
SAMPLES, = glob_wildcards(f"{RAW_DATA_DIR}/{{sample}}_1.fastq.gz")

# 로그 디렉토리 자동 생성
os.makedirs(f"{LOGS_DIR}/fastqc", exist_ok=True)
os.makedirs(f"{LOGS_DIR}/cutadapt", exist_ok=True)
os.makedirs(f"{LOGS_DIR}/star", exist_ok=True)

# 디버깅: 경로 및 샘플 정보 출력
print(f"=" * 80)
print(f"PIPELINE CONFIGURATION:")
print(f"  Standard Structure: {USE_STANDARD}")
if USE_STANDARD:
    print(f"  Project ID: {PROJECT_ID}")
    print(f"  Pipeline Type: {PIPELINE_TYPE}")
    print(f"  Base Results Dir: {BASE_RESULTS}")
    print(f"  Project Dir: {PROJECT_DIR}")
else:
    print(f"  Legacy Mode - Results Dir: {RESULTS_DIR}")
print(f"  Data Dir: {RAW_DATA_DIR}")
print(f"  Logs Dir: {LOGS_DIR}")
print(f"  Found {len(SAMPLES)} samples")
if SAMPLES:
    print(f"  Sample list: {SAMPLES[:5]}{'...' if len(SAMPLES) > 5 else ''}")
else:
    print(f"  WARNING: No samples found in {RAW_DATA_DIR}")
print(f"=" * 80)

# 경로 변수 (config.yaml에서 로드)
STAR_INDEX = config["star_index"]
ANNOTATION_GTF = config["annotation_gtf"]


# --- 2. 파이프라인의 최종 목표 정의 (Rule all) ---
# QC 리포트 생성 여부에 따라 최종 목표 설정
def get_all_targets():
    """파이프라인 최종 목표 파일 목록 반환"""
    if USE_STANDARD:
        # 표준 구조: 샘플별 manifest + 프로젝트 요약
        targets = []
        for sample in SAMPLES:
            # 각 샘플의 manifest
            targets.append(f"{PROJECT_DIR}/{sample}/rna-seq/final_outputs/manifest.json")
        
        # 프로젝트 전체 요약
        targets.append(f"{COUNTS_DIR}/counts_matrix.txt")
        targets.append(f"{COUNTS_DIR}/counts_matrix_clean.csv")
        
        if config.get("generate_multiqc", True):
            targets.append(f"{RESULTS_DIR}/multiqc_report.html")
        
        return targets
    else:
        # Legacy 구조
        targets = []
        # Raw FastQC
        targets.extend(expand(f"{QC_DIR}/{{sample}}_{{read}}_fastqc.html", sample=SAMPLES, read=[1, 2]))
        
        # FastQC auto-evaluation
        if config.get('fastqc_evaluation', {}).get('enabled', True):
            targets.append(f"{QC_DIR}/{config.get('fastqc_evaluation', {}).get('evaluation_report', 'fastqc_evaluation.txt')}")
        
        # Trimmed reads
        targets.extend(expand(f"{TRIMMED_DIR}/{{sample}}_{{read}}.fastq.gz", sample=SAMPLES, read=[1, 2]))
        
        # Aligned BAM files
        targets.extend(expand(f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES))
        
        # Count matrices
        targets.append(f"{COUNTS_DIR}/counts_matrix.txt")
        targets.append(f"{COUNTS_DIR}/counts_matrix_clean.csv")
        
        # MultiQC report
        if config.get("generate_multiqc", True):
            targets.append(f"{RESULTS_DIR}/multiqc_report.html")
        
        # QC report
        if config.get("generate_qc_report", True):
            targets.append(f"{RESULTS_DIR}/{config.get('qc_report_filename', 'qc_report.html')}")
        
        return targets

rule all:
    input:
        get_all_targets()


rule fastqc_raw:
    input:
        f"{RAW_DATA_DIR}/{{sample}}_{{read}}.fastq.gz"
    output:
        html=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.html",
        zip=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.zip"
    log:
        f"{LOGS_DIR}/fastqc/{{sample}}_{{read}}_raw.log"
    threads: config.get("fastqc_threads", 2)
    shell:
        f"fastqc {{input}} -o {QC_DIR}/ > {{log}} 2>&1"


# --- 3-1. FastQC 결과 자동 평가 (Raw 데이터) ---
rule evaluate_fastqc_raw:
    input:
        fastqc_zips=expand(f"{QC_DIR}/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=[1, 2])
    output:
        report=f"{QC_DIR}/{config.get('fastqc_evaluation', {}).get('evaluation_report', 'fastqc_evaluation.txt')}",
        json=f"{QC_DIR}/{config.get('fastqc_evaluation', {}).get('evaluation_json', 'fastqc_evaluation.json')}"
    params:
        qc_dir=QC_DIR,
        config_params=config.get('fastqc_evaluation', {})
    log:
        f"{LOGS_DIR}/fastqc/evaluate_raw.log"
    shell:
        """
        python3 src/qc/evaluate_fastqc.py \
            {params.qc_dir} \
            -o {output.report} \
            --json {output.json} > {log} 2>&1
        """


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
        bam=f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam",
        log_final=f"{ALIGNED_DIR}/{{sample}}/Log.final.out"
    log:
        f"{LOGS_DIR}/star/{{sample}}.log"
    threads: config["star_threads"]
    resources:
        mem_gb=config.get("star_memory_gb", 35)  # Memory limit per STAR job (GB)
    params:
        out_prefix=lambda wildcards: f"{ALIGNED_DIR}/{wildcards.sample}/"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {STAR_INDEX} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --limitBAMsortRAM {config[star_sort_memory_bytes]} \
             --outBAMsortingThreadN {threads} > {log} 2>&1
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
        python3 src/quantification/convert_counts_matrix.py {input} {output} > {log} 2>&1
        """

# --- 7. MultiQC 리포트 생성 규칙 ---
rule multiqc:
    input:
        # FastQC results
        fastqc_raw=expand(f"{QC_DIR}/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        # Cutadapt logs
        cutadapt_logs=expand(f"{LOGS_DIR}/cutadapt/{{sample}}.log", sample=SAMPLES),
        # STAR logs
        star_logs=expand(f"{ALIGNED_DIR}/{{sample}}/Log.final.out", sample=SAMPLES),
        # featureCounts summary
        fc_summary=f"{COUNTS_DIR}/counts_matrix.txt.summary"
    output:
        f"{RESULTS_DIR}/multiqc_report.html"
    params:
        results_dir=RESULTS_DIR,
        qc_dir=QC_DIR,
        logs_dir=LOGS_DIR,
        aligned_dir=ALIGNED_DIR,
        counts_dir=COUNTS_DIR
    log:
        f"{LOGS_DIR}/multiqc.log"
    shell:
        """
        multiqc {params.qc_dir} {params.logs_dir} {params.aligned_dir} {params.counts_dir} \
            -o {params.results_dir} \
            -n multiqc_report.html \
            --force \
            --title "RNA-seq QC Report" \
            --comment "Comprehensive quality control report for RNA-seq analysis" \
            > {log} 2>&1
        """

# --- 8. QC 리포트 생성 규칙 ---
rule generate_qc_report:
    input:
        counts=f"{COUNTS_DIR}/counts_matrix.txt",
        counts_summary=f"{COUNTS_DIR}/counts_matrix.txt.summary",
        cutadapt_logs=expand(f"{LOGS_DIR}/cutadapt/{{sample}}.log", sample=SAMPLES),
        star_logs=expand(f"{ALIGNED_DIR}/{{sample}}/Log.final.out", sample=SAMPLES)
    output:
        f"{RESULTS_DIR}/{config.get('qc_report_filename', 'qc_report.html')}"
    params:
        top_genes=config.get("qc_top_genes", 10),
        logs_dir=LOGS_DIR,
        results_dir=RESULTS_DIR,
        data_dir=DATA_DIR
    log:
        f"{LOGS_DIR}/qc_report.log"
    script:
        "src/qc/generate_qc_report.py"


# ========================================================================
# Phase 3: Standard Structure Rules
# ========================================================================

# --- 9. BAM Index 생성 (samtools index) ---
rule index_bam:
    input:
        bam=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam" if USE_STANDARD else f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        bai=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam.bai" if USE_STANDARD else f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"
    log:
        f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/samtools_index.log" if USE_STANDARD else f"{LOGS_DIR}/samtools/{{sample}}_index.log"
    conda:
        "environment.yaml"
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """


# --- 10. QC Summary 생성 (표준 구조 전용) ---
rule generate_qc_summary:
    input:
        star_log=f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/star_final.log" if USE_STANDARD else f"{ALIGNED_DIR}/{{sample}}/Log.final.out",
        fc_summary=f"{COUNTS_DIR}/counts_matrix.txt.summary",
        # 표준 구조에서는 BAM 복사를 기다림
        bam=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam" if USE_STANDARD else []
    output:
        qc_json=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/qc/qc_summary.json" if USE_STANDARD else f"{QC_DIR}/{{sample}}_qc_summary.json"
    params:
        sample_id="{sample}"
    log:
        f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/qc_summary.log" if USE_STANDARD else f"{LOGS_DIR}/qc/{{sample}}_summary.log"
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p $(dirname {output.qc_json})
        mkdir -p $(dirname {log})
        python3 scripts/standardization/generate_qc_summary.py \
            --sample-id {params.sample_id} \
            --star-log {input.star_log} \
            --featurecounts {input.fc_summary} \
            -o {output.qc_json} > {log} 2>&1
        """


# --- 11. Manifest 생성 (표준 구조 전용) ---
rule generate_manifest:
    input:
        bam=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam" if USE_STANDARD else f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam",
        bai=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam.bai" if USE_STANDARD else f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",
        qc_summary=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/qc/qc_summary.json" if USE_STANDARD else f"{QC_DIR}/{{sample}}_qc_summary.json"
    output:
        manifest=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/manifest.json" if USE_STANDARD else f"{RESULTS_DIR}/{{sample}}_manifest.json"
    params:
        sample_id="{sample}",
        sample_dir=f"{PROJECT_DIR}/{{sample}}/rna-seq" if USE_STANDARD else RESULTS_DIR,
        project_id=PROJECT_ID,
        pipeline_type=PIPELINE_TYPE
    log:
        f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/manifest.log" if USE_STANDARD else f"{LOGS_DIR}/manifest/{{sample}}.log"
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p $(dirname {output.manifest})
        mkdir -p $(dirname {log})
        python3 scripts/standardization/generate_manifest.py \
            --sample-dir {params.sample_dir} \
            --sample-id {params.sample_id} \
            --project-id {params.project_id} \
            --pipeline-type {params.pipeline_type} > {log} 2>&1
        """


# --- 12. 표준 구조로 BAM 복사 (star_align 후처리) ---
if USE_STANDARD:
    rule copy_bam_to_standard:
        input:
            bam=f"{ALIGNED_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam",
            log_final=f"{ALIGNED_DIR}/{{sample}}/Log.final.out"
        output:
            bam=f"{PROJECT_DIR}/{{sample}}/rna-seq/final_outputs/bam/aligned.sorted.bam",
            log_final=f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/star_final.log"
        log:
            f"{PROJECT_DIR}/{{sample}}/rna-seq/intermediate/logs/copy_bam.log"
        shell:
            """
            mkdir -p $(dirname {output.bam})
            mkdir -p $(dirname {output.log_final})
            mkdir -p $(dirname {log})
            cp {input.bam} {output.bam} > {log} 2>&1
            cp {input.log_final} {output.log_final} >> {log} 2>&1
            """
