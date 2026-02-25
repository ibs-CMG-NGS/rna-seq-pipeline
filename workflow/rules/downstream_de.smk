# ============================================================================
# RNA-seq Pipeline → DE/GO Analysis Integration
# ============================================================================
# 
# Add this to the end of your main Snakefile to automatically trigger
# DE/GO analysis after RNA-seq preprocessing completes.
#
# Usage:
#   snakemake --cores 12 --until trigger_de_analysis
#

# Configuration
DE_PIPELINE_DIR = config.get("de_pipeline_dir", "../RNA-Seq_DE_GO_analysis")
DE_CONFIG_FILE = config.get("de_config_file", "config_mouse_chd8.yml")

# ============================================================================
# Rule: Bridge to DE/GO Pipeline
# ============================================================================

rule bridge_to_de_pipeline:
    """
    Prepare counts matrix and metadata for DE/GO analysis pipeline.
    """
    input:
        counts = f"{PROJECT_RESULTS_DIR}/counts/counts_matrix.txt",
        summary = f"{PROJECT_RESULTS_DIR}/project_summary.json",
        sample_sheet = f"config/samples/{PROJECT_ID}.tsv"
    output:
        de_counts = f"{DE_PIPELINE_DIR}/data/raw/{PROJECT_ID}_counts.csv",
        bridge_log = f"{PROJECT_RESULTS_DIR}/bridge_to_de.log"
    log:
        f"{PROJECT_RESULTS_DIR}/logs/bridge_to_de.log"
    shell:
        """
        echo "Bridging RNA-seq results to DE/GO pipeline..." > {output.bridge_log}
        
        # Copy counts matrix
        cp {input.counts} {output.de_counts}
        echo "✅ Counts matrix copied" >> {output.bridge_log}
        
        # Check QC status
        python3 scripts/standardization/agent_query.py \
            --project-summary {input.summary} \
            --query status >> {output.bridge_log}
        
        echo "✅ Bridge preparation complete" >> {output.bridge_log}
        echo "Next: Configure {DE_PIPELINE_DIR}/config.yml and run DE analysis" >> {output.bridge_log}
        """

rule trigger_de_analysis:
    """
    Automatically trigger DE/GO analysis pipeline (optional).
    Only use if DE pipeline is ready with proper config.
    """
    input:
        de_counts = rules.bridge_to_de_pipeline.output.de_counts
    output:
        de_results = f"{DE_PIPELINE_DIR}/output/{PROJECT_ID}/final_de_results.csv"
    log:
        f"{PROJECT_RESULTS_DIR}/logs/de_analysis.log"
    params:
        de_dir = DE_PIPELINE_DIR,
        config_file = DE_CONFIG_FILE
    conda:
        "snakemake_env"  # Snakemake execution environment
    shell:
        """
        echo "Triggering DE/GO analysis pipeline..."
        
        cd {params.de_dir}
        
        # Run DE/GO pipeline
        snakemake \
            --configfile {params.config_file} \
            --cores 4 \
            --use-conda \
            2>&1 | tee {log}
        
        echo "✅ DE/GO analysis complete"
        """

# ============================================================================
# Rule: All (extended with DE/GO)
# ============================================================================

rule all_with_de:
    """
    Complete RNA-seq preprocessing + DE/GO analysis.
    
    Usage:
        snakemake --cores 12 all_with_de
    """
    input:
        # RNA-seq preprocessing outputs
        expand(
            f"{PROJECT_RESULTS_DIR}/{{sample}}/rna-seq/final_outputs/manifest.json",
            sample=SAMPLES
        ),
        f"{PROJECT_RESULTS_DIR}/project_summary.json",
        
        # DE/GO analysis outputs
        f"{DE_PIPELINE_DIR}/output/{PROJECT_ID}/final_de_results.csv"
    message:
        "✅ Complete pipeline (RNA-seq → DE/GO) finished!"
