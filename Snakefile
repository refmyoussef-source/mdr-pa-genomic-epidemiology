# ============================================
# Snakefile: Project 5 MDR-PA Pipeline
# VERSION 3.1 (Includes Download + QC + Trimming + Post-QC) # <-- V3.1
# ============================================

# --- 1. Configuration (Load metadata) ---
import pandas as pd
metadata = pd.read_csv("results/metadata/metadata_clean.csv")
SAMPLES = metadata["sample_id"].tolist()
SAMPLE_TO_RUN = metadata.set_index('sample_id')['run_id'].to_dict()

# ============================================
# --- 2. Define All Final Files (Our "Goal" Rule) ---
# ============================================

# --- Our main goal is still the set of 192 TRIMMED reads ---
def get_trimmed_reads(wildcards):
    return expand(
        "results/trimmed_reads/{sample_id}_{read_pair}.fastq.gz",
        sample_id=SAMPLES,
        read_pair=["1", "2"]
    )

rule all:
    input:
        # Our new goal is the list of all trimmed fastq files
        get_trimmed_reads
        
        # --- (NEW V3.1) We ALSO add our new reports as "goals" ---
        # This way, a simple "snakemake" will build them if missing,
        # but "rule all" is still primarily about the reads.
        # This is optional, but good practice.
        # "results/qc/multiqc_report_fastp.html" # You can add this if you want

# ============================================
# --- 3. "Worker" Rules (How to build things) ---
# ============================================

# --- Rule to run fastp Trimming ---
rule fastp_trimming:
    input:
        # It needs the RAW reads to start
        r1 = "data/raw_reads/{sample_id}_1.fastq",
        r2 = "data/raw_reads/{sample_id}_2.fastq"
    output:
        # It produces CLEANED (trimmed) reads
        r1_trimmed = "results/trimmed_reads/{sample_id}_1.fastq.gz",
        r2_trimmed = "results/trimmed_reads/{sample_id}_2.fastq.gz",
        # It also produces its own QC reports
        html = "results/qc/fastp/{sample_id}.html",
        json = "results/qc/fastp/{sample_id}.json"
    params:
        out_dir = "results/qc/fastp" # Directory for reports
    shell:
        """
        echo "--- (STEP C) Running fastp Trimming on {wildcards.sample_id} ---"
        # Ensure the output directories exist
        mkdir -p results/trimmed_reads
        mkdir -p {params.out_dir}
        
        # Run fastp
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1_trimmed} \
            -O {output.r2_trimmed} \
            -h {output.html} \
            -j {output.json} \
            --detect_adapter_for_pe # This is the "Diagnosis" we saw!
        """

# --- (NEW RULE V3.1) Rule to run MultiQC on FASTP results ---
rule multiqc_fastp:
    input:
        # We need all 96 JSON reports from fastp
        # (MultiQC prefers the .json files from fastp)
        expand("results/qc/fastp/{sample_id}.json", sample_id=SAMPLES)
    output:
        # A new, clearly named master report
        html = "results/qc/multiqc_report_fastp.html"
    params:
        report_dir = "results/qc/fastp",  # The directory where fastp reports live
        out_dir = "results/qc"            # Where to put the final multiqc report
    shell:
        """
        echo "--- (STEP D) Aggregating 96 fastp reports with MultiQC ---"
        multiqc {params.report_dir} \
            --outdir {params.out_dir} \
            --filename multiqc_report_fastp.html \
            --title "Project 5: Post-Trimming QC (fastp)"
        """

# --- (OLD RULE) Rule to run MultiQC (on FastQC results) ---
# (We leave this here, it's our "Pre-QC" report)
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html", sample_id=SAMPLES, read_pair=["1", "2"])
    output:
        "results/qc/multiqc_report.html"
    params:
        qc_dir = "results/qc/fastqc",
        out_dir = "results/qc"
    shell:
        "multiqc {params.qc_dir} --outdir {params.out_dir} --filename multiqc_report.html" # Added --filename for safety

# --- (OLD RULE) Rule to run FastQC ---
rule fastqc:
    input:
        fastq = "data/raw_reads/{sample_id}_{read_pair}.fastq"
    output:
        html = "results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html",
        zip = "results/qc/fastqc/{sample_id}_{read_pair}_fastqc.zip"
    params:
        out_dir = "results/qc/fastqc"
    shell:
        "fastqc {input.fastq} -o {params.out_dir}"

# --- (OLD RULE) Rule to download data ---
rule download_sra_data:
    output:
        r1 = "data/raw_reads/{sample_id}_1.fastq",
        r2 = "data/raw_reads/{sample_id}_2.fastq"
    params:
        run_id = lambda wildcards: SAMPLE_TO_RUN[wildcards.sample_id]
    shell:
        """
        echo "--- (STEP A) Downloading SRR ID: {params.run_id} for Sample: {wildcards.sample_id} ---"
        fasterq-dump --split-files -O data/raw_reads -p {params.run_id}
        
        echo "--- Renaming {params.run_id} to {wildcards.sample_id} ---"
        mv data/raw_reads/{params.run_id}_1.fastq {output.r1}
        mv data/raw_reads/{params.run_id}_2.fastq {output.r2}
        """