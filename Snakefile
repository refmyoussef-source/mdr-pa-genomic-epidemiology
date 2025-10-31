# ============================================
# Snakefile: Project 5 MDR-PA Pipeline
# VERSION 2.1 (Includes Download + QC - Syntax Fixed)
# ============================================

# --- 1. Configuration (Load metadata) ---
import pandas as pd
metadata = pd.read_csv("results/metadata/metadata_clean.csv")
SAMPLES = metadata["sample_id"].tolist()
SAMPLE_TO_RUN = metadata.set_index('sample_id')['run_id'].to_dict()

# ============================================
# --- 2. Define All Final Files (Our "Goal" Rule) ---
# ============================================

# This is the "target" rule. We tell Snakemake what we want in the end.
# Our NEW final goal is the MultiQC report
FINAL_QC_REPORT = "results/qc/multiqc_report.html"

# This helper function generates the list of 192 HTML files (96 samples * 2 reads)
def get_fastqc_reports(wildcards):
    return expand(
        "results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html",
        sample_id=SAMPLES,
        read_pair=["1", "2"]
    )

rule all:
    input:
        # Our new goal is the final MultiQC report
        FINAL_QC_REPORT

# ============================================
# --- 3. "Worker" Rules (How to build things) ---
# ============================================

# --- (Rule 3a) Rule to run MultiQC ---
# This rule collects all FastQC reports into one summary
rule multiqc:
    input:
        get_fastqc_reports # Wait for all 192 fastqc reports
    output:
        FINAL_QC_REPORT     
    params:
        # The folder where the reports are
        qc_dir = "results/qc/fastqc",
        # Where to put the multiqc report
        out_dir = "results/qc"
    shell:
        """
        echo "--- (STEP C) Aggregating all FastQC reports with MultiQC ---"
        multiqc {params.qc_dir} --outdir {params.out_dir}
        """

# --- (Rule 3b) Rule to run FastQC ---
# This rule runs FastQC on *one* FASTQ file
rule fastqc:
    input:
        # {read_pair} will be "1" or "2"
        fastq = "data/raw_reads/{sample_id}_{read_pair}.fastq"
    output:
        # It creates an HTML file
        html = "results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html",
        # And a zip file
        zip = "results/qc/fastqc/{sample_id}_{read_pair}_fastqc.zip"
    params:
        # The output directory
        out_dir = "results/qc/fastqc"
    shell:
        """
        echo "--- (STEP B) Running FastQC on {input.fastq} ---"
        # -o specifies the output directory
        fastqc {input.fastq} -o {params.out_dir}
        """

# --- (Rule 3c) Rule to download data ---
# (This rule is the same, but now FastQC *depends* on it)
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