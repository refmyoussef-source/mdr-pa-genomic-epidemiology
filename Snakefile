# ============================================
# Snakefile: Project 5 MDR-PA Pipeline
# ============================================

# --- 1. Configuration (Load metadata) ---
# We load our clean sample map here
import pandas as pd
metadata = pd.read_csv("results/metadata/metadata_clean.csv")

# Get a list of all sample IDs (e.g., PA097, PA096, ... HOU1)
SAMPLES = metadata["sample_id"].tolist()

# Get a list of all run IDs (e.g., SRR16632095, SRR16632096, ...)
RUN_IDS = metadata["run_id"].tolist()

# Create the critical "dictionary" to map Sample ID to Run ID
# { "PA097": "SRR16632095", "PA096": "SRR16632096", ... }
SAMPLE_TO_RUN = metadata.set_index('sample_id')['run_id'].to_dict()


# --- 2. Define All Final Files (Our "Goal" Rule) ---
# This is the "target" rule. We tell Snakemake what we want in the end.
# We want all 96 raw read files (R1 and R2 for each sample).

rule all:
    input:
        # expand() is a Snakemake function that builds the full list of files
        # It takes our list 'SAMPLES' and puts each sample_id into the string
        expand("data/raw_reads/{sample_id}_1.fastq", sample_id=SAMPLES),
        expand("data/raw_reads/{sample_id}_2.fastq", sample_id=SAMPLES)


# --- 3. Define How to Get The Files (The "Worker" Rule) ---
# This rule teaches Snakemake how to create the FASTQ files
# if they are missing.

rule download_sra_data:
    output:
        # The files this rule *creates*
        # {sample_id} is a "wildcard" that Snakemake fills in
        r1 = "data/raw_reads/{sample_id}_1.fastq",
        r2 = "data/raw_reads/{sample_id}_2.fastq"
    
    params:
        # Get the correct Run ID (e.g., SRR...) for the given Sample ID (e.g., PA097)
        run_id = lambda wildcards: SAMPLE_TO_RUN[wildcards.sample_id]
    
    shell:
        # This is the corrected shell block, enclosed in one """ string
        """
        echo "--- (STEP A) Downloading SRR ID: {params.run_id} for Sample: {wildcards.sample_id} ---"
        # This creates files named using the Run ID (e.g., SRR..._1.fastq)
        fasterq-dump --split-files -O data/raw_reads -p {params.run_id}
        
        echo "--- (STEP B) Renaming {params.run_id} to {wildcards.sample_id} ---"
        # This is critical for our pipeline's organization.
        mv data/raw_reads/{params.run_id}_1.fastq {output.r1}
        mv data/raw_reads/{params.run_id}_2.fastq {output.r2}
        """