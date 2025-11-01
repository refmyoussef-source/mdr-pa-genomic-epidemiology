# ============================================
# Snakefile: Project 5 MDR-PA Pipeline
# VERSION 4.3 (Includes Mapping Stats QC)
# ============================================

# --- 1. Configuration (Load metadata) ---
import pandas as pd
metadata = pd.read_csv("results/metadata/metadata_clean.csv")
SAMPLES = metadata["sample_id"].tolist()
SAMPLE_TO_RUN = metadata.set_index('sample_id')['run_id'].to_dict()

# --- [V4.2] Define Reference Genome ---
REF_GENOME = "data/reference_genome/PAO1_reference.fna"

# ============================================
# --- 2. Define All Final Files (Our "Goal" Rule) ---
# ============================================

# --- [V4.2] Our main goal is the set of 96 INDEXED BAM files ---
def get_bam_indices(wildcards):
    return expand(
        "results/mapped_reads/{sample_id}.sorted.bam.bai",
        sample_id=SAMPLES
    )

rule all:
    input:
        get_bam_indices,
        # --- [NEW V4.3] We ALSO add our new Mapping QC report as a "goal" ---
        "results/qc/mapping_stats/multiqc_mapping_stats.html"

# ============================================
# --- 3. "Worker" Rules (How to build things) ---
# ============================================

# --- [NEW RULE V4.3] Aggregate samtools stats with MultiQC ---
# (This is Step F.2 - It feeds Notebook 04)
rule multiqc_samtools_stats:
    input:
        expand("results/qc/mapping_stats/{sample_id}.stats", sample_id=SAMPLES)
    output:
        "results/qc/mapping_stats/multiqc_mapping_stats.html"
    params:
        stats_dir = "results/qc/mapping_stats",
        out_dir = "results/qc/mapping_stats"
    shell:
        """
        echo "--- (STEP F.2) Aggregating 96 Mapping Stats reports ---"
        multiqc {params.stats_dir} \
            --outdir {params.out_dir} \
            --filename multiqc_mapping_stats.html \
            --title "Project 5: Mapping QC (Samtools Stats)"
        """

# --- [NEW RULE V4.3] Generate samtools stats for each BAM ---
# (This is Step F.1)
rule samtools_stats:
    input:
        bam = "results/mapped_reads/{sample_id}.sorted.bam"
    output:
        stats = "results/qc/mapping_stats/{sample_id}.stats"
    shell:
        """
        echo "--- (STEP F.1) Generating Mapping Stats for {wildcards.sample_id} ---"
        # Create the directory first (Snakemake doesn't do this automatically for shell)
        mkdir -p results/qc/mapping_stats
        
        # 'stats' command generates a text-based report
        samtools stats {input.bam} > {output.stats}
        """

# --- [RULE V4.2] - (From R&D Test 4) ---
# Creates the .bai index for each sorted BAM
rule samtools_index:
    input:
        bam = "results/mapped_reads/{sample_id}.sorted.bam"
    output:
        bai = "results/mapped_reads/{sample_id}.sorted.bam.bai"
    shell:
        """
        echo "--- (STEP E.2) Indexing BAM for {wildcards.sample_id} ---"
        samtools index {input.bam}
        """

# --- [RULE V4.2] - (From R&D Test 2+3) ---
# This is the main "Pipe" rule
rule bwa_map_and_sort:
    input:
        # It needs the clean reads
        r1 = "results/trimmed_reads/{sample_id}_1.fastq.gz",
        r2 = "results/trimmed_reads/{sample_id}_2.fastq.gz",
        # And it needs the *indexed* reference genome (from rule bwa_index)
        ref_fasta = REF_GENOME,
        ref_index = expand(f"{REF_GENOME}.{{ext}}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    output:
        # It produces one sorted BAM file
        bam = "results/mapped_reads/{sample_id}.sorted.bam"
    params:
        # This is the Read Group from Notebook 03
        read_group = f"@RG\\tID:{{sample_id}}\\tSM:{{sample_id}}\\tPL:ILLUMINA"
    shell:
        """
        echo "--- (STEP E.1) Mapping & Sorting {wildcards.sample_id} ---"
        
        # This is the exact "Pipe" command from Notebook 03 (Cell 7.2)
        bwa-mem2 mem -R '{params.read_group}' {input.ref_fasta} {input.r1} {input.r2} \
            | samtools view -bS - \
            | samtools sort -o {output.bam} -
        """

# --- [RULE V4.2] - (From R&D Test 1) ---
# This rule creates the BWA index.
rule bwa_index:
    input:
        REF_GENOME
    output:
        expand(f"{REF_GENOME}.{{ext}}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    shell:
        """
        echo "--- (STEP E.0) Creating BWA-MEM2 Index for Reference Genome ---"
        bwa-mem2 index {input}
        """

# ============================================
# --- 4. "Legacy" Rules (Phase 1-4) ---
# ============================================

# --- (Rule V3.1) Rule to run MultiQC on FASTP results ---
rule multiqc_fastp:
    input:
        expand("results/qc/fastp/{sample_id}.json", sample_id=SAMPLES)
    output:
        html = "results/qc/multiqc_report_fastp.html"
    params:
        report_dir = "results/qc/fastp",
        out_dir = "results/qc"
    shell:
        """
        multiqc {params.report_dir} \
            --outdir {params.out_dir} \
            --filename multiqc_report_fastp.html \
            --title "Project 5: Post-Trimming QC (fastp)"
        """

# --- (Rule V3.0) Rule to run fastp Trimming ---
rule fastp_trimming:
    input:
        r1 = "data/raw_reads/{sample_id}_1.fastq",
        r2 = "data/raw_reads/{sample_id}_2.fastq"
    output:
        r1_trimmed = "results/trimmed_reads/{sample_id}_1.fastq.gz",
        r2_trimmed = "results/trimmed_reads/{sample_id}_2.fastq.gz",
        html = "results/qc/fastp/{sample_id}.html",
        json = "results/qc/fastp/{sample_id}.json"
    shell:
        """
        mkdir -p results/trimmed_reads
        mkdir -p results/qc/fastp
        
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1_trimmed} \
            -O {output.r2_trimmed} \
            -h {output.html} \
            -j {output.json} \
            --detect_adapter_for_pe
        """

# --- (OLD RULE) Rule to run MultiQC (on FastQC results) ---
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html", sample_id=SAMPLES, read_pair=["1", "2"])
    output:
        "results/qc/multiqc_report.html"
    params:
        qc_dir = "results/qc/fastqc",
        out_dir = "results/qc"
    shell:
        "multiqc {params.qc_dir} --outdir {params.out_dir} --filename multiqc_report.html"

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
        fasterq-dump --split-files -O data/raw_reads -p {params.run_id}
        mv data/raw_reads/{params.run_id}_1.fastq {output.r1}
        mv data/raw_reads/{params.run_id}_2.fastq {output.r2}
        """