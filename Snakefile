# ============================================
# Snakefile: Project 5 MDR-PA Pipeline
# (Phase 7: Merging & Filtering Automation)
# ============================================

# --- 1. Configuration (Load metadata) ---
import pandas as pd

# --- [V5.0] We manage TWO lists (Rule 4) ---

# The "Legacy" list of all 96 samples (for QC rules)
metadata_all = pd.read_csv("results/metadata/metadata_clean.csv")
ALL_SAMPLES = metadata_all["sample_id"].tolist()
SAMPLE_TO_RUN = metadata_all.set_index('sample_id')['run_id'].to_dict()

# The "Clean" list of 93 samples (for Phase 6 and beyond)
metadata_clean = pd.read_csv("results/metadata/metadata_final_cohort.csv")
CLEAN_SAMPLES = metadata_clean["sample_id"].tolist()

# ---  Define Reference Genome ---
REF_GENOME = "data/reference_genome/PAO1_reference.fna"

# ============================================
# --- 2. Define All Final Files (Our "Goal" Rule) ---
# ============================================

# ---  Our NEW (and FINAL) goal is the ONE Analysis-Ready VCF ---
# (This is the file for the Tree)
rule all:
    input:
        "results/variant_calling/ANALYSIS_READY.vcf.gz"

# ============================================
# --- 3. "Worker" Rules (How to build things) ---
# (Newest rules go on top)
# ============================================

# --- [NEW RULE V6.0] - (From R&D Notebook 06, Cell 6) ---
# This filters the Master VCF for high-quality SNPs.
rule filter_vcf:
    input:
        master_vcf = "results/variant_calling/MASTER_MERGED.vcf.gz"
    output:
        analysis_vcf = "results/variant_calling/ANALYSIS_READY.vcf.gz"
    shell:
        """
        echo "--- (STEP H.2) Filtering Master VCF for high-quality SNPs ---"
        
        # This is the "Filter Recipe" from Notebook 06
        bcftools view -i 'TYPE=="snp" & QUAL > 30' \
            -O z -o {output.analysis_vcf} {input.master_vcf}
        """

# --- [NEW RULE V6.0] - (From R&D Notebook 06, Cell 4) ---
# This merges all 93 VCFs into one Master VCF.
rule merge_vcf:
    input:
        # 1. The list of VCFs (This is a "source file" from Notebook 06/Git)
        vcf_list = "results/variant_calling/vcf_list_for_merge.txt",
        # 2. All 93 VCFs (This forces Snakemake to run Phase 6 first)
        vcf_files = expand("results/variant_calling/{sample_id}.vcf.gz", sample_id=CLEAN_SAMPLES),
        vcf_indices = expand("results/variant_calling/{sample_id}.vcf.gz.csi", sample_id=CLEAN_SAMPLES)
    output:
        master_vcf = "results/variant_calling/MASTER_MERGED.vcf.gz"
    shell:
        """
        echo "--- (STEP H.1) Merging 93 VCF files into Master VCF ---"
        
        # This is the "Merge Recipe" from Notebook 06
        # Note: The paths in vcf_list_for_merge.txt are already Root-Relative
        bcftools merge -l {input.vcf_list} -O z -o {output.master_vcf}
        """

# --- [RULE V5.0] - (From R&D Notebook 05) ---
# This is the full Haploid Variant Calling pipe.
rule bcftools_call_haploid:
    input:
        ref = REF_GENOME,
        # It depends on the sorted AND indexed BAM file
        bam = "results/mapped_reads/{sample_id}.sorted.bam",
        bai = "results/mapped_reads/{sample_id}.sorted.bam.bai"
    output:
        # We create the VCF file...
        vcf = "results/variant_calling/{sample_id}.vcf.gz",
        # ...and its index (.csi)
        idx = "results/variant_calling/{sample_id}.vcf.gz.csi"
    shell:
        """
        echo "--- (STEP G.1) Haploid Variant Calling for {wildcards.sample_id} ---"
        
        # This is the exact "Haploid Recipe" from Notebook 05 (Cell 6)
        bcftools mpileup -f {input.ref} {input.bam} \
            | bcftools call --ploidy 1 -mv -O z -o {output.vcf} -
        
        echo "--- (STEP G.2) Indexing VCF for {wildcards.sample_id} ---"
        bcftools index {output.vcf}
        """
        
# ============================================
# --- 4. "Legacy" Rules (Phase 1-5 QC) ---
# (We keep them here for reproducibility)
# ============================================

# --- (Rule V4.3) Aggregate samtools stats with MultiQC ---
rule multiqc_samtools_stats:
    input:
        expand("results/qc/mapping_stats/{sample_id}.stats", sample_id=ALL_SAMPLES)
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

# --- (Rule V4.3) Generate samtools stats for each BAM ---
rule samtools_stats:
    input:
        bam = "results/mapped_reads/{sample_id}.sorted.bam"
    output:
        stats = "results/qc/mapping_stats/{sample_id}.stats"
    shell:
        """
        echo "--- (STEP F.1) Generating Mapping Stats for {wildcards.sample_id} ---"
        mkdir -p results/qc/mapping_stats
        samtools stats {input.bam} > {output.stats}
        """

# --- (Rule V4.2) Creates the .bai index for each sorted BAM
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

# --- (Rule V4.2) This is the main "Pipe" rule
rule bwa_map_and_sort:
    input:
        r1 = "results/trimmed_reads/{sample_id}_1.fastq.gz",
        r2 = "results/trimmed_reads/{sample_id}_2.fastq.gz",
        ref_fasta = REF_GENOME,
        ref_index = expand(f"{REF_GENOME}.{{ext}}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    output:
        bam = "results/mapped_reads/{sample_id}.sorted.bam"
    params:
        read_group = f"@RG\\tID:{{sample_id}}\\tSM:{{sample_id}}\\tPL:ILLUMINA"
    shell:
        """
        echo "--- (STEP E.1) Mapping & Sorting {wildcards.sample_id} ---"
        bwa-mem2 mem -R '{params.read_group}' {input.ref_fasta} {input.r1} {input.r2} \
            | samtools view -bS - \
            | samtools sort -o {output.bam} -
        """

# ---  This rule creates the BWA index.
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

# --- (Rule V3.1) Rule to run MultiQC on FASTP results ---
rule multiqc_fastp:
    input:
        # (Using the FULL 96 list: ALL_SAMPLES)
        expand("results/qc/fastp/{sample_id}.json", sample_id=ALL_SAMPLES)
    output:
        html = "results/qc/multiqc_report_fastp.html"
    params:
        report_dir = "results/qc/fastp",
        out_dir = "results/qc"
    shell:
        """
        echo "--- (STEP D) Aggregating 96 fastp reports with MultiQC ---"
        multiqc {params.report_dir} \
            --outdir {params.out_dir} \
            --filename multiqc_report_fastp.html \
            --title "Project 5: Post-Trimming QC (fastp)"
        """

# ---  Rule to run fastp Trimming ---
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
        echo "--- (STEP C) Running fastp Trimming on {wildcards.sample_id} ---"
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
        # (Using the FULL 96 list: ALL_SAMPLES)
        expand("results/qc/fastqc/{sample_id}_{read_pair}_fastqc.html", sample_id=ALL_SAMPLES, read_pair=["1", "2"])
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
        # We use the full lookup dictionary (96 samples)
        run_id = lambda wildcards: SAMPLE_TO_RUN[wildcards.sample_id]
    shell:
        """
        echo "--- (STEP A) Downloading SRR ID: {params.run_id} for Sample: {wildcards.sample_id} ---"
        fasterq-dump --split-files -O data/raw_reads -p {params.run_id}
        
        echo "--- Renaming {params.run_id} to {wildcards.sample_id} ---"
        mv data/raw_reads/{params.run_id}_1.fastq {output.r1}
        mv data/raw_reads/{params.run_id}_2.fastq {output.r2}
        """
