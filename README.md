# Project 5: High-Throughput Genomic Epidemiology of 96 MDR P. aeruginosa Strains

**Goal:** This project implements a fully automated and reproducible Snakemake pipeline (V6.0) to analyze the 96 MDR *P. aeruginosa* strains from BioProject PRJNA771342. The primary scientific goal is to reproduce the phylogenetic tree and identify known epidemic clusters (e.g., ST2855).

---

## Final Phylogenetic Tree (The "Showroom")

The final analysis-ready tree (built from 93 high-quality samples, excluding 3 outliers) is visualized below. This tree was generated using `bcftools` (Haploid Calling) and `IQ-TREE` (Standard Tree Search).

![Final Phylogenetic Tree](tree_graph.png)

---

## The Automated Pipeline (The "Scientific Story")

This pipeline is managed by `Snakemake` and ensures 100% reproducibility. The "R&D" (testing) was performed in Jupyter Notebooks, and the "Production" (heavy lifting) was automated in the `Snakefile`.

### Phase 1: Download (SRA -> FASTQ)
* **Action:** `rule download_sra_data`
* **Result:** Downloaded 113GB of raw `.fastq` data for 96 samples.

### Phase 2: Pre-QC (Diagnosis)
* **Action:** `rule fastqc` & `rule multiqc`
* **Result:** **Critical Discovery:** The `multiqc_report.html` (documented in `Notebook 01`) revealed significant Adapter Contamination.

### Phase 3: Trimming (The "Cure")
* **Action:** `rule fastp_trimming`
* **Result:** Cleaned 113GB of data using `fastp` to remove adapters.

### Phase 4: Post-QC (Verification)
* **Action:** `rule multiqc_fastp`
* **Result:** **Verification:** The `multiqc_report_fastp.html` (documented in `Notebook 02`) confirmed that all adapter contamination was successfully removed.

### Phase 5: Mapping & QC (The "Outlier" Decision)
* **Action:** `rule bwa_map_and_sort` & `rule samtools_stats`
* **Result:** **Critical Decision:** The `multiqc_mapping_stats.html` (documented in `Notebook 04`) revealed 3 low-quality outlier samples (`PA033`, `PA045`, `PA070`).
* **Action:** We created `metadata_final_cohort.csv` (documented in `Notebook 05`) to **exclude** these 3 samples, leaving a "clean cohort" of 93 samples.

### Phase 6: Variant Calling (The "Haploid" Fix)
* **Action:** `rule bcftools_call_haploid`
* **Result:** **Critical Discovery:** The R&D test (documented in `Notebook 05`) revealed the default "diploid" assumption was scientifically wrong.
* **Action:** The pipeline was fixed to use `--ploidy 1` (Haploid) for bacterial genomics, generating 93 clean `.vcf.gz` files.

### Phase 7: Merging & Filtering (The "Golden File")
* **Action:** `rule merge_vcf` & `rule filter_vcf`
* **Result:** Created the final `ANALYSIS_READY.vcf.gz` (15M) file, containing only high-confidence SNPs (`QUAL > 30`) for our 93 clean samples.

### Phase 8: Phylogenetics (The Final Goal)
* **Action:** `vcf2phylip.py` (local script) & `iqtree` (R&D in `Notebook 07`)
* **Result:** Converted the final VCF to `ANALYSIS_READY.min4.phy` (11M) and generated the final phylogenetic tree (`.treefile`) visualized above.

---

## How to Run This Project (Reproducibility)

1.  **Clone the repo:**
    ```bash
    git clone https://github.com/refmyoussef-source/mdr-pa-genomic-epidemiology.git
    cd mdr-pa-genomic-epidemiology
    ```
2.  **Create the environment (V5.3):**
    ```bash
    mamba env create -f environment.yml
    conda activate popgen_env
    ```
3.  **Run the entire pipeline (The "Factory"):**
    ```bash
    # This will run all steps (Download, Trim, Map, Call, Merge, Filter)
    # WARNING: This will take ~2 days and ~150GB of disk space.
    snakemake --cores 8
    ```
