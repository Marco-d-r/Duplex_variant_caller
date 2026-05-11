# 02 SNV Calling Pipeline

### Overview
This is the backbone of the variant caller. There are 4 distinct steps to this pipeline:
1. Raw pileup and feature extraction
2. Annotation with GnomAD and panel of normals (PoN)
3. Threshold tuning
4. Application of thresholds and final VCF production

> **Execution Note:** All of the scripts here are meant to be run **locally**, after downloading the consensus `DCS.bam` and `SSCS.bam` from the HPC. 

### Environment & Resources
* **Environment:** The environment that supports all this section is `02_snv_calling_environment.yaml`.
* **Reference Genome:** GRCh38
* **Blacklist:** `/resources/master_blacklist.bed`
* **Panel BED File:** `/resources/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed`
* **GnomAD:** `/resources/af-only-gnomad.hg38_common_variants.vcf.gz` (and index)
* **Panel of Normals (PoN):** `/resources/snv_cb_pon_files/`

---

## 1) Raw Pileup and Feature Extraction (`02.1_snv_extraction/`)

The aim of this section is to make the 'rawest' possible pileup from the `DCS.bam`, that is, identify every difference between the reference and the reads. Whenever a variant is found, a breadth of metrics are extracted from the `DCS.bam`, and a few from the `SSCS.bam` (for features such as strand bias). 

This step takes under 2 hours per sample and is already significantly parallelised, such that processing more samples simultaneously is not faster. The output of this step is a VCF, with all the metrics/features stored in the INFO field.

**Contents:**
* `extended_pileup_extractor_e2.sh`: This is the Snakemake pipeline that accepts the `DCS.bam` and `SSCS.bam` and produces the verbose VCF files. The main output is `{sample}.final_combined.vcf.gz`.
* `extended_variant_extracting_config.yaml`: This is an example of the config file to pass to the Snakemake pipeline. The parameters (`n_chunks`) throttle the parallelisation of data extraction. There are several entries for different types of blacklists which are legacy from earlier versions; the path can indicate the `master_blacklist.bed` for all (this is suboptimal but costs virtually 0 in terms of compute efficiency).

---

## 2) Annotation with GnomAD and PoN (`02.2_snv_annotation/`)

The aim of this section is to annotate the output of the previous section with GnomAD and PoN information. It is also in this section that a CSV file containing the necessary metadata regarding each sample is constructed. The metadata is necessary as features such as age, cell type, and sequencing depth are required to tune the parameters in section 02.3 and calculate SNV burden in section 02.4.

This step is computationally very light and takes a few minutes in total.

**Contents:**
* `annotate_pon_gnomad.sh`: A bash script run at the CLI that annotates the `{sample}.final_combined.vcf.gz` from step 02.1 and produces `../Sample_files/{sample}.final_combined.vcf.gz`. It requires a `pon_list.csv`, which is just a list of paths to each VCF forming the PoN (an example is present in the directory). 
    * *Note:* For the moment, this has to be run manually on samples 1 by 1 at the CLI, but there clearly is the possibility to include this as the last step of 02.1.
* `prepare_csv.py`: This script takes as input a CSV file with the metadata (`pre_processed_test_samples.csv`) containing the following columns: `sample_id, vcf_path, sample_type, age, mean_dcs_depth`. The script also requires the path to the `panel.bed` and the `blacklist.bed`. It calculates the effective territory (panel.bed - blacklist.bed) and produces an estimation of the overall number of callable bases per sample. The output is `processed_test_samples.csv`.
    * *Note:* This has to be run manually at the CLI after putting together the `pre_processed_test_samples.csv`. Due to the highly heterogeneous nature of the metadata and given its importance, it is preferable to keep this as a manual step.

---

## 3) Threshold Tuning (`02.3_snv_training/`)

This section is where the magic happens. A carefully crafted loss function estimates the error between the expected SNV burden (in mutations per bp) and the observed SNV burden on training data. The Optuna library is used to tune thresholds and reduce error. Good results can be obtained with a few hundred cycles of training; optimal performance requires ~3000 cycles, which takes a few hours. 

**Inputs & Outputs:**
* **Input:** `processed_test_samples.csv` (which itself contains the path to the annotated VCF files).
* **Output:** A `.json` file with the optimal thresholds which can then be loaded for filtering of test samples.

**Contents:**
* `snv_training_thresholds.ipynb`: The training script.
* `snv_best_thresholds.json`: Current best thresholds.

---

## 4) Application of Thresholds and Final VCF Production (`02.4_snv_test_and_plot/`)

This script accepts the exact same data as step 02.3: a `metadata.csv` file which itself contains the path to the VCF. It also requires the `snv_best_thresholds.json`. This script applies the thresholds to each VCF and outputs a filtered VCF containing only the rows that have passed all thresholds. 

> **Manual Tweak:** A key threshold that I often manually tweak after training is the maximum allowed VAF, which I bring down to 20% for a maximally conservative approach to guard against germline variants.

This script also outputs an `observed_vs_expected.csv` file with the expected and actual numbers of SNVs per sample as a function of depth, age, and cell type, to plot the age gradient.

**Contents:**
* `snv_test_and_plot.ipynb`: The script that performs the filtering. A second cell loads the `observed_vs_expected.csv` and plots burden vs. age.