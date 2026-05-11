# Duplex Variant Caller

A complete, end-to-end pipeline for calling somatic variants (SNVs and indels) from **duplex sequencing** data. The pipeline takes raw FASTQ files through consensus calling, pileup and feature extraction, annotation, threshold tuning, and final filtered VCF production. A validation/QC notebook is included to assess pipeline performance and mutational signatures.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Workflow](#workflow)
- [Environments](#environments)
- [Data](#data)
- [Notebooks](#notebooks)
- [Quickstart](#quickstart)

---

## Overview

Duplex sequencing achieves ultra-low error rates by forming consensus reads from both strands of a DNA duplex. This pipeline leverages that property to call **rare somatic variants** with high specificity, prioritising a low false-positive rate over sensitivity — consistent with single-molecule variant calling objectives.

The pipeline covers:

1. **Consensus calling** — FASTQ → consensus BAM (DCS + SSCS) via `fgbio`, run on HPC
2. **SNV calling** — pileup, annotation, threshold tuning (Optuna), and filtered VCF production, run locally
3. **Indel calling** — mirrors the SNV pipeline with indel-specific feature extraction and a lower expected mutation rate
4. **Validation** — mutation burden, trinucleotide signatures (cosine similarity 0.88 vs gold standard), and dN/dS analysis (~1.0, confirming no germline contamination)

---

## Repository Structure

```
Duplex_variant_caller/
├── data/
│   ├── snvs/
│   │   ├── filtered_3_13_vcf/      # Post-threshold filtered VCFs + master_list.csv
│   │   └── Unfiltered_vcf/         # Raw pileup VCFs from step 02.1
│   ├── indels/                     # Indel VCFs (mirroring snvs/)
│   └── machado_trinucleotide/      # External reference SNV data (Machado 2022)
│
├── docs/
│   └── rationale.txt               # Design decisions and biological rationale
│
├── environments/
│   ├── 01.1_fgbio_pipeline_environment.yaml
│   ├── 02_snv_calling_environment.yaml
│   ├── r_dndscv_env.yaml
│   └── vep_env.yaml
│
├── notebooks/
│   ├── pre_denominator_correction/
│   │   ├── snvs/                   # SNV validation & QC (burden, signatures, dN/dS)
│   │   └── indels/                 # Indel validation & QC
│   └── post_denominator_correction/
│
└── workflow/
    ├── 01_consensus_calling/       # FASTQ → DCS/SSCS BAMs (HPC, Snakemake + fgbio)
    ├── 02_snv_calling/             # SNV calling: extraction → annotation → training → filtering
    └── 03_indel_calling/           # Indel calling: same 4-step architecture as SNV
```

---

## Workflow

The pipeline is split into three sequential workflow stages. Steps 02 and 03 are independent of each other and can be run in parallel after step 01.

```
FASTQ files
    │
    ▼
[01] Consensus Calling  (HPC — Snakemake + fgbio)
    │   → DCS.bam, SSCS.bam, raw.mapped.bam, key_metrics
    │
    ├──────────────────────┐
    ▼                      ▼
[02] SNV Calling        [03] Indel Calling   (local)
  02.1 Extraction         03.1 Extraction
  02.2 Annotation         03.2 Annotation
  02.3 Training           03.3 Training
  02.4 Filtering          03.4 Filtering
    │                      │
    └──────────┬───────────┘
               ▼
        Filtered VCFs
               │
               ▼
    [Notebooks] Validation & QC
```

### Step 01 — Consensus Calling

Processes FASTQs (including multi-lane samples) into consensus BAMs using `fgbio`. Run on HPC with a Slurm profile; alignment of the raw BAM uses the SL2 (paying) account due to runtime. See [`workflow/01_consensus_calling/README.md`](workflow/01_consensus_calling/README.md).

### Step 02 — SNV Calling

Four sequential sub-steps, all run **locally** after downloading BAMs from HPC:

| Sub-step | Script | Description |
|---|---|---|
| `02.1_snv_extraction` | `extended_pileup_extractor_e2.sh` | Raw pileup; extracts features into a verbose VCF |
| `02.2_snv_annotation` | `annotate_pon_gnomad.sh` + `prepare_csv.py` | Annotates with GnomAD and Panel of Normals; builds metadata CSV |
| `02.3_snv_training` | `snv_training_thresholds.ipynb` | Optimises filtering thresholds via Optuna (~3000 cycles for best results) |
| `02.4_snv_test_and_plot` | `snv_test_and_plot.ipynb` | Applies thresholds; outputs filtered VCFs and burden vs. age plots |

See [`workflow/02_snv_calling/README.md`](workflow/02_snv_calling/README.md).

### Step 03 — Indel Calling

Exactly mirrors the 4-step SNV architecture. Key technical differences include indel-specific feature extraction (anchor base quality, indel length, adjusted NM calculation) and a lower expected mutation rate in the loss function. See [`workflow/03_indel_calling/README.md`](workflow/03_indel_calling/README.md).

---

## Environments

All Conda environments are in `environments/`. Each workflow step specifies which environment it requires.

| File | Used by |
|---|---|
| `01.1_fgbio_pipeline_environment.yaml` | Step 01 — consensus calling |
| `02_snv_calling_environment.yaml` | Steps 02 and 03 — SNV and indel calling |
| `r_dndscv_env.yaml` | Validation notebooks (dN/dS analysis) |
| `vep_env.yaml` | VEP annotation (used in notebooks) |

```bash
conda env create -f environments/<env_file>.yaml
```

---

## Data

`data/` holds the VCF files that are inputs and outputs of the pipeline. It is not under active version control for large binary files (`.vcf.gz`, `.tbi`).

- `data/snvs/Unfiltered_vcf/` — raw pileup VCFs from step 02.1, one per sample
- `data/snvs/filtered_3_13_vcf/` — post-threshold filtered VCFs + `master_list.csv` (sample metadata)
- `data/machado_trinucleotide/` — external SNV reference data from Machado *et al.* 2022, used for mutational signature comparison

Sample naming follows the convention `{cohort}_{ID}` (e.g. `Cordblood_H1`, `UKCTOCS_UDI3`, `CUH832_TILS_CD8`).

---

## Notebooks

Validation and exploratory notebooks live in `notebooks/`. They are intended as **post hoc QC**, not part of the core pipeline.

### `notebooks/pre_denominator_correction/snvs/`

Three key analyses used to validate pipeline performance:

**1. Mutation Burden** — estimates SNV rate across cord blood, T-cells, and PBMCs. Rates are slightly conservative (specificity is prioritised over sensitivity), but show the expected age gradient.

**2. Trinucleotide Signatures** — profiles called SNVs across 96 trinucleotide contexts and compares to single-cell colony WGS data.
> Cosine similarity = **0.88** vs. gold standard (PBMC WGS). Residual divergence attributed to WES vs. WGS differences.

**3. dN/dS Analysis** — uses `dNdScv` to test for selection pressure.
> dN/dS ≈ **1.0**, consistent with neutral somatic drift and confirming the absence of germline contamination (germline variants typically yield dN/dS of 0.1–0.3).

Environment: `r_dndscv_env.yaml`

---

## Quickstart

```bash
# 1. Clone the repo
git clone https://github.com/Marco-d-r/Duplex_variant_caller.git
cd Duplex_variant_caller

# 2. Set up the relevant environment (example: SNV calling)
conda env create -f environments/02_snv_calling_environment.yaml
conda activate 02_snv_calling

# 3. Run consensus calling on HPC (see workflow/01_consensus_calling/README.md)
snakemake \
  --snakefile workflow/01_consensus_calling/01.1_fgbio_pipeline.sh \
  --configfile workflow/01_consensus_calling/example_fgbio_pipeline_config.yaml \
  --profile workflow/01_consensus_calling/fgbio_profile \
  --jobs 50 --rerun-incomplete --use-conda

# 4. Download DCS.bam and SSCS.bam from HPC, then run SNV calling locally
#    Follow the sub-steps in workflow/02_snv_calling/README.md
```

---

## Reference

Machado *et al.* 2022 — used as external reference for trinucleotide mutational signature comparison (`data/machado_trinucleotide/`).
