# Indel Calling Validation & QC

## Overview

This notebook contains the validation and quality control analyses used to assess the performance of the indel calling pipeline. The aim is to demonstrate that the variant calls:

- Show a biologically plausible mutation burden with a clear age gradient
- Exhibit the expected indel size and type distributions (insertions vs. deletions)
- Are not dominated by technical artefacts (homopolymer context, strand bias, read position)
- Match known somatic indel mutational signatures (COSMIC ID83 framework)

> This is not part of the core pipeline, but rather a **post hoc validation step** used throughout development and for sanity checking final outputs.

---

## Inputs

The notebook operates on the **filtered VCFs** produced in `03.4_indel_test_and_plot/`, along with associated metadata.

| Input | Description |
|---|---|
| Filtered VCF files | Post-thresholding indel variant calls (`.indel_candidates_annotated.filtered.vcf.gz`) |
| `master.csv` | Metadata (sample ID, age, sample type, sequencing depth, effective territory) pointing to filtered VCFs |

> **Location:** `data/indels/filtered_vcf/`

> **Environment:** `02_snv_calling_environment.yaml`

---

## Key Analyses

### 1 · Mutation Burden Assessment

Estimation of the indel mutation rate after filtering, expressed as **indels per megabase**, across three sample types: **cord blood**, **T-cells**, and **PBMCs**. Burden is plotted as a function of age with a linear regression overlay (slope, intercept, and R²).

A separate bar chart compares the **observed cord blood indel rate** against the **expected biological reference rate** (2×10⁻⁹ indels per base), as a direct sanity check on pipeline calibration.

> **Note:** One cord blood sample (`Cordblood_H1`) is excluded from burden plots; this is due to an issue in library prep affecting depth and suspected to have introduced additional errors

---

### 2 · Indel Diagnostic Suite (QC)

A five-panel diagnostic figure assessing the technical quality of called indels:

| Panel | Description | What to look for |
|---|---|---|
| Insertion vs. Deletion Ratio | Global count of INS vs. DEL | Roughly balanced, or biologically justified skew |
| Size Distribution | Indel length (1–4 bp and 5+ bp), split by type | Exponential decay (shorter indels should dominate) |
| Homopolymer Context | Flanking homopolymer length per indel type | Excess long homopolymers may indicate polymerase slippage artefacts |
| Normalised Read Position | Distribution of indels along the read (0–1) | Spikes near 0 or 1 suggest soft-clipping or alignment artefacts |
| Strand Bias | SSCS strand balance score | Distribution centred at ~0.5; extremes indicate technical bias |

---

### 3 · Mutational Signature Extraction (ID83)

Extraction and attribution of somatic indel signatures using the **SigProfiler** framework (ID83 context).

Steps:
1. Filtered VCFs are decompressed and aggregated into a single input directory
2. `SigProfilerMatrixGenerator` constructs the ID83 mutation matrix (GRCh38)
3. `SigProfilerExtractor` performs NMF-based signature extraction (1–5 signatures; 10 replicates for exploration, 100 for final results)
4. COSMIC ID83 activity scores are loaded and visualised as a stacked bar chart showing the proportional contribution of each signature across the cohort

> **Note:** Indel signature extraction is computationally heavier than SNV extraction. Use `nmf_replicates=10` for exploratory runs and `nmf_replicates=100` for final results.

---

### 4 · Recurrence Analysis & Hotspot Annotation

Assessment of variant recurrence across the cohort to distinguish true somatic hotspots from technical artefacts.

Steps:
1. Each indel is assigned a unique identifier (`CHROM:POS_REF>ALT`)
2. Recurrence is calculated as the number of samples sharing the exact same variant
3. A histogram of recurrence counts is plotted to assess the cohort-wide distribution
4. The **top 20 recurrent variants** are annotated with gene names via the **Ensembl REST API** (GRCh38)
5. A blacklist check flags variants in known artefact-prone genes (`MUC16`, `TTN`, `HLA-A/B/C`, `MUC19`, `CDC27`)

> **Note:** Recurrent variants in blacklisted genes should be investigated and potentially excluded from downstream analyses.

---

## Outputs

| Output | Use |
|---|---|
| Indel burden vs. age plots | Burden assessment per sample type with age regression |
| Cord blood rate bar chart | Direct comparison of observed vs. expected indel rate |
| Five-panel QC diagnostic figure | Artefact identification and quality assessment |
| SigProfiler ID83 matrix | Input to signature extraction |
| Stacked signature activity bar chart | Mutational signature decomposition visualisation |
| Top 20 recurrent hotspot table | Candidate hotspot identification and artefact flagging |

---

## Notes

- The `Cordblood_H1` sample is excluded from burden plotting; this exclusion is intentional, there was an issue in library prep leading to low depth and suspected additional error.
- SigProfiler outputs are written to `sigprofiler_output/` and the ID83 matrix input to `sigprofiler_input_indels/`; neither directory is under version control
- The Ensembl REST API calls in the recurrence section require an internet connection and introduce a `0.1s` sleep between requests to respect rate limits
- Sample type standardisation (mapping `LEGACY` and `UKCTOCS` to `PBMCs`, `cord` to `cordblood`) is applied at load time and must be consistent with the SNV notebook to allow cross-variant-type comparisons

---

## Workflow Context

This notebook sits at the end of the indel calling workflow as a **final sanity check**:

```
Pipeline  →  Filtering (03.4)  →  Validation (this notebook)
```

The analyses mirror those in `notebooks/pre_denominator_correction/snvs/`, adapted for indel-specific biology (lower expected mutation rate, ID83 signature context, and indel-specific QC metrics).

---

## Future Improvements

- Formalise COSMIC ID83 signature decomposition and compare to published blood/ageing indel signatures (e.g. ID1, ID2 for APOBEC; ID8 for DSB repair)
- Extend recurrence analysis to flag variants shared with the Panel of Normals (PoN) that may have escaped filtering
- Cross-compare indel burden age gradients with the SNV burden to assess the relative contribution of each variant class per sample
