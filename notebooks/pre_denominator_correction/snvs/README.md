Here's a cleaned-up version with better formatting — tighter hierarchy, consistent style, and cleaner callout blocks:

---

# SNV Calling Validation & QC

## Overview

This notebook contains the validation and quality control analyses used to assess the performance of the SNV calling pipeline. The aim is to demonstrate that the variant calls:

- Have a low and well-controlled error rate
- Show biologically plausible mutational patterns
- Are not dominated by technical artefacts

> This is not part of the core pipeline, but rather a **post hoc validation step** used throughout development and for sanity checking final outputs.

---

## Inputs

The notebook operates on the **filtered VCFs** produced in `02.4_snv_test_and_plot/`, along with associated metadata.

| Input | Description |
|---|---|
| Filtered VCF files | Post-thresholding variant calls |
| `master.csv` | Metadata (age, sample type, depth, etc.) pointing to filtered VCFs |

> **Environment:** runs on a kernel supported by `r_dndscv_env`

---

## Key Analyses

### 1 · Mutation Burden Assessment

Estimation of the mutation rate after filtering across three sample types: **cord blood**, **T-cells**, and **PBMCs** — each with their own expected mutation rate.

Plots compare mutation burden as a function of cell type and age, showing an age gradient and biologically plausible mutation rates.

> **Note:** Filtering may be conservative — rates are slightly but consistently lower than expected. This is acceptable: the pipeline prioritises **specificity over sensitivity**, in keeping with the single-molecule calling objective.

---

### 2 · Trinucleotide Context / SBS Signatures

Analysis of the trinucleotide context of called SNVs. Mutations are grouped into one of **96 trinucleotide contexts** (e.g. `T[C>A]G`) and profiles are compared quantitatively to known profiles from single-cell derived colony data.

> **Result:** Cosine similarity = **0.88** vs. gold standard.

The gold standard is derived from PBMC WGS data. Our cohort includes T-cells and PBMCs assayed by WES, which likely accounts for some of the divergence — an important avenue for follow-up investigation.

---

### 3 · dN/dS Analysis

Assessment of selection pressure using dN/dS ratios via **dNdScv**.

> **Result:** Overall dN/dS ≈ **1.0**

This is highly reassuring:
- Germline variants typically yield dN/dS ratios of 0.1–0.3, so a ratio of ~1 indicates **no germline contamination**
- A value of ~1 is consistent with **neutral selection**, as expected — at our sequencing depth (200×+), the bulk of captured mutations should reflect drift and hitchhiking rather than positive selection

---

## Outputs

| Output | Use |
|---|---|
| Mutation burden plots | Burden vs. age, depth, and other covariates |
| Trinucleotide mutation spectra | Mutational signature visualisation |
| dN/dS estimates | Selection pressure summary |
| Diagnostic plots | Artefact identification |

These figures are used for pipeline validation during development, and for inclusion in reports and presentations (e.g. viva).

---

## Notes

- This notebook is intentionally exploratory in parts, but key plots have been stabilised over time
- Some thresholds and plotting choices may be manually adjusted depending on the dataset

---

## Workflow Context

This step acts as a **final sanity check** after variant calling:

```
Pipeline  →  Filtering  →  Validation (this notebook)
```

---

## Future Improvements

- Formalise SBS signature fitting (e.g. COSMIC decomposition)
- Assess mutation rate and mutational signatures as **WES vs. WES** rather than WES vs. WGS using single-cell derived colony data