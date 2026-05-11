# 03 Indel Calling Pipeline

### Overview
This section is for indel calling. It perfectly mirrors the architecture and workflow of the SNV calling pipeline. The 4 sections are the same, contain the same (analogous) files, the usage is the same, and the environment is the same (`02_snv_calling_environment.yaml`).

There is no practical difference to the use, only technical differences in the code to accommodate indels vs SNVs. Because of the high similarity between the 2 processes, it is conceivable that they could be merged as one stream, though that would be quite delicate for limited benefit.

---

### Key Differences from SNV Calling
While the workflow execution is identical, there are critical technical differences in the underlying code to handle the distinct biology of indels:

* **Mutation Rate:** Indels have a much lower mutation rate than SNVs, and this is reflected in the loss function of `03.3_indel_training`.
* **Feature Extraction:** There are features that do not translate between SNVs and indels. For example:
    * **Base Quality:** Base quality doesn't mean anything for deletions (we use anchor base quality instead), and we use the average base quality for insertions.
    * **New Metrics:** Metrics such as indel length are added.
    * **NM Calculation:** There are changes to the NM calculation (the distance between the read and the reference).

---

### Workflow Sections
*(Note: These sections perfectly mirror the architecture of the `02_snv_calling` pipeline. Please refer to that README for detailed execution instructions.)*

1. **`03.1_indel_extraction/`**: Raw pileup and feature extraction.
2. **`03.2_indel_annotation/`**: Annotation with GnomAD and panel of normals (PoN).
3. **`03.3_indel_training/`**: Threshold tuning.
4. **`03.4_indel_test_and_plot/`**: Application of thresholds and final VCF production.