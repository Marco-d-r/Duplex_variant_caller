# 01 Consensus Calling Pipeline

### Overview
Pipeline to process fastq produced by the sequencing facility into usable bams. This is Venetia's latest pipeline, without 3' or 5' clipping of reads, as this can be performed in filtering.

The pipeline is designed to accept multiple files per sample, as we often sequence the same sample on multiple lanes (this is the first step: concat fastqs).

### Requirements & Resources
* **Configuration:** The config file shows the type of input that is necessary.
* **External Files:** The pipeline requires a **BED file** (if doing panel sequencing) and a **reference genome** for alignment with `bwa`.
* **Software:** This pipeline mostly uses the `fgbio` tools to process the fastqs to bams. 
    * [Link to fgbio help page](https://github.com/fulcrumgenomics/fgbio)

### Key Outputs
1. **Duplex Consensus BAM (`DCS.bam`):** This is the bam file with consensuses from both top and bottom strands, which should be used for variant calling.
2. **Single Strand Consensus BAM (`SSCS.bam`):** This file contains consensuses before duplexing (so consensus on top strand and bottom strand are separate). It is of limited utility as the error rate is higher than `DCS.bam`.
3. **Raw Mapped BAM (`raw.mapped.bam`):** This is the file without consensus; there has been very little filtering, mostly conversion from fastq to bam and alignment of the reads. 
    * **Use Case:** Can be useful to determine depth accurately (as many reads that fail to form consensus are lost and not found in `DCS.bam`, leading to biased depth assessment). 
    * > **WARNING:** This can be an extremely large file depending on input.
4. **Key Metrics:** This useful file provides a summary of the `SSCS.bam` and `DCS.bam`, with information such as average depth and off-target fraction.

---

### HPC Execution
This pipeline is meant to be run on the HPC. It defaults to the **SL3** account (hardcoded in the profile) but uses the paying **SL2** account for the alignment of the raw bam, as this is often longer than 12 hours. 

**Setup Requirements:**
* A slurm profile (an example profile, `fgbio_profile`, is included in this directory).
* An environment supporting Snakemake v7 (see `01.1_fgbio_pipeline_environment.yaml` for dependencies).

**Typical running CLI command:**
```bash
snakemake \
  --snakefile /path/to/01.1_fgbio_pipeline.sh \
  --configfile /path/to/example_fgbio_pipeline_config.yaml \
  --profile /path/to/fgbio_profile \
  --jobs 50 \
  --rerun-incomplete \
  --use-conda