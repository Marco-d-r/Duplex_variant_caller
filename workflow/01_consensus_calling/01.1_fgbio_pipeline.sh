# Title: Snakefile for running Fgbio pipeline on HPC
# Author: Venetia D'Arcy
# Date: 09-02-2026

# This Snakefile processes paired-end FASTQ.GZ files to generate SSCS and DCS consensus reads using fgbio.
# It is designed to run on the HPC with each rule submitted as an individual job.

# It corresponds to the fgbio pipeline V4 which is running on the server.
# Version 4 removes clipping entirely.
# This version removes previous issue where the max no call fraction was set to =1. 
# It does not include variant calling or ASXS filtering (these will be run separately)

# Required software or dependencies:
#   snakemake
#   bwa
#   samtools v1.14 or higher
#   fgbio v2.0.0 of higher ==> itself requiring java 1.8 or higher

# Run with:
# snakemake \
#   --snakefile Snakefile \
#   --configfile config_file.yaml \
#   --profile fgbio_snakemake \
#   --jobs 50 \
#   --use-conda

###########################################
# Load config file
###########################################
VERSION = config["version"]
REFERENCE = config["reference"]
fgbio_jar = config["fgbio_jar"]
TEMP_OUTPUT_DIR = config["temp_output_dir"]
FINAL_OUTPUT_DIR = config["final_output_dir"]
CONCAT_DIR = config["concat_dir"]

SAMPLES = [s["id"] for s in config["samples"]]

SAMPLE_TO_R1  = {s["id"]: s["read1"] for s in config["samples"]}
SAMPLE_TO_R2  = {s["id"]: s["read2"] for s in config["samples"]}
SAMPLE_TO_PANEL = {s["id"]: s["library"] for s in config["samples"]}

###########################################
# Target output files
###########################################
rule all:
    input:
        expand(FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.raw.bam", sample=SAMPLES, version=VERSION),
        expand(FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.dcs.bam", sample=SAMPLES, version=VERSION),
        expand(FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.sscs.bam", sample=SAMPLES, version=VERSION),
        expand(FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.key_metrics.csv", sample=SAMPLES, version=VERSION)

###########################################
# Concatenate FASTQs if needed
###########################################
rule concat_fastqs:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample],
        r2=lambda wc: SAMPLE_TO_R2[wc.sample]
    output:
        r1=temp(CONCAT_DIR + "/{sample}.R1.concat.fastq.gz"),
        r2=temp(CONCAT_DIR + "/{sample}.R2.concat.fastq.gz")
    threads: 2
    resources:
        mem_mb=2000,
        time="00:30:00",
        partition="icelake"
    log:
        "logs/concat.{sample}.log"
    shadow: "minimal"
    shell:
        """
        set -o pipefail
        mkdir -p /home/mvd25/rds/rds-general_use-Q5jpz9ZAQ7k/Marco/CUH832_snakemake_7/temporary_output_folder
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """


###########################################
# FGBIO pipeline
###########################################

# Step 1: FASTQ -> uBam
rule fastq_to_ubam:
    """
    Generates a uBam from concatenated FASTQs.
    """
    input:
        r1=CONCAT_DIR + "/{sample}.R1.concat.fastq.gz",
        r2=CONCAT_DIR + "/{sample}.R2.concat.fastq.gz"
    output:
        bam=temp("/home/mvd25/rds/rds-general_use-Q5jpz9ZAQ7k/Marco/CUH832_snakemake_7/temporary_output_folder/{sample}.{version}.unmapped.raw.bam")
    params:
        rs1=config["r1_read_structure"],
        rs2=config["r2_read_structure"],
        library=lambda wc: SAMPLE_TO_PANEL[wc.sample],
        sequencing_center=config.get("sequencing_center", "CRUK"),
        platform=config.get("platform", "illumina"),
        platform_model=config.get("platform_model", "NovaSeqX_25B")
    threads: 16
    resources:
        mem_mb=225000,
        time="12:00:00",
        partition="icelake-himem"
    log:
        "logs/fastq_to_ubam.{sample}.{version}.log"
    benchmark:
        "benchmarks/fastq_to_ubam.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        set -o pipefail
        ulimit -n 8000
        /usr/bin/time -v java -Xmx20G -Djava.io.tmpdir=$PWD \
            -jar {fgbio_jar} --compression 1 --async-io FastqToBam \
            --input {input.r1} {input.r2} \
            --read-structures {params.rs1} {params.rs2} \
            --sample {wildcards.sample} \
            --library {params.library} \
            --sequencing-center {params.sequencing_center} \
            --platform {params.platform} \
            --platform-model {params.platform_model} \
            --output {output.bam} \
        &> {log}
        """


### Step 2: uBam -> Mapped BAM
# Align the data with bwa, then recover headers and tags from the unmapped BAM.
#This step uses a combination of `samtools` (to convert the uBam back to FASTQ), `bwa mem` to do the actual alignment and `fgbio ZipperBams` to carry metadata and UMIs forward from the unmapped BAM.  

rule align_raw_bam:
    """
    Takes an unmapped BAM, aligns it using BWA, and adds metadata back with ZipperBams
    """
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.raw.bam",
        ref=config["reference"]
    output:
        bam = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.raw.bam"
    threads: 20
    resources:
        mem_mb=48000,
        time="24:00:00",
        partition="icelake",
        slurm_account="JBLUNDELL-MVD25-SL2-CPU"
    log:
        "logs/align_bam.{sample}.{version}.log"
    benchmark:
        "benchmarks/align_raw_bam.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        set -o pipefail && \
        (
        ulimit -n 8000
        /usr/bin/time -v samtools fastq -@ {threads} {input.bam} \
        | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - \
        | java -Xmx44g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 --async-io ZipperBams \
            --unmapped {input.bam} \
            --ref {input.ref} \
            --output {output.bam} \
        ) &> {log}
        """

# Some important `bwa mem` options used here are:
# - `-K 150000000` to tell `bwa` to process reads/bases in fixed size chunks regardless of the number of threads; this makes `bwa` deterministic and reproducible with different thread counts
# - `-Y` to use soft-clipping (instead of hard clipping) in supplemental alignments
# Under no circumstances should `-M` be used to `mark shorter split hits as secondary` as this will cause problems with downstream tools.


### Step 3: Mapped BAM -> Grouped BAM
#This step identifies reads or read pairs that originate from the same source molecule based on genomic positions and UMIs.

# Group the reads by position and UMI

rule group_reads_by_umi:
    """Group the raw reads by UMIs and position ready for consensus calling."""
    input:
        bam = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.raw.bam",
    output:
        bam = temp(TEMP_OUTPUT_DIR + "/{sample}.{version}.grouped.bam"),
        family_sizes = FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.tag-family-sizes.txt",
        grouping_metrics = FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.grouping-metrics.txt"
    params:
        allowed_edits=config["UMI_edits"],
        min_map_q=config["min_map_q"]
    threads: 5
    resources:
        mem_mb=34000,
        time="12:00:00",
        partition="icelake"
    log:
        "logs/group_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/group_reads_by_umi.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
       """
       (
        ulimit -n 8000
        /usr/bin/time -v java -Xmx32g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 --async-io GroupReadsByUmi \
            --input {input.bam} \
            --strategy paired \
            --edits {params.allowed_edits} \
            --output {output.bam} \
            --family-size-histogram {output.family_sizes} \
            --grouping-metrics {output.grouping_metrics} \
            --min-map-q {params.min_map_q} \
        ) &> {log}
        """

## Step 4A: SSCS calling from the grouped bam file

# Step 4A.1: GroupedBam -> SSCS uBam --> Filtered SSCS uBam
#This step generates unmapped SSCS reads from the grouped reads and immediately filters them.  

rule call_sscs:
    """Generate unmapped single strand consensus from grouped bam and immediately filter this."""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.grouped.bam",
        ref=config["reference"]
    output:
        bam = temp(TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.sscs.bam"),
        sscs_metrics = FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.sscs.metrics.txt"
    threads: 1
    resources:
        mem_mb=8000,
        time="12:00:00",
        partition="sapphire"
    params:
        min_sscs_reads=config["min_sscs_reads"],
        min_input_base_qual=config["min_input_base_quality"]
    log:
        "logs/call_sscs_consensus_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/call_sscs.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        ulimit -n 8000
        set -o pipefail && \
        (
         /usr/bin/time -v java -Xmx30g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 CallMolecularConsensusReads \
            --input {input.bam} \
            --output {output.bam} \
            --min-reads {params.min_sscs_reads} \
            --min-input-base-quality {params.min_input_base_qual} \
            --stats {output.sscs_metrics} \
        ) &> {log}
        """

rule filter_sscs:
    """Generate unmapped single strand consensus from grouped bam and immediately filter this."""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.sscs.bam",
        ref=config["reference"]
    output:
        bam = temp(TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.sscs.filtered.bam"),
    threads: 5
    resources:
        mem_mb=32000,
        time="12:00:00",
        partition="icelake"
    params:
        min_sscs_reads=config["min_sscs_reads"],
        min_filter_base_qual=config["min_filter_base_quality"],
        max_base_error_rate=config["max_base_error_rate"],
        max_no_calls=config["max_no_calls"],
    log:
        "logs/filter_sscs_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/filter_sscs.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        ulimit -n 8000
        set -o pipefail && \
        (
         /usr/bin/time -v java -Xmx30g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 FilterConsensusReads \
            --input {input.bam} \
            --output {output.bam} \
            --ref {input.ref} \
            --min-reads {params.min_sscs_reads} \
            --min-base-quality {params.min_filter_base_qual} \
            --max-base-error-rate {params.max_base_error_rate} \
            --max-no-calls {params.max_no_calls} \
        ) &> {log}
        """

#Step 4A.2: Filtered SSCS uBam -> Mapped SSCS BAM

rule map_sscs:
    """Map, sort, and index the single strand consensus"""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.sscs.filtered.bam",
        ref=config["reference"]
    output:
        bam = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.sscs.bam",
        bai = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.sscs.bam.bai"
    threads: 16
    resources:
        mem_mb=180000,
        time="12:00:00",
        partition="icelake-himem"
    log:
        "logs/map_sscs_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/map_sscs.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        ulimit -n 8000
        set -o pipefail && \
        (/usr/bin/time -v samtools fastq -@ {threads} {input.bam} \
            | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - \
            | java -Xmx90g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 0 --async-io ZipperBams \
               --unmapped {input.bam} \
               --ref {input.ref} \
               --tags-to-reverse Consensus \
               --tags-to-revcomp Consensus \
            | samtools sort -m 3G --threads {threads} -T $PWD/sort_tmp -o {output.bam}##idx##{output.bam}.bai --write-index \
        ) &> {log} 
        """


#Step 4B: Generating filtered DCS reads from the reads grouped by UMI (same input as phase 4A)
#This step calls duplex consensus reads:

rule call_dcs:
    """Generate unmapped duplex consensus from grouped bam"""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.grouped.bam"
    output:
        bam = temp(TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.dcs.bam"),
        dcs_metrics = FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.dcs.metrics.txt"
    threads: 1
    resources:
        mem_mb=16000,
        time="12:00:00",
        partition="sapphire"
    params:
       min_reads=config["min_dcs_reads"],
    log:
        "logs/call_duplex_consensus.{sample}.{version}.log"
    benchmark:
        "benchmarks/call_dcs.{sample}.{version}.txt"
    shell:
        """
        (
        ulimit -n 8000
        set -o pipefail && \
        /usr/bin/time -v java -Xmx14g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 CallDuplexConsensusReads \
          --input {input.bam} \
          --output {output.bam} \
          --min-reads {params.min_reads} \
          --stats {output.dcs_metrics} \
        ) &> {log}
        """
#Can add --trim option to the above to quality trim input reads in addition to masking low Q bases - need to trial this option to see how it affects the output

#Step 4B.2: Unfiltered DCS uBam -> Filtered DCS uBam
rule filter_dcs:
    """Filter the duplex consensus"""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.dcs.bam",
        ref=config["reference"]
    output:
        bam = temp(TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.dcs.filtered.bam")
    threads: 3
    resources:
        mem_mb=16000,
        time="12:00:00",
        partition="icelake"
    params:
        min_dcs_reads=config["min_dcs_reads"],
        min_filter_base_qual=config["min_filter_base_quality"],
        max_base_error_rate=config["max_base_error_rate"],
        max_no_calls=config["max_no_calls"],
    log:
        "logs/filter_duplex_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/filter_dcs.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        (
        ulimit -n 8000
        /usr/bin/time -v java -Xmx14g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 1 FilterConsensusReads \
        --input {input.bam} \
        --output {output.bam} \
        --ref {input.ref} \
        --min-reads {params.min_dcs_reads} \
        --min-base-quality {params.min_filter_base_qual} \
        --max-base-error-rate {params.max_base_error_rate} \
        --max-no-calls {params.max_no_calls} \
        ) &> {log}
        """

#Step 4B.3: Filtered DCS uBam -> Mapped DCS BAM
rule map_dcs:
    """Map, sort, and index the duplex reads"""
    input:
        bam = TEMP_OUTPUT_DIR + "/{sample}.{version}.unmapped.dcs.filtered.bam",
        ref=config["reference"]
    output:
        bam = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.dcs.bam",
        bai = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.dcs.bam.bai"
    threads: 12
    resources:
        mem_mb=180000,
        time="12:00:00",
        partition="icelake-himem"
    log:
        "logs/map_duplex_reads.{sample}.{version}.log"
    benchmark:
        "benchmarks/map_dcs.{sample}.{version}.txt"
    shadow: "minimal"
    shell:
        """
        ulimit -n 8000
        set -o pipefail && \
        (/usr/bin/time -v samtools fastq -@ {threads} {input.bam} \
            | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - \
            | java -Xmx90g -Djava.io.tmpdir=$PWD -jar {fgbio_jar} --compression 0 --async-io ZipperBams \
               --unmapped {input.bam} \
               --ref {input.ref} \
               --tags-to-reverse Consensus \
               --tags-to-revcomp Consensus \
            | samtools sort -m 2G --threads {threads} -T $PWD/sort_tmp -o {output.bam}##idx##{output.bam}.bai --write-index \
        ) &> {log} 
        """

#Step 5: Producing metrics
rule key_metrics_summary:
    """
    Extract key metrics from SSCS and DCS BAMs restricted to BED targets.
    Handles empty BAMs and avoids Bash strict mode exit issues.
    """
    input:
        sscs = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.sscs.bam",
        dcs  = FINAL_OUTPUT_DIR + "/{sample}.{version}.mapped.dcs.bam",
        bed = lambda wc: config[f"{SAMPLE_TO_PANEL[wc.sample]}_bed"] #this loads the relevant bed file for the sample based on the library/panel specified in the config file
    output:
        csv = FINAL_OUTPUT_DIR + "/metrics/{sample}.{version}.key_metrics.csv"
    log:
        "logs/key_metrics_summary.{sample}.{version}.log"
    benchmark:  
        "benchmarks/key_metrics_summary.{sample}.{version}.txt"    
    threads: 
        2
    resources:
        mem_mb=1220,
        time="05:00:00",
        partition="icelake"
    shell:
        r"""
        set -euo pipefail

        bed="{input.bed}"
        sscs="{input.sscs}"
        dcs="{input.dcs}"
        csv="{output.csv}"

        metrics_sorted() {{
          BAM=$1
          LABEL=$2
          BED=$3

          mean_depth=$(samtools depth -b "$BED" "$BAM" | awk '{{sum+=$3;n++}}END{{if(n>0)print sum/n; else print 0}}')

          read cov1 cov10 cov50 <<<$(samtools depth -b "$BED" "$BAM" | \
            awk '{{c1+=($3>=1);c10+=($3>=10);c50+=($3>=50);n++}}END{{if(n>0)printf "%f %f %f\n",c1/n*100,c10/n*100,c50/n*100; else print "0 0 0"}}')

          mean_len=$(samtools view -L "$BED" "$BAM" | awk '{{sum+=length($10);n++}}END{{if(n>0)print sum/n; else print 0}}')

          total=$(samtools view -c -L "$BED" "$BAM")
          proper=$(samtools view -c -f 2 -L "$BED" "$BAM")
          if [ "$total" -gt 0 ]; then
            proper_prop=$(awk -v p=$proper -v t=$total 'BEGIN{{print p/t*100}}')
          else
            proper_prop=0
          fi

          mean_ins=$(samtools view -f 2 -L "$BED" "$BAM" | awk '{{sum+=($9>=0?$9:-$9);n++}}END{{if(n>0)print sum/n; else print 0}}')

          on=$(samtools view -c -L "$BED" "$BAM")
          all=$(samtools view -c "$BAM")
          if [ "$all" -gt 0 ]; then
            off_frac=$(awk -v on=$on -v all=$all 'BEGIN{{print (1-on/all)*100}}')
          else
            off_frac=0
          fi

          echo "$LABEL,$mean_depth,$cov1,$cov10,$cov50,$mean_len,$proper_prop,$mean_ins,$off_frac"
        }}

        # header
        echo "Type,MeanDepth,Cov>=1x(%),Cov>=10x(%),Cov>=50x(%),MeanReadLength,ProperPairProp(%),MeanInsertSize,OffTargetFraction(%)" > "$csv"

        # metrics
        metrics_sorted "$sscs" SSCS "$bed" >> "$csv"
        metrics_sorted "$dcs"  DCS  "$bed" >> "$csv"

        echo "Metrics written to $csv"
        """
