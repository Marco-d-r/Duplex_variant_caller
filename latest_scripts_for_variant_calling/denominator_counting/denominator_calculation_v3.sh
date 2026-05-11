configfile: "config.yaml"

SAMPLES = config["samples"].keys()
CHROMS = config["chroms"]

rule all:
    input:
        "results/master_denominator_report.csv"

rule preprocess_sample_bed:
    input:
        wes = config["wes_bed"],
        blacklist = config["blacklist_bed"],
        vcf = lambda wildcards: config["samples"][wildcards.sample]["vcf"]
    output:
        bed = "results/{sample}/{sample}_clean_targets.bed",
        stats = "results/{sample}/{sample}_mask_stats.json"
    script:
        "scripts/preprocess_bed.py"

rule count_callable_by_chrom:
    """The heavy lifting: calculates counts for one chromosome."""
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["bam"],
        bed = "results/{sample}/{sample}_clean_targets.bed"
    output:
        json = "results/{sample}/chrom_counts/{chrom}.json"
    threads: 1 # Snakemake will run 14 of these in parallel
    script:
        "scripts/count_callable.py"

rule aggregate_sample:
    input:
        # These must match the names used in your python script
        chrom_jsons = expand("results/{{sample}}/chrom_counts/{chrom}.json", chrom=config["chroms"]),
        mask_stats = "results/{sample}/{sample}_mask_stats.json"
    output:
        csv = "results/{sample}/total_counts.csv"
    script:
        "scripts/aggregate_results.py"

rule aggregate_all_samples:
    input:
        sample_csvs = expand("results/{sample}/total_counts.csv", sample=config["samples"].keys())
    output:
        master_csv = "results/master_denominator_report.csv"
    run:
        import pandas as pd
        
        # Read all individual CSVs into a list of DataFrames
        df_list = [pd.read_csv(f) for f in input.sample_csvs]
        
        # Concatenate them vertically (stacking rows)
        master_df = pd.concat(df_list, ignore_index=True)
        # Inside rule aggregate_all_samples run block:
        master_df['molecular_depth'] = master_df['final_callable_bases'] / master_df['final_target_size']
        master_df.to_csv(output.master_csv, index=False)
        
        # Save to the final destination
        master_df.to_csv(output.master_csv, index=False)