# Title: Snakefile for running SNV and indel calling pipeline on HPC
# Author: Venetia D'Arcy
# Date: 09-03-2026

# This Snakefile processes DCS and SSCS bam files to generate an annotated VCF file containing indel information.
# It is designed to run on the HPC with each rule submitted as an individual job.
# Version 2 has longer time limits set for the chunked jobs (so it can be submitted with a config file that has 6 chunks instead of 60, to avoid scheduling issues)

# It corresponds to Marco's scripts:
#       Extended_pileup_extractor_V2.sh
#       Extended_indel_extractor_V1.sh

# Required software or dependencies:
#   snakemake

# Run with:
# snakemake \
#   --snakefile SNV_and_indel_calling_pipeline.V1.sh \
#   --configfile config_file.yaml \
#   --profile fgbio_snakemake \
#   --jobs 50 \
#   --use-conda
#   --rerun-incomplete \
#   --restart-times 0


###########################################
# Load config file
###########################################
SAMPLES = list(config["samples"].keys())
REF = config["PATHS"]["ref_fasta"]
CONTIGS_HDR = config["PATHS"]["contigs_hdr"]
BED = config["PATHS"]["bed_targets"]
BED_CHUNKS_DIR = config["PATHS"]["bed_chunks_dir"]
PROCESSED_DATA_DIR = config["PATHS"]["processed_data_dir"]
N_CHUNKS = config["PARAMS"]["n_chunks_dcs"] # For DCS splitting
SSCS_CHUNKS = config["PARAMS"]["n_chunks_sscs"] # For SSCS splitting
DCS_CHUNK_IDS = [f"chunk_{i}" for i in range(1, N_CHUNKS + 1)]
SSCS_CHUNK_IDS = [f"chunk_{i}" for i in range(1, SSCS_CHUNKS + 1)]
GNOMAD_VCF = config["PATHS"]["gnomad_vcf"]
SNV_PON_CSV = config["PATHS"]["snv_pon_list"]
INDEL_PON_CSV = config["PATHS"]["indel_pon_list"]


localrules: 
    all, 
    split_bed_regions, 
    gather_mutect2_raw, 
    extract_indels_for_snv_processing,
    gather_dcs_snv_global_data

###########################################
# Target output files
###########################################
rule all:
    input:
        expand(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined.vcf.gz", sample=SAMPLES),
        expand(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_with_header.vcf.gz", sample=SAMPLES),
        expand(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined_annotated.vcf.gz", sample=SAMPLES),
        expand(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_annotated.vcf.gz", sample=SAMPLES)

###########################################
# PRE_PROCESSING FOR BOTH SNV AND INDEL calling
###########################################

# STEP 1.1: SPLIT BED INTO CHUNKS
# the idea here is to use the exome panel bed file to create small regions to parallelize the variant calling 
# the key metrics rule can run in parallel with this step 

rule split_bed_regions:
    input:
        BED
    output:
        chunks = temp(expand(f"{BED_CHUNKS_DIR}/chunk_{{i}}.bed", i=range(1, N_CHUNKS + 1)))
    params:
        n = N_CHUNKS,
        outdir = BED_CHUNKS_DIR
    threads: 1
    resources:
        mem_mb=2000,
        time="00:30:00"
    log:
        "logs/split_bed_regions.log"
    shell:
        """
        (
        set -euo pipefail
        mkdir -p {BED_CHUNKS_DIR}

        # 1. Count valid data lines (excluding headers)
        total_lines=$(grep -c -v "^#" {input} || true)
        
        # 2. Calculate lines per chunk (rounding up)
        lines_per_chunk=$(( (total_lines + {params.n} - 1) / {params.n} ))
        
        # 3. Clean existing chunks
        rm -f {BED_CHUNKS_DIR}/chunk_*.bed
        rm -f {BED_CHUNKS_DIR}/split_tmp_*
        
        # 4. Split the file 
        split -a 4 -l $lines_per_chunk {input} {BED_CHUNKS_DIR}/split_tmp_
        
        # 5. Rename and add .bed extension
        i=1
        for f in {BED_CHUNKS_DIR}/split_tmp_*; do
            mv "$f" "{BED_CHUNKS_DIR}/chunk_$i.bed"
            i=$((i + 1))
        done
        
        echo "Split {input} into $i chunks." 
        ) &> {log}
        """

# STEP 1.2: CALCULATE DEPTH PERCENTILES

rule calculate_depth_lookup:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        bed = config["PATHS"]["bed_targets"] 
    output:
        depth_lookup = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/depth_percentile_lookup.pkl"
    threads: 3
    resources:
        mem_mb=6000,
        time="00:30:00"
    log:
        "logs/calculate_depth_lookup.{sample}.log"
    run:
        import subprocess, pickle, sys

        with open(log[0], "w") as log_f:
        # -b limits the calculation to regions in your BED file
        # -a ensures we still output '0' values for target bases that were missed (crucial for accurate stats)
            cmd = f"samtools depth -a -b {input.bed} {input.bam} | awk '{{count[$3]++}} END {{for (d in count) print d, count[d]}}'"
            
            depth_counts = {}
            total_bases = 0
            
            # Stream the output of the command
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, text=True)
            
            for line in process.stdout:
                log_f.write(line)
                try:
                    depth, count = map(int, line.strip().split())
                    depth_counts[depth] = count
                    total_bases += count
                except ValueError:
                    log_f.write(f"Error parsing line: {line}")
                    continue
            
            process.stdout.close()
            retcode = process.wait()
            if retcode != 0:
                raise RuntimeError(f"samtools depth failed, see log: {log[0]}")

            # 2. Convert Histogram to Percentile Lookup
            sorted_depths = sorted(depth_counts.keys())
            cumulative_count = 0
            depth_to_percentile = {}

            for depth in sorted_depths:
                cumulative_count += depth_counts[depth]
                percentile = int((cumulative_count / total_bases) * 100)
                depth_to_percentile[depth] = min(percentile, 100) 

            # 3. Save as Pickle
            with open(output.depth_lookup, "wb") as f:
                pickle.dump(depth_to_percentile, f)


###########################################
# SNV calling
###########################################
#Step 2.1: PROCESS THE DCS ###
# Using the scattered BED chunks to parallelize variant calling at the DCS level 

rule mutect2_dcs_scatter_snv:
    input:
        ref = REF,
        bam = lambda wc: config["samples"][wc.sample]["dcs_bam"],
        bed_chunk = lambda wc: f"{BED_CHUNKS_DIR}/{wc.dcs_chunk}.bed"
    output:
        raw_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_mutect2.{{dcs_chunk}}.vcf.gz")
    threads: 2
    resources:
        mem_mb=6000,
        time="12:00:00",
        tmpdir=f"{PROCESSED_DATA_DIR}/variants/tmp"
    log:
        "logs/mutect2_dcs_scatter_snv.{sample}.{dcs_chunk}.log"
    shell:
        r"""
        (
        ulimit -n 8192
        mkdir -p {resources.tmpdir}

        # Define temp file path
        TEMP_VCF="{resources.tmpdir}/raw_mutect2.{wildcards.sample}.{wildcards.dcs_chunk}.tmp.vcf"

        gatk --java-options "-Xmx5G -Djava.io.tmpdir={resources.tmpdir}" Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed_chunk} \
            --native-pair-hmm-threads {threads} \
            -O $TEMP_VCF

        # Compress to output
        bgzip -c $TEMP_VCF > {output.raw_vcf}
        bcftools index -t {output.raw_vcf}

        # Clean up the uncompressed intermediate file
        rm $TEMP_VCF
        ) &> {log}
        """

rule gather_mutect2_raw:
    input:
        raw_vcf=lambda wildcards: expand(f"{PROCESSED_DATA_DIR}/variants/{wildcards.sample}/raw_mutect2.{{dcs_chunk}}.vcf.gz", dcs_chunk=DCS_CHUNK_IDS)
    output:
        gathered_raw_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_mutect2_gathered.vcf.gz")
    threads: 1  # bcftools concat is mostly single-threaded I/O bound
    resources:
        mem_mb=8000,
        time="12:00:00"
    log:
        "logs/gather_mutect2_raw.{sample}.log"
    shell:
        """
        (
        bcftools concat -a -D -O z --threads {threads} -o {output.gathered_raw_vcf} {input.raw_vcf}
        bcftools index -t {output.gathered_raw_vcf}
        ) &> {log}
        """ 

rule extract_indels_for_snv_processing:
    input:
        # Input VCF: The gathered raw Mutect2 VCF (SNVs + Indels)
        vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_mutect2_gathered.vcf.gz" 
    output:
        # Temporary output: The Indel VCF needed for distance calculation
        indel_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}_mutect2_indels_only.vcf.gz"),
        idx = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}_mutect2_indels_only.vcf.gz.tbi")
    threads: 1  
    resources:
        mem_mb=4000,
        time="02:00:00"
    log:
        "logs/extract_indels_for_snv_processing.{sample}.log"
    shell:
        r"""
        (
        # Select only INDELs/MNPs/Complex, sort, compress, and index.
        # This VCF contains only the POSITIONS of indels for distance calculation.
        bcftools view -v indels,mnps,other {input.vcf} | \
        bcftools sort -Oz -o {output.indel_vcf}
        bcftools index -t {output.indel_vcf}
        ) &> {log}
        """ 

rule precalculate_dcs_snv_chunk_data:
    """
    (SCATTER) Pre-calculates global data structures for a single BED chunk.
    This replaces the monolithic pre-calculation rule.
    """
    input:
        dcs_bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        ref = REF,
        bed_chunk = lambda wc: f"{BED_CHUNKS_DIR}/{wc.dcs_chunk}.bed"
    output:
        # Temporary output, one per chunk
        global_data_chunk = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/chunk_data/{{dcs_chunk}}.pkl")
    threads: 3
    resources:
        mem_mb=6000,
        time="9:00:00"
    log:
        "logs/precalculate_dcs_snv_chunk_data.{sample}.{dcs_chunk}.log"
    run:
        import pysam, pandas as pd, pickle, sys, contextlib
        from intervaltree import IntervalTree
        from collections import defaultdict
        
        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):
            # --- Helper functions (identical to your original rule) ---
            def build_interval_trees(bed_path):
                bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0,1,2],
                                names=["chrom","start","end"], dtype={"chrom": str})
                trees = {}
                for chrom in bed['chrom'].unique():
                    tree = IntervalTree()
                    for _, row in bed[bed['chrom'] == chrom].iterrows():
                        tree.addi(int(row.start), int(row.end)+1, True)
                    trees[chrom] = tree
                return trees

            def get_variant_positions_and_vafs_bed(samfile, bed_tree, ref_fasta, min_alt_count=1):
                positions, vafs = defaultdict(list), defaultdict(list)
                for chrom in bed_tree.keys():
                    for interval in bed_tree[chrom]:
                        for pileupcolumn in samfile.pileup(chrom, int(interval.begin), int(interval.end), truncate=True):
                            pos = pileupcolumn.reference_pos + 1
                            ref_base = ref_fasta.fetch(str(chrom), pos-1, pos).upper()
                            alt_count = total_count = 0
                            for pr in pileupcolumn.pileups:
                                if pr.is_del or pr.is_refskip or pr.query_position is None: continue
                                total_count += 1
                                base = pr.alignment.query_sequence[pr.query_position].upper()
                                if base != ref_base: alt_count += 1
                            if total_count == 0 or alt_count < min_alt_count: continue
                            positions[chrom].append(pos)
                            vafs[chrom].append(alt_count/total_count)
                return positions, vafs

            def build_read_to_variant_positions(samfile, bed_tree, ref_fasta):
                tmp = defaultdict(lambda: defaultdict(set))
                for chrom in bed_tree.keys():
                    for interval in bed_tree[chrom]:
                        for pileupcolumn in samfile.pileup(chrom, int(interval.begin), int(interval.end), truncate=True):
                            pos = pileupcolumn.reference_pos + 1
                            ref_base = ref_fasta.fetch(str(chrom), pos-1, pos).upper()
                            for pr in pileupcolumn.pileups:
                                if pr.is_del or pr.is_refskip or pr.query_position is None: 
                                    continue
                                base = pr.alignment.query_sequence[pr.query_position].upper()
                                if base != ref_base:
                                    tmp[chrom][pr.alignment.query_name].add(pos)
                read_to_variants = defaultdict(dict)
                for chrom in tmp:
                    for rname, posset in tmp[chrom].items():
                        read_to_variants[chrom][rname] = sorted(posset)
                return read_to_variants

            # --- MAIN EXECUTION (runs only on the chunk) ---
            
            # This now builds an interval tree for just the single input chunk
            print(f"Building interval trees for {input.bed_chunk}...")
            trees = build_interval_trees(input.bed_chunk) 
            
            samfile = pysam.AlignmentFile(input.dcs_bam, "rb")
            ref_fasta = pysam.FastaFile(input.ref)

            # 1. Get VAFs/Positions for this chunk
            # Note: We do NOT sort/unique them here. We do that once in the gather step.
            print("Extracting VAFs/variant positions...")
            chrom_variant_positions, chrom_variant_vafs = get_variant_positions_and_vafs_bed(samfile, trees, ref_fasta)

            # 2. Build the read-to-variant map for this chunk
            print("Building read-to-variant map...")
            read_to_variants = build_read_to_variant_positions(samfile, trees, ref_fasta)

            # 3. Save the results for this chunk
            print(f"Saving results to {output.global_data_chunk}...")
            chunk_data = {
                "chrom_variant_positions": chrom_variant_positions,
                "chrom_variant_vafs": chrom_variant_vafs,
                "read_to_variants": read_to_variants
            }
            with open(output.global_data_chunk, "wb") as f:
                pickle.dump(chunk_data, f)

rule gather_dcs_snv_global_data:
    """
    (GATHER) Merges the pre-calculated chunk data into a single global data file.
    """
    input:
        # This will collect all chunk.pkl files for this {sample}
        chunks = lambda wc: expand(f"{PROCESSED_DATA_DIR}/variants/{wc.sample}/chunk_data/{{dcs_chunk}}.pkl", dcs_chunk=DCS_CHUNK_IDS)
    output:
        # The final, non-temporary file the rest of the workflow expects
        global_data = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/dcs_global_data.pkl"
    threads: 1
    resources:
        mem_mb=16000,
        time="02:00:00"
    log:
        "logs/gather_dcs_snv_global_data.{sample}.log"
    run:
        import pickle, sys, contextlib
        from collections import defaultdict

        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):

            print("Starting merge of chunk data...")
            # These will hold the merged data from all chunks
            final_positions = defaultdict(list)
            final_vafs = defaultdict(list)
            
            # This map needs special merging to combine variant sets for the same read
            final_read_map_sets = defaultdict(lambda: defaultdict(set))

            # Loop over all chunk pickles and load/merge their data
            for chunk_file in input.chunks:
                with open(chunk_file, "rb") as f:
                    chunk_data = pickle.load(f)
                
                # 1. Merge positions and VAFs (simple extend)
                for chrom, pos_list in chunk_data["chrom_variant_positions"].items():
                    final_positions[chrom].extend(pos_list)
                    final_vafs[chrom].extend(chunk_data["chrom_variant_vafs"][chrom])

                # 2. Merge read-to-variant map
                # We use sets to handle reads that span chunks
                for chrom, reads in chunk_data["read_to_variants"].items():
                    for read_name, pos_list in reads.items():
                        final_read_map_sets[chrom][read_name].update(pos_list)

            # --- Post-process the merged data ---
            print("Post-processing merged data...")

            # 1. Clean and sort the VAF/Position pairs (as in your original rule)
            cleaned_positions = defaultdict(list)
            cleaned_vafs = defaultdict(list)
            for chrom in final_positions.keys():
                # Apply the sorting and unique-ing logic *after* gathering
                pos_vaf_pairs = sorted(set(zip(final_positions[chrom], final_vafs[chrom])))
                cleaned_positions[chrom] = [p for p,v in pos_vaf_pairs]
                cleaned_vafs[chrom] = [v for p,v in pos_vaf_pairs]
                
            # 2. Convert the read map sets back to sorted lists
            final_read_map_sorted = defaultdict(dict)
            for chrom, reads in final_read_map_sets.items():
                for rname, posset in reads.items():
                    final_read_map_sorted[chrom][rname] = sorted(posset)

            # 3. Save the final, merged data
            final_global_data = {
                "chrom_variant_positions": cleaned_positions,
                "chrom_variant_vafs": cleaned_vafs,
                "read_to_variants": final_read_map_sorted
            }
            
            print(f"Saving final merged data to {output.global_data}...")
            with open(output.global_data, "wb") as f:
                pickle.dump(final_global_data, f)
            
            print("Merge complete.")

rule make_vcf_snv_scatter:
    input:
        global_data  = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/dcs_global_data.pkl", 
        depth_lookup = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/depth_percentile_lookup.pkl", 
        dcs_bam      = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        bed_chunk    = lambda wildcards: f"{BED_CHUNKS_DIR}/{wildcards.dcs_chunk}.bed",
        ref          = REF,
        encode       = config["PATHS"]["encode_blacklist"],
        repeatmasker = config["PATHS"]["repeatmasker_blacklist"],
        indel_vcf    = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}_mutect2_indels_only.vcf.gz"
    output:
        feature_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/feature_rich.{{dcs_chunk}}.vcf")
    threads: 4
    resources:
        mem_mb=16000,
        time="10:00:00"
    log:
        "logs/make_vcf_smv_scatter.{sample}.{dcs_chunk}.log"
    run:
        import pysam, pandas as pd, math, statistics, bisect, gzip, pickle
        from intervaltree import IntervalTree
        from collections import defaultdict, Counter
        import contextlib

        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):

            print("Starting feature extraction for chunk:", input.bed_chunk)

            # --- HELPERS ---
            def build_interval_trees(bed_path):
                bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0,1,2],
                                names=["chrom","start","end"], dtype={"chrom": str})
                trees = {}
                for chrom in bed['chrom'].unique():
                    tree = IntervalTree()
                    for _, row in bed[bed['chrom'] == chrom].iterrows():
                        tree.addi(int(row.start), int(row.end)+1, True)
                    trees[chrom] = tree
                return trees

            def in_bed(chrom, pos, trees):
                return chrom in trees and trees[chrom].overlaps(pos)

            def get_homopolymer_length(ref_fasta, chrom, pos):
                base = ref_fasta.fetch(str(chrom), pos-1, pos).upper()
                seq_start = max(pos-20, 0)
                seq_end = pos+20
                seq = ref_fasta.fetch(str(chrom), seq_start, seq_end).upper()
                center = pos - seq_start - 1
                left = 0
                for i in range(center-1, -1, -1):
                    if seq[i] == base: left += 1
                    else: break
                right = 0
                for i in range(center+1, len(seq)):
                    if seq[i] == base: right += 1
                    else: break
                return left + 1 + right

            def calc_shannon_entropy(ref_fasta, chrom, pos, window=40):
                half_window = window // 2
                start = max(pos - 1 - half_window, 0)
                end = pos - 1 + half_window
                seq = ref_fasta.fetch(str(chrom), start, end).upper()
                if not seq: return 0.0
                counts = Counter(seq)
                total = len(seq)
                entropy = 0.0
                for base, count in counts.items():
                    p = count / total
                    entropy -= p * math.log2(p)
                return entropy

            def count_nearby_variants(chrom, pos, variant_positions_dict, distance):
                if chrom not in variant_positions_dict: return 0
                positions = variant_positions_dict[chrom]
                left = bisect.bisect_left(positions, pos - distance)
                right = bisect.bisect_right(positions, pos + distance)
                count = 0
                for i in range(left, right):
                    if positions[i] != pos:
                        count += 1
                return count

            def get_softclip_metrics(read, query_pos):
                cigar = read.cigartuples
                if not cigar: return 0, -1
                total_clips = 0
                left_clip_len = 0
                right_clip_len = 0
                if cigar[0][0] == 4:
                    left_clip_len = cigar[0][1]
                    total_clips += left_clip_len
                if cigar[-1][0] == 4:
                    right_clip_len = cigar[-1][1]
                    total_clips += right_clip_len
                if total_clips == 0: return 0, 1000
                
                dist_left = 1000
                dist_right = 1000
                if left_clip_len > 0: dist_left = max(0, query_pos - left_clip_len)
                if right_clip_len > 0:
                    rlen = read.query_length
                    boundary = rlen - right_clip_len
                    dist_right = max(0, boundary - query_pos)
                
                valid_dists = []
                if left_clip_len > 0: valid_dists.append(dist_left)
                if right_clip_len > 0: valid_dists.append(dist_right)
                return total_clips, min(valid_dists)

            def build_indel_positions(indel_vcf_path):
                indel_positions = defaultdict(list)
                open_func = gzip.open if indel_vcf_path.endswith(".gz") else open
                with open_func(indel_vcf_path, "rt") as f:
                    for line in f:
                        if line.startswith("#"): continue
                        fields = line.strip().split("\t")
                        chrom, pos = fields[0], int(fields[1])
                        indel_positions[chrom].append(pos)
                for chrom in indel_positions: indel_positions[chrom].sort()
                return indel_positions

            def distance_to_closest_indel(chrom, pos, indel_positions):
                if chrom not in indel_positions or not indel_positions[chrom]: return -1
                positions = indel_positions[chrom]
                idx = bisect.bisect_left(positions, pos)
                distances = []
                if idx < len(positions): distances.append(abs(positions[idx]-pos))
                if idx > 0: distances.append(abs(positions[idx-1]-pos))
                return min(distances) if distances else -1

            def get_best_alt_allele(pileupcolumn, ref_base):
                alt_counts = defaultdict(int)
                total_reads = 0
                for pr in pileupcolumn.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None: continue
                    base = pr.alignment.query_sequence[pr.query_position].upper()
                    total_reads += 1
                    if base != ref_base and base in ['A','C','G','T']:
                        alt_counts[base] += 1
                if not alt_counts: return None, 0, total_reads
                best_alt = max(alt_counts, key=alt_counts.get)
                return best_alt, alt_counts[best_alt], total_reads

            # --- LOAD DATA ---
            print("Loading global data...")
            with open(input.global_data, "rb") as f:
                global_data = pickle.load(f)
                chrom_variant_positions = global_data["chrom_variant_positions"]
                read_to_variants = global_data["read_to_variants"]
            
            print("Loading depth lookup...")
            with open(input.depth_lookup, "rb") as f:
                depth_to_percentile = pickle.load(f)
            
            def get_depth_percentile(depth, lookup):
                if depth in lookup: return lookup[depth]
                if not lookup: return 0
                if depth > max(lookup.keys()): return 100
                return lookup.get(depth, 100)

            encode_blacklist = build_interval_trees(input.encode)
            repeatmasker_blacklist = build_interval_trees(input.repeatmasker)
            indel_positions = build_indel_positions(input.indel_vcf)
            samfile = pysam.AlignmentFile(input.dcs_bam, "rb")
            ref_fasta = pysam.FastaFile(input.ref)
            bed_chunk_trees = build_interval_trees(input.bed_chunk) 

            regions_to_scan = []
            for chrom in bed_chunk_trees.keys():
                for interval in bed_chunk_trees[chrom]:
                    regions_to_scan.append((chrom, interval.begin, interval.end))

            # --- WRITE OUTPUT ---
            print(f"Writing VCF to {output.feature_vcf}...")
            with open(output.feature_vcf, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write(f"##consensus_bam={input.dcs_bam}\n")
                
                # --- INFO HEADERS ---
                f.write('##INFO=<ID=N_TOTAL,Number=1,Type=Integer,Description="Total reads covering position">\n')
                f.write('##INFO=<ID=N_ALT,Number=1,Type=Integer,Description="Reads supporting alt base">\n')
                f.write('##INFO=<ID=DEPTH_PERCENTILE,Number=1,Type=Integer,Description="Genomic depth percentile (0-100) of this locus">\n')
                
                # Base & Read Quality
                f.write('##INFO=<ID=MEAN_ALT_BASEQ,Number=1,Type=Float,Description="Mean base quality of the specific variant base">\n')
                f.write('##INFO=<ID=SD_ALT_BASEQ,Number=1,Type=Float,Description="SD of base quality of the specific variant base">\n')
                f.write('##INFO=<ID=AVG_ALT_READ_MEAN_BQ,Number=1,Type=Float,Description="Average of the mean base qualities of the whole reads supporting alt">\n')
                f.write('##INFO=<ID=MEAN_ALT_MAPQ,Number=1,Type=Float,Description="Mean mapping quality of alt-supporting reads">\n')
                f.write('##INFO=<ID=AVG_ALT_ASXS,Number=1,Type=Float,Description="Average ASXS score across alt reads">\n')
                f.write('##INFO=<ID=MEAN_REF_MAPQ,Number=1,Type=Float,Description="Mean mapping quality of non-alt reads">\n')
                f.write('##INFO=<ID=AVG_REF_ASXS,Number=1,Type=Float,Description="Average ASXS score across non-alt reads">\n')

                # Read Hygiene & Variants
                f.write('##INFO=<ID=AVG_ALT_NM,Number=1,Type=Float,Description="Average NM tag (edit distance) of alt reads">\n')
                f.write('##INFO=<ID=AVG_ALT_N_COUNT,Number=1,Type=Float,Description="Average count of N bases in alt reads">\n')
                f.write('##INFO=<ID=AVG_ALT_SOFTCLIP,Number=1,Type=Float,Description="Average number of softclipped bases in alt reads">\n')
                f.write('##INFO=<ID=MEDIAN_DIST_TO_SOFTCLIP,Number=1,Type=Float,Description="Median distance of variant to nearest softclipped edge (1000 if none)">\n')
                
                # --- RESTORED VARIANT FEATURES ---
                f.write('##INFO=<ID=AVG_ALTREAD_VARIANTS,Number=1,Type=Float,Description="Average number of other variants on alt reads">\n')
                f.write('##INFO=<ID=AVG_NONALTREAD_VARIANTS,Number=1,Type=Float,Description="Average number of variant positions on reads NOT supporting the alt allele">\n')
                
                # Geometry & Context
                f.write('##INFO=<ID=MEAN_INSERT,Number=1,Type=Float,Description="Mean insert size of alt-supporting reads">\n')
                f.write('##INFO=<ID=SD_INSERT,Number=1,Type=Float,Description="SD of insert size of alt-supporting reads">\n')
                f.write('##INFO=<ID=AVG_READ_POSITION,Number=1,Type=Float,Description="Average absolute position of the variant in alt-supporting reads (1-based)">\n')
                
                # --- RESTORED LENGTH FEATURE ---
                f.write('##INFO=<ID=AVG_READ_LENGTH,Number=1,Type=Float,Description="Average length of alt-supporting reads">\n')
                
                f.write('##INFO=<ID=REF_ENTROPY_40BP,Number=1,Type=Float,Description="Shannon entropy of reference sequence in +/- 20bp window">\n')
                f.write('##INFO=<ID=VARIANTS_20BP,Number=1,Type=Integer,Description="Count of other variants within 20bp">\n')
                f.write('##INFO=<ID=VARIANTS_250BP,Number=1,Type=Integer,Description="Count of other variants within 250bp">\n')
                f.write('##INFO=<ID=HOMOPOLYMER_LEN,Number=1,Type=Integer,Description="Homopolymer length at variant position">\n')
                f.write('##INFO=<ID=GC_CONTENT,Number=1,Type=Float,Description="GC content in +/- 20bp">\n')
                f.write('##INFO=<ID=DIST_TO_INDEL,Number=1,Type=Integer,Description="Distance to closest indel">\n')
                f.write('##INFO=<ID=IN_ENCODE_BLACKLIST,Number=0,Type=Flag,Description="Variant lies within ENCODE blacklist region">\n')
                f.write('##INFO=<ID=IN_REPEATMASKER,Number=0,Type=Flag,Description="Variant lies within RepeatMasker region">\n')
                f.write('##INFO=<ID=PAIRED_VARIANT,Number=0,Type=Flag,Description="Alt allele is on reads that also contain another variant within 30bp">\n')

                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                for chrom, start, end in regions_to_scan:
                    for pileupcolumn in samfile.pileup(chrom, start, end, truncate = True):
                        pos = pileupcolumn.reference_pos + 1
                        try:
                            ref_base = ref_fasta.fetch(str(chrom), pos-1, pos).upper()
                        except:
                            continue 
                        
                        alt_base, n_alt, n_total = get_best_alt_allele(pileupcolumn, ref_base)
                        
                        if alt_base is None or n_alt == 0: 
                            continue
                        
                        # --- METRIC COLLECTIONS ---
                        alt_base_quals_specific = [] 
                        alt_read_mean_bqs = []       
                        alt_mapqs = []
                        alt_inserts = []
                        alt_asxs = []
                        alt_nms = []
                        alt_n_counts = []
                        alt_softclip_amts = []
                        alt_dists_to_clip = []
                        alt_read_positions = []
                        alt_reads_objs = [] 
                        alt_read_lengths = [] # RESTORED

                        ref_mapqs = []
                        ref_asxs = []
                        ref_reads_objs = [] # RESTORED FOR NON-ALT VARIANT COUNTING

                        for pr in pileupcolumn.pileups:
                            if pr.is_del or pr.is_refskip or pr.query_position is None: continue
                            
                            read = pr.alignment
                            read_base = read.query_sequence[pr.query_position].upper()
                            
                            try: as_score = read.get_tag("AS")
                            except KeyError: as_score = 0
                            
                            if read_base == alt_base:
                                alt_base_quals_specific.append(read.query_qualities[pr.query_position])
                                
                                if read.query_qualities:
                                    mean_whole_read_bq = sum(read.query_qualities) / len(read.query_qualities)
                                    alt_read_mean_bqs.append(mean_whole_read_bq)
                                
                                alt_mapqs.append(read.mapping_quality)
                                alt_inserts.append(abs(read.template_length))
                                alt_asxs.append(as_score)
                                
                                try: alt_nms.append(read.get_tag("NM"))
                                except KeyError: alt_nms.append(0)
                                
                                alt_n_counts.append(read.query_sequence.count("N"))
                                
                                clip_amt, clip_dist = get_softclip_metrics(read, pr.query_position)
                                alt_softclip_amts.append(clip_amt)
                                if clip_dist < 999: alt_dists_to_clip.append(clip_dist)

                                alt_read_positions.append(pr.query_position + 1)
                                alt_read_lengths.append(read.query_length) # RESTORED
                                alt_reads_objs.append(read)
                                
                            elif read_base == ref_base:
                                ref_mapqs.append(read.mapping_quality)
                                ref_asxs.append(as_score)
                                ref_reads_objs.append(read) # RESTORED

                        # --- AGGREGATION ---
                        mean_alt_bq_spec = statistics.mean(alt_base_quals_specific) if alt_base_quals_specific else 0.0
                        sd_alt_bq_spec = statistics.stdev(alt_base_quals_specific) if len(alt_base_quals_specific) > 1 else 0.0
                        avg_alt_read_mean_bq = statistics.mean(alt_read_mean_bqs) if alt_read_mean_bqs else 0.0

                        mean_alt_mapq = statistics.mean(alt_mapqs) if alt_mapqs else 0.0
                        avg_alt_asxs = statistics.mean(alt_asxs) if alt_asxs else 0.0
                        mean_insert = statistics.mean(alt_inserts) if alt_inserts else 0.0
                        sd_insert = statistics.stdev(alt_inserts) if len(alt_inserts) > 1 else 0.0
                        
                        avg_alt_nm = statistics.mean(alt_nms) if alt_nms else 0.0
                        avg_alt_n = statistics.mean(alt_n_counts) if alt_n_counts else 0.0
                        avg_alt_softclip = statistics.mean(alt_softclip_amts) if alt_softclip_amts else 0.0
                        median_dist_clip = statistics.median(alt_dists_to_clip) if alt_dists_to_clip else 1000.0
                        
                        avg_read_pos = statistics.mean(alt_read_positions) if alt_read_positions else 0.0
                        avg_read_len = statistics.mean(alt_read_lengths) if alt_read_lengths else 0.0 # RESTORED

                        mean_ref_mapq = statistics.mean(ref_mapqs) if ref_mapqs else 0.0
                        avg_ref_asxs = statistics.mean(ref_asxs) if ref_asxs else 0.0

                        # Context
                        ref_entropy = calc_shannon_entropy(ref_fasta, chrom, pos)
                        vars_20 = count_nearby_variants(chrom, pos, chrom_variant_positions, 20)
                        vars_250 = count_nearby_variants(chrom, pos, chrom_variant_positions, 250)
                        homopolymer = get_homopolymer_length(ref_fasta, chrom, pos)
                        
                        start_gc, end_gc = max(pos-21, 0), pos+20
                        seq_gc = ref_fasta.fetch(str(chrom), start_gc, end_gc).upper()
                        gc_content = (seq_gc.count("G")+seq_gc.count("C"))/len(seq_gc) if seq_gc else 0.0
                        
                        dist_to_indel = distance_to_closest_indel(chrom, pos, indel_positions)
                        
                        # Read variant counts (ALT)
                        avg_altread_vars = 0.0
                        paired_flag = False
                        if alt_reads_objs:
                            var_counts = []
                            for r in alt_reads_objs:
                                vps = read_to_variants.get(chrom, {}).get(r.query_name, [])
                                var_counts.append(len(vps))
                                if not paired_flag:
                                    if any(abs(vp - pos) <= 30 and vp != pos for vp in vps):
                                        paired_flag = True
                            avg_altread_vars = statistics.mean(var_counts)
                        
                        # Read variant counts (REF) -- RESTORED
                        avg_nonaltread_vars = 0.0
                        if ref_reads_objs:
                            var_counts_ref = []
                            for r in ref_reads_objs:
                                vps = read_to_variants.get(chrom, {}).get(r.query_name, [])
                                var_counts_ref.append(len(vps))
                            avg_nonaltread_vars = statistics.mean(var_counts_ref)

                        # Lookups
                        percentile = get_depth_percentile(n_total, depth_to_percentile)
                        in_enc = in_bed(chrom, pos, encode_blacklist)
                        in_rep = in_bed(chrom, pos, repeatmasker_blacklist)
                        
                        flags_str = []
                        if in_enc: flags_str.append("IN_ENCODE_BLACKLIST")
                        if in_rep: flags_str.append("IN_REPEATMASKER")
                        if paired_flag: flags_str.append("PAIRED_VARIANT")
                        flags_out = ";" + ";".join(flags_str) if flags_str else ""

                        info = (
                            f"N_TOTAL={n_total};N_ALT={n_alt};DEPTH_PERCENTILE={percentile};"
                            f"MEAN_ALT_BASEQ={mean_alt_bq_spec:.2f};SD_ALT_BASEQ={sd_alt_bq_spec:.2f};"
                            f"AVG_ALT_READ_MEAN_BQ={avg_alt_read_mean_bq:.2f};"
                            f"MEAN_ALT_MAPQ={mean_alt_mapq:.2f};AVG_ALT_ASXS={avg_alt_asxs:.2f};"
                            f"MEAN_REF_MAPQ={mean_ref_mapq:.2f};AVG_REF_ASXS={avg_ref_asxs:.2f};"
                            f"AVG_ALT_NM={avg_alt_nm:.2f};AVG_ALT_N_COUNT={avg_alt_n:.2f};"
                            f"AVG_ALT_SOFTCLIP={avg_alt_softclip:.2f};MEDIAN_DIST_TO_SOFTCLIP={median_dist_clip:.1f};"
                            f"AVG_ALTREAD_VARIANTS={avg_altread_vars:.2f};AVG_NONALTREAD_VARIANTS={avg_nonaltread_vars:.2f};"
                            f"MEAN_INSERT={mean_insert:.1f};SD_INSERT={sd_insert:.1f};"
                            f"AVG_READ_POSITION={avg_read_pos:.1f};AVG_READ_LENGTH={avg_read_len:.1f};"
                            f"REF_ENTROPY_40BP={ref_entropy:.3f};VARIANTS_20BP={vars_20};VARIANTS_250BP={vars_250};"
                            f"HOMOPOLYMER_LEN={homopolymer};GC_CONTENT={gc_content:.2f};DIST_TO_INDEL={dist_to_indel}"
                            f"{flags_out}"
                        )

                        f.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\t.\t{info}\n")
                        
                print("Feature extraction complete!")

# Step 2.2: GATHER THE CHUNKED VCFs INTO A SINGLE VCF PER SAMPLE
# The dcs info has been called in parallel, now we need to gather the scattered VCFs into a single VCF per sample

rule clean_dcs_snv_scatter:
    input:
        feature_vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/feature_rich.{{dcs_chunk}}.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # Temporary output that has a clean header, ready for concatenation
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/cleaned_scatter.{{dcs_chunk}}.vcf.gz"),
        tbi = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/cleaned_scatter.{{dcs_chunk}}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=6000,
        time="06:00:00"
    log:
        "logs/clean_dcs_snv_scatter.{sample}.{dcs_chunk}.log"
    shell:
        r"""
        (
        set -euo pipefail
        mkdir -p $(dirname {output.vcf_gz})

        # --- Step 1: Annotate Contig Header FIRST ---
        # This reads the header-deficient input file and adds the
        # contig lines from the contigs_hdr file.
        # The output is piped to the rest of the chain.
        bcftools annotate -h {input.contigs_hdr} {input.feature_vcf} | \

        # --- Step 2: Normalize ---
        # This command now receives a VCF with a valid header.
        bcftools norm -m - -f {input.ref} - | \
        
        # --- Step 3: Drop Genotypes ---
        bcftools view -G -O v | \
        
        # --- Step 4: Clean FILTER column (Awk) ---
        awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{$7="."; print}}' | \
        
        # --- Step 5: Sort, Compress, and Index ---
        bcftools sort -Oz -o {output.vcf_gz}
        
        # --- Step 6: Index ---
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """

rule gather_snv_vcf:
    input:
        scattered_vcf = lambda wc: expand(f"{PROCESSED_DATA_DIR}/variants/{wc.sample}/cleaned_scatter.{{dcs_chunk}}.vcf.gz", dcs_chunk=DCS_CHUNK_IDS),
        scattered_tbi = lambda wc: expand(f"{PROCESSED_DATA_DIR}/variants/{wc.sample}/cleaned_scatter.{{dcs_chunk}}.vcf.gz.tbi", dcs_chunk=DCS_CHUNK_IDS)
    output:
        gathered_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/full_feature_rich.vcf")
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00",
        partition="icelake"
    log:
        "logs/gather_snv_vcf.{sample}.log"
    shell:
        """
        (
        # Ensure the VCFs are concatenated correctly, preserving header order.
        # -a ensures that the header only comes from the first file.
        # -D suppresses duplicate INFO/FORMAT fields, useful if headers were slightly different.
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        ) &> {log}
        """

# The final list of DCS candidates will guide the SSCS extraction.


rule clean_dcs_snv_candidates:
    input:
        # Input VCF: Contains all features gathered from parallel runs (no filtering applied yet)
        gathered_vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/full_feature_rich.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # The output is the final, indexed guide VCF for SSCS
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_dcs_candidates.vcf.gz"),
        idx = temp(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_dcs_candidates.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/clean_dcs_snv_candidates.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail
        mkdir -p $(dirname {output.vcf_gz})
        # --- Step 1: Normalize and Drop Genotypes (-G) ---
        # Normalize/Left-align variants (-m -f) and drop large FORMAT/SAMPLE columns (-G)
        # The output of this step is a clean VCF without sample data.
        bcftools norm -m - -f {input.ref} {input.gathered_vcf} | \
        bcftools view -G -O v | \
        
        # --- Step 2: Clean FILTER column (Awk) ---
        # Replicates your Python script's core function: setting $7 (FILTER column) to '.'
        awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{$7="."; print}}' | \
        
        # --- Step 3: Annotate Contig Header ---
        # Re-introduces all necessary ##contig lines from your reference header file.
        bcftools annotate -h {input.contigs_hdr} -O v | \
        
        # --- Step 4: Sort, Compress, and Index ---
        # Sorts the VCF, compresses it with bgzip (-Oz), and pipes to the final output.
        bcftools sort -Oz -o {output.vcf_gz}
        
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """ 

# STEP 2.3: COLLECTING SSCS DATA USING THE DCS CANDIDATES GUIDE VCF


rule split_snv_candidate_vcf:
    input:
        dcs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_dcs_candidates.vcf.gz"
    output:
       chunks=temp(expand("{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}/chunk_{i}.vcf.gz", PROCESSED_DATA_DIR=PROCESSED_DATA_DIR, i=range(1, SSCS_CHUNKS + 1)))
    params:
        n = SSCS_CHUNKS,
        outdir = f"{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}",
        split_suffix = 2
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/split_snv_candidate_vcf.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail

        mkdir -p {params.outdir}

        bcftools view -h {input.dcs_vcf} > {params.outdir}/header.vcf

        DATA_LINES=$(bcftools view -H {input.dcs_vcf} | wc -l)

        if [ "$DATA_LINES" -eq 0 ]; then
            LINES_PER_CHUNK=1
        else
            LINES_PER_CHUNK=$(( (DATA_LINES + {params.n} - 1) / {params.n} ))
        fi

        bcftools view -H {input.dcs_vcf} > {params.outdir}/all_variants.tmp

        for i in $(seq 1 {params.n}); do

            start=$(( (i - 1) * LINES_PER_CHUNK + 1 ))
            end=$(( i * LINES_PER_CHUNK ))

            output_file="{params.outdir}/chunk_${{i}}.vcf.gz"

            sed -n "${{start}},${{end}}p" {params.outdir}/all_variants.tmp > {params.outdir}/chunk_${{i}}.vcf

            cat {params.outdir}/header.vcf {params.outdir}/chunk_${{i}}.vcf | \
                bcftools sort -Oz -o $output_file

            tabix -p vcf $output_file

            rm {params.outdir}/chunk_${{i}}.vcf
        done

        rm {params.outdir}/header.vcf
        rm {params.outdir}/all_variants.tmp

        ) &> {log}
        """

rule snv_format_sscs_counts_scatter:
    input:
        # The small chunk of DCS candidates (the positions we need SSCS data for)
        candidate_vcf_chunk = f"{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}/{{sscs_chunk}}.vcf.gz",
        # The SSCS BAM file
        sscs_bam = lambda wildcards: config["samples"][wildcards.sample]["sscs_bam"],
        sscs_bai = lambda wildcards: config["samples"][wildcards.sample]["sscs_bam"] + ".bai",
    output:
        # Output VCF chunk containing SSCS counts and new metrics (e.g., RP_MED_SSCS)
        sscs_chunk_vcf = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/metrics.{{sscs_chunk}}.vcf")
    threads: 1
    resources:
        mem_mb=2000,
        time="05:00:00"
    log:
        "logs/snv_format_sscs_counts_scatter.{sample}.{sscs_chunk}.log"
    run:
        import pysam, statistics
        from collections import defaultdict
        import contextlib

        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):

            print("Starting SSCS metrics calculation for chunk:", input.candidate_vcf_chunk)

            # Helper function to calculate read-level metrics
            def calculate_sscs_metrics(pileupcolumn, ref_base, alt_base):
                rel_positions, alt_fwd, alt_rev = [], 0, 0
                n_alt = 0
                
                for pr in pileupcolumn.pileups:
                    # Filter out deletions, reference skips, and reads without a query position
                    if pr.is_del or pr.is_refskip or pr.query_position is None: continue
                    
                    read_base = pr.alignment.query_sequence[pr.query_position].upper()
                    
                    if read_base == alt_base: 
                        n_alt += 1
                        
                        # 1. Read Position (RP)
                        rel_positions.append(pr.query_position / pr.alignment.query_length)
                        
                        # 2. Strand Bias (SB)
                        if pr.alignment.is_reverse:
                            alt_rev += 1
                        else:
                            alt_fwd += 1
                
                n_total = len([pr for pr in pileupcolumn.pileups if not pr.is_del and not pr.is_refskip])
                
                # Calculate final metrics
                # Defaults provided if n_alt == 0 to avoid crashes
                rp_med = statistics.median(rel_positions) if rel_positions else 0.5
                rp_sd = statistics.stdev(rel_positions) if len(rel_positions) > 1 else 0.0
                
                total_alt_reads = alt_fwd + alt_rev
                strand_bias = alt_fwd / total_alt_reads if total_alt_reads > 0 else 0.5

                return n_alt, n_total, rp_med, rp_sd, strand_bias

            # --- Main Execution ---
            print("Opening BAM and VCF files...")
            sscs_bam = pysam.AlignmentFile(input.sscs_bam, "rb")
            vcf_in = pysam.VariantFile(input.candidate_vcf_chunk)
            
            try:
                print("Writing SSCS-enhanced VCF to", output.sscs_chunk_vcf)
                with open(output.sscs_chunk_vcf, "w") as fout:
                    # Prepare VCF Header
                    fout.write("##fileformat=VCFv4.2\n")
                    fout.write("##INFO=<ID=N_TOTAL_SSCS,Number=1,Type=Integer,Description=\"Total depth (SSCS) at position\">\n")
                    fout.write("##INFO=<ID=N_ALT_SSCS,Number=1,Type=Integer,Description=\"Alternate depth (SSCS) at position\">\n")
                    fout.write("##INFO=<ID=RP_MED_SSCS,Number=1,Type=Float,Description=\"Median read-relative position of alt base (SSCS)\">\n")
                    fout.write("##INFO=<ID=RP_SD_SSCS,Number=1,Type=Float,Description=\"SD of read-relative positions of alt base (SSCS)\">\n")
                    fout.write("##INFO=<ID=STRAND_BIAS_SSCS,Number=1,Type=Float,Description=\"Fraction of alt reads on forward strand (SSCS)\">\n")
                    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                    for record in vcf_in:
                        chrom = record.chrom
                        pos = record.pos 
                        ref_base = record.ref
                        
                        # Skip complex/multi-allelic variants
                        if len(record.alts) != 1:
                            continue
                        
                        alt_base = record.alts[0]
                        
                        # Fetch pileup for the single position
                        # truncate=True ensures we don't grab reads that merely overlap the window but not the specific base
                        try:
                            iter_pileup = sscs_bam.pileup(chrom, pos - 1, pos, truncate=True, min_mapping_quality=0, min_base_quality=0)
                            pileupcolumn = next(iter_pileup, None)
                        except ValueError:
                            pileupcolumn = None

                        if pileupcolumn is None:
                            # No coverage: Output 0 support explicitly
                            n_alt, n_total, rp_med, rp_sd, strand_bias = 0, 0, 0.5, 0.0, 0.5
                        else:
                            n_alt, n_total, rp_med, rp_sd, strand_bias = calculate_sscs_metrics(pileupcolumn, ref_base, alt_base)

                        # --- CRITICAL CHANGE ---
                        # We NO LONGER skip if n_alt == 0. 
                        # We write the line so the final VCF clearly records "0" SSCS support.
                        
                        new_info = (f"N_TOTAL_SSCS={n_total};N_ALT_SSCS={n_alt};"
                                    f"RP_MED_SSCS={rp_med:.4f};RP_SD_SSCS={rp_sd:.4f};"
                                    f"STRAND_BIAS_SSCS={strand_bias:.4f}")
                        
                        # Write the new VCF line
                        out_line = f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\t.\t{new_info}\n"
                        fout.write(out_line)

            finally:
                # Ensure file handles are closed safely
                sscs_bam.close()
                vcf_in.close()
                
            print("SSCS metrics calculation complete for chunk:", input.candidate_vcf_chunk)

rule snv_clean_sscs_scatter:
    input:
        feature_vcf = f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/metrics.{{sscs_chunk}}.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # Temporary output that has a clean header, ready for concatenation
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/cleaned_scatter.{{sscs_chunk}}.vcf.gz"),
        tbi = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/cleaned_scatter.{{sscs_chunk}}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=1000,
        time="06:00:00",
        partition="icelake"
    log:
        "logs/snv_clean_sscs_scatter.{sample}.{sscs_chunk}.log"
    shell:
        r"""
        (
        set -euo pipefail
        mkdir -p $(dirname {output.vcf_gz})

        # --- Step 1: Annotate Contig Header FIRST ---
        # This reads the header-deficient input file and adds the
        # contig lines from the contigs_hdr file.
        # The output is piped to the rest of the chain.
        bcftools annotate -h {input.contigs_hdr} {input.feature_vcf} | \

        # --- Step 2: Normalize ---
        # This command now receives a VCF with a valid header.
        bcftools norm -m - -f {input.ref} - | \
        
        # --- Step 3: Drop Genotypes ---
        bcftools view -G -O v | \
        
        # --- Step 4: Clean FILTER column (Awk) ---
        awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{$7="."; print}}' | \
        
        # --- Step 5: Sort, Compress, and Index ---
        bcftools sort -Oz -o {output.vcf_gz}
        
        # --- Step 6: Index ---
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """ 


rule snv_gather_sscs_counts:
    input:
        scattered_vcf = lambda wc: [f"{PROCESSED_DATA_DIR}/sscs_scatter/{wc.sample}/cleaned_scatter.{chunk}.vcf.gz" for chunk in SSCS_CHUNK_IDS],
        idx = lambda wc: [f"{PROCESSED_DATA_DIR}/sscs_scatter/{wc.sample}/cleaned_scatter.{chunk}.vcf.gz.tbi" for chunk in SSCS_CHUNK_IDS]
    output:
        gathered_vcf = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/full_sscs_metrics.vcf.gz")
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    log:
        "logs/snv_gather_sscs_counts.{sample}.log"
    shell:
        """
        (
        # Concatenate VCFs, ensuring headers are merged correctly.
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        ) &> {log}
        """

rule snv_clean_sscs_vcf:
    input:
        gathered_vcf = f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/full_sscs_metrics.vcf.gz",
        contigs_hdr = CONTIGS_HDR,
    output:
        # Final clean, compressed, and indexed SSCS VCF
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_sscs_metrics_clean.vcf.gz"),
        idx = temp(f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_sscs_metrics_clean.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/snv_clean_sscs_vcf.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail
        mkdir -p $(dirname {output.vcf_gz})
        # Annotate header (ensures correct contig lines) and sort
        bcftools annotate -h {input.contigs_hdr} {input.gathered_vcf} | \
        bcftools sort -Oz -o {output.vcf_gz}
        
        # Index the VCF
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """ 

# STEP 2.4: COMBINE USE THE SSCS DATA TO ANNOTATE THE DCS VCF, CREATING THE FINAL VARIANT LIST WITH ALL METRICS


rule snv_final_annotate:
    input:
        dcs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_dcs_candidates.vcf.gz",
        sscs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_sscs_metrics_clean.vcf.gz",
        dcs_vcf_idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_dcs_candidates.vcf.gz.tbi",
        sscs_vcf_idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_sscs_metrics_clean.vcf.gz.tbi"
    output:
      vcf_gz = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined.vcf.gz",
        idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined.vcf.gz.tbi"
    params:
        annotation_fields='INFO/N_ALT_SSCS,INFO/N_TOTAL_SSCS,INFO/STRAND_BIAS_SSCS,INFO/RP_MED_SSCS,INFO/RP_SD_SSCS'
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/snv_final_annotate.{sample}.log"
    shell:
        r"""
        (
        bcftools annotate \
          -a {input.sscs_vcf} \
          -c "{params.annotation_fields}" \
          -Oz -o {output.vcf_gz} \
          {input.dcs_vcf}
          
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """ 


###########################################
# Indel calling - this actually proceeds alongside SNV calling
###########################################

#Step 3.1: PROCESS THE DCS ###
# Using the scattered BED chunks to parallelize variant calling at the DCS level 

rule scan_raw_indel_signals:
    """
    (SCATTER) Scans the BAM chunk to build a 'Map of the Terrain'.
    Output: A pickle file containing:
      1. Indel Targets: Positions where we will attempt to call variants.
      2. Noise Maps: Sorted lists of messy positions (SNVs/Indels) for context metrics.
    """
    input:
        dcs_bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        bed_chunk = lambda wc: f"{BED_CHUNKS_DIR}/{wc.dcs_chunk}.bed",
        ref = REF
    output:
        scan_data = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_signals/{{dcs_chunk}}.pkl")
    threads: 2
    resources:
        mem_mb=6000,
        time="12:00:00"
    log:
        "logs/scan_raw_indel_signals.{sample}.{dcs_chunk}.log"
    run:
        import pysam, pandas as pd, pickle, contextlib
        from intervaltree import IntervalTree
        from collections import defaultdict

        # 1. Load BED Regions
        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):
            bed = pd.read_csv(input.bed_chunk, sep="\t", header=None, usecols=[0,1,2],
                            names=["chrom","start","end"], dtype={"chrom": str})
            trees = {}
            for chrom in bed['chrom'].unique():
                tree = IntervalTree()
                for _, row in bed[bed['chrom'] == chrom].iterrows():
                    tree.addi(int(row.start), int(row.end)+1, True)
                trees[chrom] = tree

            # 2. Open Files
            samfile = pysam.AlignmentFile(input.dcs_bam, "rb")
            fasta = pysam.FastaFile(input.ref)

            # 3. Initialize Maps
            # Use sets for fast unique insertion
            candidate_pos = defaultdict(set)
            snv_noise_pos = defaultdict(set)

            # 4. The Scan Loop
            for chrom in trees.keys():
                # Prefetch reference sequence for the whole chromosome (or chunk) to speed up SNV checks
                # Note: If chromosomes are huge, fetch per-interval, but for chunks this is fine.
                # We will just fetch per-site to be memory safe.
                
                for interval in trees[chrom]:
                    # truncate=True prevents crashing at ends
                    for col in samfile.pileup(chrom, int(interval.begin), int(interval.end), truncate=True):
                        
                        pos = col.reference_pos + 1 # 1-based coordinate
                        
                        # Fetch Ref Base for SNV check (Optimization: cache this if slow)
                        try:
                            ref_base = fasta.fetch(chrom, pos-1, pos).upper()
                        except (KeyError, ValueError):
                            continue # Skip if outside ref bounds

                        has_indel = False
                        has_mismatch = False

                        for pr in col.pileups:
                            # Skip skips (N) and unmapped parts
                            if pr.is_refskip or pr.query_position is None: continue

                            # A. Check for Indel (Target)
                            # pr.indel != 0 means insertion/deletion starts *after* this base
                            if pr.indel != 0:
                                has_indel = True

                            # B. Check for SNV (Noise)
                            # Only check if it's NOT a deletion at this specific base
                            if not pr.is_del:
                                read_base = pr.alignment.query_sequence[pr.query_position].upper()
                                if read_base != ref_base:
                                    has_mismatch = True
                            
                            # Optimization: If we found both, we can break early for this site
                            if has_indel and has_mismatch:
                                break
                        
                        if has_indel:
                            candidate_pos[chrom].add(pos)
                        if has_mismatch:
                            snv_noise_pos[chrom].add(pos)

            # 5. Convert to Sorted Lists (for Binary Search in next step)
            final_data = {
                "candidates": {k: sorted(list(v)) for k, v in candidate_pos.items()},
                "noise_snv": {k: sorted(list(v)) for k, v in snv_noise_pos.items()},
                "noise_indel": {k: sorted(list(v)) for k, v in candidate_pos.items()} 
            }

            with open(output.scan_data, "wb") as f:
                pickle.dump(final_data, f)
            print(f"Scanned {input.dcs_bam} for {input.bed_chunk}: Found {sum(len(v) for v in candidate_pos.values())} indel targets and {sum(len(v) for v in snv_noise_pos.values())} SNV noise positions.")

rule make_raw_indel_vcf:
    """
    (SCATTER) Generates a raw VCF for Indels with exhaustive feature extraction.
    Replaces Mutect2 to provide a 'rawest possible' pileup view.
    """
    input:
        scan_data = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_signals/{{dcs_chunk}}.pkl",
        dcs_bam   = lambda wc: config["samples"][wc.sample]["dcs_bam"],
        ref       = REF,
        depth_lookup = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/depth_percentile_lookup.pkl",
        # Optional: Add specific blacklists if you have them, otherwise use a generic one
        blacklist = config["PATHS"]["encode_blacklist"]
    output:
        vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_indels.{{dcs_chunk}}.vcf")
    threads: 2
    resources:
        mem_mb=16000,
        time="12:00:00"
    log:
        "logs/make_raw_indel_vcf.{sample}.{dcs_chunk}.log"
    run:
        import pysam, pickle, statistics, math, bisect, pandas as pd
        from collections import defaultdict, Counter
        from intervaltree import IntervalTree
        import contextlib

        # --- 1. HELPER FUNCTIONS ---

        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):

            def get_entropy(seq):
                """Calculates Shannon entropy of a sequence string."""
                if not seq: return 0.0
                counts = Counter(seq)
                total = len(seq)
                return -sum((c/total) * math.log2(c/total) for c in counts.values())

            def get_gc(seq):
                if not seq: return 0.0
                return (seq.count('G') + seq.count('C')) / len(seq)

            def get_homopolymer_len(ref_obj, chrom, pos):
                """Returns length of the homopolymer run at this Reference position."""
                try:
                    base = ref_obj.fetch(chrom, pos-1, pos).upper()
                    run = 0
                    # Look ahead up to 20bp
                    for i in range(20):
                        # fetch pos to pos+1 (1-based logic in pysam fetch is start:end)
                        if ref_obj.fetch(chrom, pos+i, pos+i+1).upper() == base: 
                            run += 1
                        else: 
                            break
                    return run
                except: return 0

            def safe_mean(l):
                return statistics.mean(l) if l else 0.0

            def safe_sd(l):
                return statistics.stdev(l) if len(l) > 1 else 0.0

            # --- 2. LOAD DATA ---

            # Load Blacklist into IntervalTree
            blacklist_trees = defaultdict(IntervalTree)
            try:
                # Assumes 3-column BED [chrom, start, end]
                df = pd.read_csv(input.blacklist, sep="\t", header=None, comment='#', 
                                usecols=[0,1,2], names=["chrom", "start", "end"], dtype={"chrom": str})
                for _, row in df.iterrows():
                    blacklist_trees[row.chrom].addi(row.start, row.end)
            except Exception as e:
                print(f"Warning: Could not load blacklist: {e}")

            # Load Pre-Scanned Maps
            with open(input.scan_data, "rb") as f:
                data = pickle.load(f)
                candidates = data["candidates"] # {chrom: [sorted_pos...]}
                noise_all = data["noise_snv"]   # {chrom: [sorted_pos...]} (Using SNV noise as proxy for context noise)
                
            with open(input.depth_lookup, "rb") as f:
                depth_percentile_map = pickle.load(f)

            samfile = pysam.AlignmentFile(input.dcs_bam, "rb")
            ref_fasta = pysam.FastaFile(input.ref)

            # --- 3. MAIN PROCESSING LOOP ---

            with open(output.vcf, "w") as fout:
                # 1. Standard Header
                fout.write("##fileformat=VCFv4.2\n")
                fout.write(f"##source=IndelScanner_v1.0\n")
                
                # 2. DEFINITIONS FOR ALL YOUR CUSTOM METRICS (Critical Fix)
                # Basic Counts
                fout.write('##INFO=<ID=INDEL_LENGTH,Number=1,Type=Integer,Description="Length of indel (positive=ins, negative=del)">\n')
                fout.write('##INFO=<ID=N_TOTAL,Number=1,Type=Integer,Description="Total reads at locus">\n')
                fout.write('##INFO=<ID=N_ALT,Number=1,Type=Integer,Description="Number of reads supporting indel">\n')
                fout.write('##INFO=<ID=DEPTH_PERCENTILE,Number=1,Type=Float,Description="Depth percentile relative to region">\n')
                
                # Quality & Mapping
                fout.write('##INFO=<ID=MEAN_INSET_QUAL,Number=1,Type=Float,Description="Mean quality of inserted bases">\n')
                fout.write('##INFO=<ID=MEAN_DEL_ANCHOR_QUAL,Number=1,Type=Float,Description="Mean quality of anchor base for deletions">\n')
                fout.write('##INFO=<ID=SD_BASEQ,Number=1,Type=Float,Description="Standard deviation of base quality">\n')
                fout.write('##INFO=<ID=AVG_ALT_READ_MEAN_BQ,Number=1,Type=Float,Description="Mean BQ of alt reads">\n')
                fout.write('##INFO=<ID=MEAN_ALT_MAPQ,Number=1,Type=Float,Description="Mean mapping quality of alt reads">\n')
                fout.write('##INFO=<ID=MEAN_REF_MAPQ,Number=1,Type=Float,Description="Mean mapping quality of ref reads">\n')
                
                # Alignment Characteristics
                fout.write('##INFO=<ID=AVG_ALT_ASXS,Number=1,Type=Float,Description="Average AS/XS score ratio for alt reads">\n')
                fout.write('##INFO=<ID=AVG_REF_ASXS,Number=1,Type=Float,Description="Average AS/XS score ratio for ref reads">\n')
                fout.write('##INFO=<ID=AVG_ALT_EXCESS_NM,Number=1,Type=Float,Description="Excess mismatches (NM) in alt reads">\n')
                fout.write('##INFO=<ID=AVG_ALT_N_COUNT,Number=1,Type=Float,Description="Average N bases in alt reads">\n')
                fout.write('##INFO=<ID=AVG_ALT_SOFTCLIP,Number=1,Type=Float,Description="Average softclipping length in alt reads">\n')
                fout.write('##INFO=<ID=MEDIAN_DIST_TO_SOFTCLIP,Number=1,Type=Float,Description="Median distance to nearest softclip">\n')
                fout.write('##INFO=<ID=AVG_ALTREAD_VARIANTS,Number=1,Type=Float,Description="Avg other variants on alt reads">\n')
                fout.write('##INFO=<ID=AVG_NONALTREAD_VARIANTS,Number=1,Type=Float,Description="Avg other variants on ref reads">\n')
                
                # Insert Size & Positioning
                fout.write('##INFO=<ID=MEAN_INSERT,Number=1,Type=Float,Description="Mean insert size">\n')
                fout.write('##INFO=<ID=SD_INSERT,Number=1,Type=Float,Description="Standard deviation of insert size">\n')
                fout.write('##INFO=<ID=AVG_READ_POSITION,Number=1,Type=Float,Description="Average relative read position (0-1)">\n')
                fout.write('##INFO=<ID=AVG_READ_LENGTH,Number=1,Type=Float,Description="Average read length">\n')
                
                # Context / Sequence
                fout.write('##INFO=<ID=REF_ENTROPY_40BP,Number=1,Type=Float,Description="Shannon entropy of reference context">\n')
                fout.write('##INFO=<ID=VARIANTS_20BP,Number=1,Type=Integer,Description="Number of other variants within 20bp">\n')
                fout.write('##INFO=<ID=VARIANTS_250BP,Number=1,Type=Integer,Description="Number of other variants within 250bp">\n')
                fout.write('##INFO=<ID=HOMOPOLYMER_LEN,Number=1,Type=Integer,Description="Length of adjacent homopolymer">\n')
                fout.write('##INFO=<ID=GC_CONTENT,Number=1,Type=Float,Description="GC content of local region">\n')
                fout.write('##INFO=<ID=DIST_TO_OTHER_INDEL,Number=1,Type=Integer,Description="Distance to nearest other indel candidate">\n')
                fout.write('##INFO=<ID=IN_BLACKLIST,Number=1,Type=String,Description="Whether variant overlaps a blacklist region (YES/NO)">\n')

                # 3. Column Header
                fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                for chrom, pos_list in candidates.items():
                    for pos in pos_list:
                        
                        # A. PILEUP ANALYSIS
                        # ------------------
                        # Fetch pileup at exactly this position
                        # pysam is 0-based, so pos-1. 
                        cols = samfile.pileup(chrom, pos-1, pos, truncate=True)
                        target_col = None
                        for col in cols:
                            if col.reference_pos == pos-1:
                                target_col = col
                                break
                        
                        if not target_col: continue

                        # Identify the Major Indel
                        indel_counter = Counter()
                        for pr in target_col.pileups:
                            if pr.is_refskip or pr.query_position is None: continue
                            if pr.indel != 0:
                                indel_counter[pr.indel] += 1
                        
                        if not indel_counter: continue
                        best_indel_len, n_alt = indel_counter.most_common(1)[0]

                        # Separate Reads
                        alt_reads = []
                        ref_reads = []
                        for pr in target_col.pileups:
                            if pr.is_refskip or pr.query_position is None: continue
                            if pr.indel == best_indel_len:
                                alt_reads.append(pr)
                            elif pr.indel == 0:
                                ref_reads.append(pr)

                        n_total = len(alt_reads) + len(ref_reads)
                        if n_total == 0: continue

                        # B. EXTRACT FEATURES
                        # -------------------
                        
                        # Storage for raw values
                        f_alt_bqs = []       # Indel-specific BQ
                        f_alt_read_bqs = []  # Whole read BQ
                        f_alt_mapq = []
                        f_alt_as = []
                        f_alt_excess_nm = [] # NM minus indel length
                        f_alt_n_count = []
                        f_alt_softclip = []
                        f_alt_dist_clip = []
                        f_alt_tlen = []
                        f_alt_pos = []       # Relative position (0-1)
                        f_alt_len = []

                        f_ref_mapq = []
                        f_ref_as = []
                        f_ref_excess_nm = [] # Pure NM (noise)

                        ref_base = ref_fasta.fetch(chrom, pos-1, pos).upper()
                        alt_base = "" 

                        # 1. Process ALT Reads
                        for pr in alt_reads:
                            aln = pr.alignment
                            qp = pr.query_position
                            
                            # -- Indel Specific BaseQ --
                            if best_indel_len > 0: # Insertion
                                # Qualities of inserted bases
                                inserted_quals = aln.query_qualities[qp+1 : qp+1+best_indel_len]
                                val = sum(inserted_quals)/len(inserted_quals) if len(inserted_quals) > 0 else 0
                                f_alt_bqs.append(val)
                                
                                # Construct VCF String
                                if not alt_base:
                                    ins_seq = aln.query_sequence[qp+1 : qp+1+best_indel_len]
                                    alt_base = ref_base + ins_seq
                            else: # Deletion
                                # Quality of anchor base
                                f_alt_bqs.append(aln.query_qualities[qp])
                                
                                # Construct VCF String
                                if not alt_base:
                                    del_seq = ref_fasta.fetch(chrom, pos, pos + abs(best_indel_len)).upper()
                                    ref_base = ref_base + del_seq
                                    alt_base = ref_base[0]

                            # -- General Stats --
                            f_alt_read_bqs.append(safe_mean(aln.query_qualities))
                            f_alt_mapq.append(aln.mapping_quality)
                            f_alt_as.append(aln.get_tag("AS") if aln.has_tag("AS") else 0)
                            f_alt_n_count.append(aln.query_sequence.count('N'))
                            f_alt_tlen.append(abs(aln.template_length))
                            f_alt_len.append(aln.query_length)
                            f_alt_pos.append(qp / aln.query_length if aln.query_length > 0 else 0.5)

                            # -- Excess NM --
                            nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
                            f_alt_excess_nm.append(max(0, nm - abs(best_indel_len)))

                            # -- Softclip --
                            cigar = aln.cigartuples or []
                            left = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
                            right = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
                            f_alt_softclip.append(left + right)
                            
                            # Distance to clip
                            dist_l = qp
                            dist_r = aln.query_length - qp
                            f_alt_dist_clip.append(min(dist_l, dist_r))

                        # 2. Process REF Reads
                        for pr in ref_reads:
                            aln = pr.alignment
                            f_ref_mapq.append(aln.mapping_quality)
                            f_ref_as.append(aln.get_tag("AS") if aln.has_tag("AS") else 0)
                            nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
                            f_ref_excess_nm.append(nm)

                        # C. AGGREGATE STATS
                        # ------------------
                        mean_inset_qual = safe_mean(f_alt_bqs) if best_indel_len > 0 else 0
                        mean_del_anchor_qual = safe_mean(f_alt_bqs) if best_indel_len < 0 else 0
                        sd_baseq = safe_sd(f_alt_bqs)
                        
                        avg_alt_read_mean_bq = safe_mean(f_alt_read_bqs)
                        mean_alt_mapq = safe_mean(f_alt_mapq)
                        avg_alt_as = safe_mean(f_alt_as)
                        mean_ref_mapq = safe_mean(f_ref_mapq)
                        avg_ref_as = safe_mean(f_ref_as)
                        avg_excess_nm = safe_mean(f_alt_excess_nm) # Avg Altread Variants
                        avg_ref_vars = safe_mean(f_ref_excess_nm)  # Avg Non-Altread Variants
                        avg_alt_n = safe_mean(f_alt_n_count)
                        avg_softclip = safe_mean(f_alt_softclip)
                        median_dist_clip = statistics.median(f_alt_dist_clip) if f_alt_dist_clip else 0
                        mean_insert = safe_mean(f_alt_tlen)
                        sd_insert = safe_sd(f_alt_tlen)
                        avg_read_pos = safe_mean(f_alt_pos)
                        avg_read_len = safe_mean(f_alt_len)

                        # D. CONTEXT FEATURES
                        # -------------------
                        percentile = depth_percentile_map.get(n_total, 50)
                        
                        # Noise in 20bp/250bp (using pre-calculated SNV noise map)
                        vars_20 = 0
                        vars_250 = 0
                        if chrom in noise_all:
                            c_list = noise_all[chrom]
                            # 20bp
                            l20 = bisect.bisect_left(c_list, pos - 20)
                            r20 = bisect.bisect_right(c_list, pos + 20)
                            vars_20 = max(0, r20 - l20 - 1)
                            # 250bp
                            l250 = bisect.bisect_left(c_list, pos - 250)
                            r250 = bisect.bisect_right(c_list, pos + 250)
                            vars_250 = max(0, r250 - l250 - 1)

                        # Dist to other Indel
                        dist_to_indel = 9999
                        if chrom in candidates:
                            idx = bisect.bisect_left(candidates[chrom], pos)
                            if idx > 0:
                                dist_to_indel = min(dist_to_indel, abs(pos - candidates[chrom][idx-1]))
                            if idx + 1 < len(candidates[chrom]):
                                dist_to_indel = min(dist_to_indel, abs(pos - candidates[chrom][idx+1]))

                        # Ref Context
                        ref_ctx_40 = ref_fasta.fetch(chrom, max(0, pos-20), pos+20).upper()
                        ref_entropy = get_entropy(ref_ctx_40)
                        gc_content = get_gc(ref_ctx_40)
                        homopolymer = get_homopolymer_len(ref_fasta, chrom, pos)

                        # E. FLAGS & BLACKLIST
                        # --------------------
                        in_blacklist = "NO"
                        if chrom in blacklist_trees:
                            if blacklist_trees[chrom].overlaps(pos, pos+1):
                                in_blacklist = "YES"

                        # Custom Flag Logic (Expandable)
                        flags_str = []
                        # Example placeholders if you add more trees:
                        # if in_encode_tree: flags_str.append("IN_ENCODE_BLACKLIST")
                        # if in_repeat_tree: flags_str.append("IN_REPEATMASKER")
                        
                        flags_out = ";" + ";".join(flags_str) if flags_str else ""

                        # F. CONSTRUCT VCF LINE
                        # ---------------------
                        info = (
                            f"INDEL_LENGTH={best_indel_len};"
                            f"N_TOTAL={n_total};N_ALT={n_alt};DEPTH_PERCENTILE={percentile};"
                            f"IN_BLACKLIST={in_blacklist};"
                            f"MEAN_INSET_QUAL={mean_inset_qual:.2f};MEAN_DEL_ANCHOR_QUAL={mean_del_anchor_qual:.2f};"
                            f"SD_BASEQ={sd_baseq:.2f};"
                            f"AVG_ALT_READ_MEAN_BQ={avg_alt_read_mean_bq:.2f};"
                            f"MEAN_ALT_MAPQ={mean_alt_mapq:.2f};AVG_ALT_ASXS={avg_alt_as:.2f};"
                            f"MEAN_REF_MAPQ={mean_ref_mapq:.2f};AVG_REF_ASXS={avg_ref_as:.2f};"
                            f"AVG_ALT_EXCESS_NM={avg_excess_nm:.2f};AVG_ALT_N_COUNT={avg_alt_n:.2f};"
                            f"AVG_ALT_SOFTCLIP={avg_softclip:.2f};MEDIAN_DIST_TO_SOFTCLIP={median_dist_clip:.1f};"
                            f"AVG_ALTREAD_VARIANTS={avg_excess_nm:.2f};AVG_NONALTREAD_VARIANTS={avg_ref_vars:.2f};"
                            f"MEAN_INSERT={mean_insert:.1f};SD_INSERT={sd_insert:.1f};"
                            f"AVG_READ_POSITION={avg_read_pos:.2f};AVG_READ_LENGTH={avg_read_len:.1f};"
                            f"REF_ENTROPY_40BP={ref_entropy:.3f};VARIANTS_20BP={vars_20};VARIANTS_250BP={vars_250};"
                            f"HOMOPOLYMER_LEN={homopolymer};GC_CONTENT={gc_content:.2f};DIST_TO_OTHER_INDEL={dist_to_indel}"
                            f"{flags_out}"
                        )

                        fout.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\t.\t{info}\n")
                        
            print("Feature extraction complete!")

#### Step 3.2: GATHER DCS ###
### The dcs info has been called in parallel, now we need to gather the scattered VCFs into a single VCF per sample ###

rule clean_raw_indel_scatter:
    input:
        # The output from the Python rule we just wrote
        raw_vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_indels.{{dcs_chunk}}.vcf",
        contigs_hdr = CONTIGS_HDR, # Essential: defines chromosome lengths for bcftools
        ref = REF
    output:
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/clean_scatter.indel.{{dcs_chunk}}.vcf.gz"),
        tbi = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/clean_scatter.indel.{{dcs_chunk}}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/clean_raw_indel_scatter.{sample}.{dcs_chunk}.log"
    shell:
        r"""
        (
        set -euo pipefail

        # 1. Annotate Header: 
        # The Python script creates a minimal header. We inject the full ##contig list here
        # so bcftools norm knows the chromosome lengths.
        bcftools annotate -h {input.contigs_hdr} {input.raw_vcf} -O v | \

        # 2. Normalize (CRITICAL FOR INDELS):
        # Left-aligns indels against the reference. 
        # -m - : Splits multiallelics (though our script produces 1 per line anyway).
        # -f : Reference fasta is required for realignment.
        bcftools norm -m - -f {input.ref} -O v | \

        # 3. Sort and Compress:
        bcftools sort -Oz -o {output.vcf_gz}

        # 4. Index
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """

rule gather_raw_indels:
    input:
        # Collects the cleaned chunks
        scattered_vcf = lambda wc: expand(f"{PROCESSED_DATA_DIR}/variants/{wc.sample}/clean_scatter.indel.{{dcs_chunk}}.vcf.gz", dcs_chunk=DCS_CHUNK_IDS),
        scattered_tbi = lambda wc: expand(f"{PROCESSED_DATA_DIR}/variants/{wc.sample}/clean_scatter.indel.{{dcs_chunk}}.vcf.gz.tbi", dcs_chunk=DCS_CHUNK_IDS)
    output:
        # The merged file containing ALL raw candidates for the sample
        gathered_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_indels.gathered.vcf")
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    log:
        "logs/gather_raw_indels.{sample}.log"
    shell:
        """
        (
        # -a: Reset header attributes (safe for gathered files)
        # -D: Remove duplicate entries if chunks slightly overlapped (safety net)
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        ) &> {log}
        """

# The final list of DCS candidates will guide the SSCS extraction.

rule indel_finalize_dcs_candidates:
    input:
        gathered_vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/raw_indels.gathered.vcf",
        # We don't strictly need contigs_hdr/ref here as they were applied in step 1,
        # but keeping ref for a sanity check sort is fine.
    output:
        # The final input for your Machine Learning model
        vcf_gz = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_raw_indel_candidates.vcf.gz",
        idx    = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_raw_indel_candidates.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/indel_finalize_dcs_candidates.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail
        
        # We just Sort and Compress.
        # The data was already normalized and cleaned in the scatter step.
        bcftools sort {input.gathered_vcf} -Oz -o {output.vcf_gz}
        
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """


### STEP 3.3: COLLECTING SSCS DATA USING THE DCS CANDIDATES GUIDE VCF ###


rule split_candidate_indel_vcf:
    input:
        dcs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_raw_indel_candidates.vcf.gz"
    output:
        chunks=temp(expand(
            "{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}/indel.chunk_{i}.vcf.gz",
            PROCESSED_DATA_DIR=PROCESSED_DATA_DIR, i=range(1, SSCS_CHUNKS + 1)
        ))
    params:
        n = SSCS_CHUNKS,
        outdir = f"{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}"
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/split_candidate_indel_vcf.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail

        mkdir -p {params.outdir}

        # Extract header once
        bcftools view -h {input.dcs_vcf} > {params.outdir}/header.vcf

        # Count variant lines
        DATA_LINES=$(bcftools view -H {input.dcs_vcf} | wc -l)

        if [ "$DATA_LINES" -eq 0 ]; then
            LINES_PER_CHUNK=1
        else
            LINES_PER_CHUNK=$(( (DATA_LINES + {params.n} - 1) / {params.n} ))
        fi

        # Extract all variant lines
        bcftools view -H {input.dcs_vcf} > {params.outdir}/all_variants.tmp

        # Split into exactly N chunks
        for i in $(seq 1 {params.n}); do
            start=$(( (i - 1) * LINES_PER_CHUNK + 1 ))
            end=$(( i * LINES_PER_CHUNK ))

            outfile="{params.outdir}/indel.chunk_${{i}}.vcf.gz"

            # Slice lines (even if empty)
            sed -n "${{start}},${{end}}p" {params.outdir}/all_variants.tmp > {params.outdir}/chunk_${{i}}.vcf

            # Combine header + data, compress, index
            cat {params.outdir}/header.vcf {params.outdir}/chunk_${{i}}.vcf | \
                bcftools sort -Oz -o $outfile

            tabix -p vcf $outfile

            rm {params.outdir}/chunk_${{i}}.vcf
        done

        # Cleanup
        rm {params.outdir}/header.vcf
        rm {params.outdir}/all_variants.tmp

        ) &> {log}
        """

rule indel_format_sscs_counts_scatter:
    input:
        candidate_vcf_chunk = f"{PROCESSED_DATA_DIR}/sscs_chunks/{{sample}}/indel.{{sscs_chunk}}.vcf.gz",
        sscs_bam = lambda wc: config["samples"][wc.sample]["sscs_bam"],
        sscs_bai = lambda wc: config["samples"][wc.sample]["sscs_bam"] + ".bai",
    output:
        sscs_chunk_vcf = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/indel.metrics.{{sscs_chunk}}.vcf")
    threads: 1
    resources:
        mem_mb=2000,
        time="05:00:00"
    log:
        "logs/indel_format_sscs_counts_scatter.{sample}.{sscs_chunk}.log"
    run:
        import pysam, statistics
        from collections import defaultdict
        import contextlib

        with open(log[0], "w") as log_f, contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):

            print("Starting SSCS metrics calculation for chunk:", input.candidate_vcf_chunk)

            def calculate_indel_metrics(pileupcolumn, target_indel_len):
                """
                pileupcolumn: pysam pileup at the ANCHOR base.
                target_indel_len: The length of the indel (+ for ins, - for del).
                """
                rel_positions = []
                alt_fwd, alt_rev = 0, 0
                
                n_exact_match = 0
                n_slippage = 0 # Reads with DIFFERENT indels
                n_ref = 0      # Reads with NO indel
                
                for pr in pileupcolumn.pileups:
                    # Skip unusable reads
                    if pr.is_del or pr.is_refskip or pr.query_position is None: 
                        continue
                    
                    # Check Indel Status
                    # pr.indel returns the length of the indel starting at the NEXT base.
                    # 0 = Ref, >0 = Ins, <0 = Del.
                    read_indel = pr.indel
                    
                    if read_indel == target_indel_len:
                        # EXACT MATCH
                        n_exact_match += 1
                        
                        # Store Metrics
                        rel_positions.append(pr.query_position / pr.alignment.query_length)
                        if pr.alignment.is_reverse:
                            alt_rev += 1
                        else:
                            alt_fwd += 1
                            
                    elif read_indel != 0:
                        # It has an indel, but the wrong length (Slippage/Noise)
                        n_slippage += 1
                    else:
                        # It matches the reference (No indel)
                        n_ref += 1

                n_total = n_exact_match + n_slippage + n_ref
                
                # Aggregate
                rp_med = statistics.median(rel_positions) if rel_positions else 0.5
                rp_sd = statistics.stdev(rel_positions) if len(rel_positions) > 1 else 0.0
                
                total_alt = alt_fwd + alt_rev
                sb = alt_fwd / total_alt if total_alt > 0 else 0.5

                return n_exact_match, n_slippage, n_ref, n_total, rp_med, rp_sd, sb

            # --- EXECUTION ---
            print("Opening BAM and VCF files...")
            sscs_bam = pysam.AlignmentFile(input.sscs_bam, "rb")
            vcf_in = pysam.VariantFile(input.candidate_vcf_chunk)
            
            with open(output.sscs_chunk_vcf, "w") as fout:
                # Header
                fout.write("##fileformat=VCFv4.2\n")
                fout.write("##INFO=<ID=N_TOTAL_SSCS,Number=1,Type=Integer,Description=\"Total usable SSCS reads\">\n")
                fout.write("##INFO=<ID=N_ALT_SSCS,Number=1,Type=Integer,Description=\"Exact indel match count (SSCS)\">\n")
                fout.write("##INFO=<ID=N_SLIPPAGE_SSCS,Number=1,Type=Integer,Description=\"Reads with different indels at this locus\">\n")
                fout.write("##INFO=<ID=RP_MED_SSCS,Number=1,Type=Float,Description=\"Median Read Position\">\n")
                fout.write("##INFO=<ID=RP_SD_SSCS,Number=1,Type=Float,Description=\"SD Read Position\">\n")
                fout.write("##INFO=<ID=STRAND_BIAS_SSCS,Number=1,Type=Float,Description=\"Strand Bias\">\n")
                fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                for record in vcf_in:
                    chrom = record.chrom
                    pos = record.pos # 1-based POS from VCF (Anchor Base)
                    ref = record.ref
                    alt = record.alts[0]
                    
                    # Calculate Target Indel Length
                    # Deletion: ref=AT, alt=A -> len(A)-len(AT) = 1-2 = -1
                    # Insertion: ref=A, alt=AT -> len(AT)-len(A) = 2-1 = +1
                    target_len = len(alt) - len(ref)
                    
                    # Fetch Pileup at Anchor Base (pos-1)
                    try:
                        # We check the exact anchor base. 
                        # If the read supports the indel, pr.indel at this base will be non-zero.
                        iter_pileup = sscs_bam.pileup(chrom, pos - 1, pos, truncate=True)
                        col = next(iter_pileup, None)
                    except: col = None
                    
                    if col:
                        n_alt, n_slip, n_ref, n_tot, rp_med, rp_sd, sb = calculate_indel_metrics(col, target_len)
                    else:
                        n_alt, n_slip, n_ref, n_tot, rp_med, rp_sd, sb = 0,0,0,0, 0.5, 0.0, 0.5
                    
                    # Info String
                    info = (
                        f"N_TOTAL_SSCS={n_tot};N_ALT_SSCS={n_alt};N_SLIPPAGE_SSCS={n_slip};"
                        f"RP_MED_SSCS={rp_med:.3f};RP_SD_SSCS={rp_sd:.3f};"
                        f"STRAND_BIAS_SSCS={sb:.3f}"
                    )
                    
                    fout.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

            print("SSCS metrics calculation complete for chunk:", input.candidate_vcf_chunk)

rule indel_clean_sscs_scatter:
    input:
        feature_vcf = f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/indel.metrics.{{sscs_chunk}}.vcf",
        contigs_hdr = CONTIGS_HDR,
        ref = REF
    output:
        vcf_gz = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/indel.clean.{{sscs_chunk}}.vcf.gz"),
        tbi = temp(f"{PROCESSED_DATA_DIR}/sscs_scatter/{{sample}}/indel.clean.{{sscs_chunk}}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=1000,
        time="10:00:00"
    log:
        "logs/indel_clean_sscs_scatter.{sample}.{sscs_chunk}.log"
    shell:
        r"""
        (
        bcftools annotate -h {input.contigs_hdr} {input.feature_vcf} | \
        bcftools norm -m - -f {input.ref} - | \
        bcftools sort -Oz -o {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """

rule indel_gather_sscs_counts:
    input:
        scattered_vcf = lambda wc: [f"{PROCESSED_DATA_DIR}/sscs_scatter/{wc.sample}/indel.clean.{chunk}.vcf.gz" for chunk in SSCS_CHUNK_IDS],
        idx = lambda wc: [f"{PROCESSED_DATA_DIR}/sscs_scatter/{wc.sample}/indel.clean.{chunk}.vcf.gz.tbi" for chunk in SSCS_CHUNK_IDS]
    output:
        gathered_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_indel_sscs_metrics.vcf.gz",
        idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_indel_sscs_metrics.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    log:
        "logs/indel_gather_sscs_counts.{sample}.log"
    shell:
        """
        (
        bcftools concat -a -D -O z -o {output.gathered_vcf} {input.scattered_vcf}
        tabix -p vcf {output.gathered_vcf}
        ) &> {log}
        """

### STEP 3.4: COMBINE USE THE SSCS DATA TO ANNOTATE THE DCS VCF, CREATING THE FINAL VARIANT CALLSET ###


rule indel_add_header:
    input:
        dcs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_raw_indel_candidates.vcf.gz",
        sscs_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_indel_sscs_metrics.vcf.gz",
        dcs_idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_raw_indel_candidates.vcf.gz.tbi",
        sscs_idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}_indel_sscs_metrics.vcf.gz.tbi"
    output:
        vcf_gz = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_with_header.vcf.gz",
        idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_with_header.vcf.gz.tbi"
    params:
        fields='INFO/N_ALT_SSCS,INFO/N_TOTAL_SSCS,INFO/N_SLIPPAGE_SSCS,INFO/STRAND_BIAS_SSCS,INFO/RP_MED_SSCS,INFO/RP_SD_SSCS'
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/indel_add_header.{sample}.log"
    shell:
        r"""
        (
        bcftools annotate \
          -a {input.sscs_vcf} \
          -c "{params.fields}" \
          -Oz -o {output.vcf_gz} \
          {input.dcs_vcf}
          
        tabix -p vcf {output.vcf_gz}
        ) &> {log}
        """

###########################################
# SNV annotation
###########################################
rule merge_pon_snv:
    input:
        snv_pon_csv = SNV_PON_CSV
    output:
        pon_list = temp("temp/pon_file_list.txt"),
        merged_vcf = temp(f"{PROCESSED_DATA_DIR}/variants/snv_merged_pon.vcf.gz"),
        merged_vcf_idx = temp(f"{PROCESSED_DATA_DIR}/variants/snv_merged_pon.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/merge_pon_snv.log"
    shell:
        """
        (
        grep -v '^$' {input.snv_pon_csv} | cut -d',' -f1 > {output.pon_list}

        bcftools merge \
            --file-list {output.pon_list} \
            --force-samples \
        | bcftools view -G -Oz -o {output.merged_vcf}

        bcftools index -t {output.merged_vcf}
        ) &> {log}
        """


rule snv_annotate_gnomad:
    input:
        vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined.vcf.gz",
        gnomad_vcf = GNOMAD_VCF
    output:
        vcf = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}.snv_gnomad_annotated.vcf.gz"),
        idx = temp(f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}.snv_gnomad_annotated.vcf.gz.tbi")
    params:
        header = "gnomad_header.txt"
    threads: 2
    resources:
        mem_mb=8000,
        time="00:30:00"
    shadow: "minimal"
    log:
        "logs/snv_annotate_gnomad.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail

        # ------------------------------------------------------------------
        # Create INFO header definition
        # ------------------------------------------------------------------

        echo '##INFO=<ID=IN_GNOMAD_COMMON,Number=1,Type=String,Description="Variant present in gnomAD common variants (Yes/No)">' > gnomad_header.txt

        # ------------------------------------------------------------------
        # Run intersection
        # ------------------------------------------------------------------

        bcftools isec -p isec {input.vcf} {input.gnomad_vcf}

        # ------------------------------------------------------------------
        # Variants NOT in gnomAD
        # ------------------------------------------------------------------

        bcftools view isec/0000.vcf \
        | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_GNOMAD_COMMON=No"; else $8=$8";IN_GNOMAD_COMMON=No"; print }}' \
        | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_no.vcf.gz

        bcftools index snv_tagged_no.vcf.gz

        # ------------------------------------------------------------------
        # Variants IN gnomAD
        # ------------------------------------------------------------------

        if [ -f isec/0002.vcf ]; then

            bcftools view isec/0002.vcf \
            | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_GNOMAD_COMMON=Yes"; else $8=$8";IN_GNOMAD_COMMON=Yes"; print }}' \
            | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_yes.vcf.gz

        else

            # create empty VCF with header if no shared variants
            bcftools view -h {input.vcf} \
            | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_yes.vcf.gz

        fi

        bcftools index snv_tagged_yes.vcf.gz

        # ------------------------------------------------------------------
        # Combine results
        # ------------------------------------------------------------------

        bcftools concat -a snv_tagged_no.vcf.gz snv_tagged_yes.vcf.gz \
        | bcftools sort -Oz -o {output.vcf}

        bcftools index -t {output.vcf}

        ) &> {log}
        """

rule snv_annotate_pon_and_finalize:
    input:
        vcf = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}.snv_gnomad_annotated.vcf.gz",
        idx = f"{PROCESSED_DATA_DIR}/variants/{{sample}}/{{sample}}.snv_gnomad_annotated.vcf.gz.tbi",
        pon = f"{PROCESSED_DATA_DIR}/variants/snv_merged_pon.vcf.gz",
        merged_PON_idx = f"{PROCESSED_DATA_DIR}/variants/snv_merged_pon.vcf.gz.tbi"
    output:
        vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined_annotated.vcf.gz",
        idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.snv_final_combined_annotated.vcf.gz.tbi"
    params:
        header = "pon_header.txt"
    threads: 2
    resources:
        mem_mb=8000,
        time="00:30:00"
    shadow: "minimal"
    log:
        "logs/snv_annotate_pon_and_finalize.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail

        # ------------------------------------------------------------
        # Create INFO header definition
        # ------------------------------------------------------------

        echo '##INFO=<ID=IN_CB_PON,Number=1,Type=String,Description="Variant present in Panel of Normals (Yes/No)">' > {params.header}

        # ------------------------------------------------------------
        # Run intersection against PoN
        # ------------------------------------------------------------

        bcftools isec -p isec {input.vcf} {input.pon}

        # ------------------------------------------------------------
        # Variants NOT in PoN
        # ------------------------------------------------------------

        bcftools view isec/0000.vcf \
        | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_CB_PON=No"; else $8=$8";IN_CB_PON=No"; print }}' \
        | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_no.vcf.gz

        bcftools index snv_tagged_no.vcf.gz

        # ------------------------------------------------------------
        # Variants IN PoN
        # ------------------------------------------------------------

        if [ -f isec/0002.vcf ]; then

            bcftools view isec/0002.vcf \
            | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_CB_PON=Yes"; else $8=$8";IN_CB_PON=Yes"; print }}' \
            | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_yes.vcf.gz

        else

            # create empty VCF with correct header if no shared variants
            bcftools view -h {input.vcf} \
            | bcftools annotate -h {params.header} --no-version -Oz -o snv_tagged_yes.vcf.gz

        fi

        bcftools index snv_tagged_yes.vcf.gz

        # ------------------------------------------------------------
        # Final combine and sort
        # ------------------------------------------------------------

        bcftools concat -a snv_tagged_no.vcf.gz snv_tagged_yes.vcf.gz \
        | bcftools sort -Oz -o {output.vcf}

        bcftools index -t {output.vcf}

        ) &> {log}
        """


###########################################
# Indel annotation
###########################################

rule prepare_gnomad_indels:
    input:
        gnomad_vcf = GNOMAD_VCF,
        ref = REF
    output:
        gnomad_norm = temp(f"{PROCESSED_DATA_DIR}/variants/gnomad_split_norm.indel.vcf.gz"),
        gnomad_norm_idx = temp(f"{PROCESSED_DATA_DIR}/variants/gnomad_split_norm.indel.vcf.gz.tbi")
    params:
        regions = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00",
        partition="icelake"
    log:
        "logs/prepare_gnomad_indels.log"
    shell:
        """
        (
        bcftools view --regions {params.regions} {input.gnomad_vcf} -Ou \
        | bcftools norm -m - -f {input.ref} -Oz -o {output.gnomad_norm}
        bcftools index -t {output.gnomad_norm}
        ) &> {log}
        """


rule prepare_pon_indels:
    input:
        indel_pon_csv = INDEL_PON_CSV,
        ref = REF
    output:
        pon_merged = temp(f"{PROCESSED_DATA_DIR}/variants/merged_pon_split.indel.vcf.gz"),
        pon_merged_idx = temp(f"{PROCESSED_DATA_DIR}/variants/merged_pon_split.indel.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    log:
        "logs/prepare_pon_indels.log"
    shell:
        r"""
        (
        set -euo pipefail
        grep -v '^$' {input.indel_pon_csv} | cut -d',' -f1 > pon_file_list.txt

        bcftools merge --file-list pon_file_list.txt --force-samples -Ou \
        | bcftools norm -m - -f {input.ref} \
        | bcftools view -G -Oz -o {output.pon_merged}

        bcftools index -t {output.pon_merged}
        ) &> {log}
        """


rule annotate_sample_indels:
    input:
        vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_with_header.vcf.gz",
        gnomad = f"{PROCESSED_DATA_DIR}/variants/gnomad_split_norm.indel.vcf.gz",
        gnomad_idx = f"{PROCESSED_DATA_DIR}/variants/gnomad_split_norm.indel.vcf.gz.tbi",
        pon = f"{PROCESSED_DATA_DIR}/variants/merged_pon_split.indel.vcf.gz",
        pon_merged_idx = f"{PROCESSED_DATA_DIR}/variants/merged_pon_split.indel.vcf.gz.tbi",
        ref = REF
    output:
        final_vcf = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_annotated.vcf.gz",
        final_idx = f"{PROCESSED_DATA_DIR}/final_results/{{sample}}/{{sample}}.indel_final_annotated.vcf.gz.tbi"
    params:
        header = "indel_headers.txt"
    threads: 2
    resources:
        mem_mb=8000,
        time="01:00:00"
    shadow: "minimal"
    log:
        "logs/annotate_sample_indels.{sample}.log"
    shell:
        r"""
        (
        set -euo pipefail

        # --------------------------------------
        # Create temporary INFO header
        # --------------------------------------
        echo '##INFO=<ID=IN_GNOMAD_COMMON,Number=1,Type=String,Description="Variant present in gnomAD common variants (Yes/No)">' > {params.header}
        echo '##INFO=<ID=IN_CB_PON,Number=1,Type=String,Description="Variant present in Panel of Normals (Yes/No)">' >> {params.header}

        # --------------------------------------
        # Normalize sample VCF
        # --------------------------------------
        SAFE_NAME=$(basename {input.vcf} | sed 's/\\./_/g')
        bcftools norm -m - -f {input.ref} {input.vcf} -Oz -o sample_norm.vcf.gz
        bcftools index -t sample_norm.vcf.gz

        # --------------------------------------
        # Annotate against gnomAD
        # --------------------------------------
        bcftools isec -p isec_gnomad sample_norm.vcf.gz {input.gnomad} -n~10

        bcftools view isec_gnomad/0000.vcf \
        | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_GNOMAD_COMMON=No"; else $8=$8";IN_GNOMAD_COMMON=No"; print }}' \
        | bcftools annotate -h {params.header} --no-version -Oz -o gnomad_tagged_no.vcf.gz

        bcftools index gnomad_tagged_no.vcf.gz

        if [ -f isec_gnomad/0002.vcf ]; then
            bcftools view isec_gnomad/0002.vcf \
            | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_GNOMAD_COMMON=Yes"; else $8=$8";IN_GNOMAD_COMMON=Yes"; print }}' \
            | bcftools annotate -h {params.header} --no-version -Oz -o gnomad_tagged_yes.vcf.gz
        else
            bcftools view -h sample_norm.vcf.gz \
            | bcftools annotate -h {params.header} --no-version -Oz -o gnomad_tagged_yes.vcf.gz
        fi
        bcftools index gnomad_tagged_yes.vcf.gz

        bcftools concat -a gnomad_tagged_no.vcf.gz gnomad_tagged_yes.vcf.gz \
        | bcftools sort -Oz -o step1.vcf.gz
        bcftools index step1.vcf.gz

        # --------------------------------------
        # Annotate against PoN
        # --------------------------------------
        bcftools isec -p isec_pon step1.vcf.gz {input.pon}

        bcftools view isec_pon/0000.vcf \
        | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_CB_PON=No"; else $8=$8";IN_CB_PON=No"; print }}' \
        | bcftools annotate -h {params.header} --no-version -Oz -o pon_tagged_no.vcf.gz

        bcftools index pon_tagged_no.vcf.gz

        if [ -f isec_pon/0002.vcf ]; then
            bcftools view isec_pon/0002.vcf \
            | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if($8==".") $8="IN_CB_PON=Yes"; else $8=$8";IN_CB_PON=Yes"; print }}' \
            | bcftools annotate -h {params.header} --no-version -Oz -o pon_tagged_yes.vcf.gz
        else
            bcftools view -h step1.vcf.gz \
            | bcftools annotate -h {params.header} --no-version -Oz -o pon_tagged_yes.vcf.gz
        fi
        bcftools index pon_tagged_yes.vcf.gz

        # --------------------------------------
        # Final concatenate + sort
        # --------------------------------------
        bcftools concat -a pon_tagged_no.vcf.gz pon_tagged_yes.vcf.gz \
        | bcftools sort -Oz -o {output.final_vcf}
        bcftools index -t {output.final_vcf}
        ) &> {log}
        """