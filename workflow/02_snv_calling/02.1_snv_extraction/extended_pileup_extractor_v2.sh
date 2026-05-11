##############################################
#  Scratch / cache setup
##############################################
import os


### Define samples and paths to input files from config globally ###

# Define SAMPLES based on the keys in the config file
SAMPLES = list(config["samples"].keys())

# Replace hardcoded paths with config lookups
REF = config["REF_FILES"]["ref_fasta"]
CONTIGS_HDR = config["REF_FILES"]["contigs_hdr"]
BED = config["REF_FILES"]["bed_targets"]
N_CHUNKS = config["PARAMS"]["n_chunks_dcs"] # For DCS splitting
SSCS_CHUNKS = config["PARAMS"]["n_chunks_sscs"] # For SSCS splitting


# Define the list of regions globally based on the splitting rule
DCS_REGION_WILDCARDS = [f"chunk_{i}" for i in range(1, N_CHUNKS + 1)]
# Define the SSCS region wildcard based on the splitting rule
SSCS_REGION_WILDCARDS = [f"chunk_{i}" for i in range(1, SSCS_CHUNKS + 1)]
# All final permanent outputs for a sample are now consolidated in 'final_results/{sample}/'
TARGETS_COMBINED_VCF = expand("final_results/{sample}/{sample}.final_combined.vcf.gz", sample=SAMPLES)
TARGETS_CLASSIFIED_VCF = expand("final_results/{sample}/{sample}.classified.vcf.gz", sample=SAMPLES)
# The plots are also a good final target
TARGETS_METRICS = expand("final_results/{sample}/{sample}.key_metrics.csv", sample=SAMPLES)

### RULE ALL ###

rule all:
    input:
        TARGETS_COMBINED_VCF
       # TARGETS_CLASSIFIED_VCF
        #TARGETS_METRICS,
        #TARGETS_PON_VCF
        

#### PART 1: SCATTER ####
### the idea here is to use the WES .bed file to create small regions to parallelize the variant calling ###
### the key metrics rule can run in parallel with this step ###

rule split_bed_regions:
    input:
        BED
    output:
        # This expands to 'bed_chunks/chunk_1.bed' ... 'bed_chunks/chunk_500.bed'
        chunks=temp(expand("bed_chunks/chunk_{i}.bed", i=range(1, N_CHUNKS + 1)))
    params:
        n = N_CHUNKS,
        outdir = "bed_chunks"
    threads: 1
    resources:
        mem_gb = 2
    shell:
        """
        set -euo pipefail
        mkdir -p bed_chunks
        
        # 1. Count valid data lines (excluding headers)
        total_lines=$(grep -c -v "^#" {input} || true)
        
        # 2. Calculate lines per chunk (rounding up)
        lines_per_chunk=$(( (total_lines + 60 - 1) / 60 ))
        
        # 3. Clean existing chunks
        rm -f bed_chunks/chunk_*.bed
        rm -f bed_chunks/split_tmp_*
        
        # 4. Split the file (macOS Compatible)
        # We removed '-d' and '--additional-suffix' which were causing the crash.
        # This will create files like: bed_chunks/split_tmp_aaaa, bed_chunks/split_tmp_aaab, etc.
        split -a 4 -l $lines_per_chunk {input} bed_chunks/split_tmp_
        
        # 5. Rename and add .bed extension
        # We iterate through the temp files (alphabetical order ensures correct sequence)
        i=1
        for f in bed_chunks/split_tmp_*; do
            mv "$f" "bed_chunks/chunk_$i.bed"
            i=$((i + 1))
        done
        
        echo "Split {input} into $i chunks."
        """

########

rule calculate_depth_lookup:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        # ADDED: Define where your bed file is located
        bed = config["REF_FILES"]["bed_targets"] 
    output:
        depth_lookup = "variants/{sample}/depth_percentile_lookup.pkl"
    threads: 3
    resources:
        mem_gb = 6
    run:
        import subprocess, pickle, sys

        # UPDATED CMD: Added '-b {input.bed}'
        # -b limits the calculation to regions in your BED file
        # -a ensures we still output '0' values for target bases that were missed 
        cmd = f"samtools depth -a -b {input.bed} {input.bam} | awk '{{count[$3]++}} END {{for (d in count) print d, count[d]}}'"
        
        depth_counts = {}
        total_bases = 0
        
        # Stream the output of the command
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, text=True)
        
        for line in process.stdout:
            try:
                depth, count = map(int, line.strip().split())
                depth_counts[depth] = count
                total_bases += count
            except ValueError:
                continue
        
        process.wait()

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

#### PART 2: PROCESS THE SCATTERED DCS ###
### Using the scattered BED chunks to parallelize variant calling at the DCS level ###


rule mutect2_dcs_scatter:
    input:
        ref = REF,
        bam = lambda wc: config["samples"][wc.sample]["dcs_bam"],
        bed_chunk = "bed_chunks/{region}.bed"
    output:
        raw_vcf = temp("variants/{sample}/raw_mutect2.{region}.vcf.gz")
    threads: 2
    resources:
        mem_gb = 6,
        # Ensure this dir exists or is created; distinct per job is safer but fixed is ok if filenames differ
        tmpdir = "tmp/mutect2_dcs_scatter" 
    shell:
        r"""
        ulimit -n 8192
        mkdir -p {resources.tmpdir}

        # Define temp file path
        TEMP_VCF="{resources.tmpdir}/raw_mutect2.{wildcards.sample}.{wildcards.region}.tmp.vcf"

        # FIX 1: Added backslashes (\) for line continuation
        gatk --java-options "-Xmx5G -Djava.io.tmpdir={resources.tmpdir}" Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed_chunk} \
            --native-pair-hmm-threads {threads} \
            -O $TEMP_VCF

        # Compress to output
        bgzip -c $TEMP_VCF > {output.raw_vcf}
        bcftools index -t {output.raw_vcf}

        # FIX 2: Clean up the uncompressed intermediate file to save disk space
        rm $TEMP_VCF
        """

rule gather_mutect2_raw:
    input:
        raw_vcf=lambda wildcards: expand("variants/{sample}/raw_mutect2.{region}.vcf.gz", sample=wildcards.sample, region=DCS_REGION_WILDCARDS)
    output:
        gathered_raw_vcf = temp("variants/{sample}/raw_mutect2_gathered.vcf.gz")
    threads: 1  # bcftools concat is mostly single-threaded I/O bound
    resources:
        mem_gb = 8
    shell:
        """
        # FIX 3: Added '-O z' to ensure the output is actually gzipped
        # Added --threads to utilize the requested core for compression
        bcftools concat -a -D -O z --threads {threads} -o {output.gathered_raw_vcf} {input.raw_vcf}
        
        # Good practice to index the gathered result immediately
        bcftools index -t {output.gathered_raw_vcf}
        """

rule extract_indels:
    input:
        # Input VCF: The gathered raw Mutect2 VCF (SNVs + Indels)
        vcf = "variants/{sample}/raw_mutect2_gathered.vcf.gz" 
    output:
        # Temporary output: The Indel VCF needed for distance calculation
        indel_vcf = temp("variants/{sample}/{sample}_mutect2_indels_only.vcf.gz"),
        idx = temp("variants/{sample}/{sample}_mutect2_indels_only.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
        # Select only INDELs/MNPs/Complex, sort, compress, and index.
        # This VCF contains only the POSITIONS of indels for distance calculation.
        bcftools view -v indels,mnps,other {input.vcf} | \
        bcftools sort -Oz -o {output.indel_vcf}
        bcftools index -t {output.indel_vcf}
        """

rule precalculate_dcs_chunk_data:
    """
    (SCATTER) Pre-calculates global data structures for a single BED chunk.
    This replaces the monolithic pre-calculation rule.
    """
    input:
        dcs_bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        ref = REF,
        bed_chunk = "bed_chunks/{region}.bed" # Input is one small chunk
    output:
        # Temporary output, one per chunk
        global_data_chunk = temp("variants/{sample}/chunk_data/{region}.pkl")
    threads: 
        1
    resources:
        mem_gb = 2  
    run:
        import pysam, pandas as pd, pickle
        from intervaltree import IntervalTree
        from collections import defaultdict
        
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
        trees = build_interval_trees(input.bed_chunk) 
        
        samfile = pysam.AlignmentFile(input.dcs_bam, "rb")
        ref_fasta = pysam.FastaFile(input.ref)

        # 1. Get VAFs/Positions for this chunk
        # Note: We do NOT sort/unique them here. We do that once in the gather step.
        chrom_variant_positions, chrom_variant_vafs = get_variant_positions_and_vafs_bed(samfile, trees, ref_fasta)

        # 2. Build the read-to-variant map for this chunk
        read_to_variants = build_read_to_variant_positions(samfile, trees, ref_fasta)

        # 3. Save the results for this chunk
        chunk_data = {
            "chrom_variant_positions": chrom_variant_positions,
            "chrom_variant_vafs": chrom_variant_vafs,
            "read_to_variants": read_to_variants
        }
        with open(output.global_data_chunk, "wb") as f:
            pickle.dump(chunk_data, f)

rule gather_dcs_global_data:
    """
    (GATHER) Merges the pre-calculated chunk data into a single global data file.
    """
    input:
        # This will collect all 500 chunk.pkl files for this {sample}
        chunks = expand("variants/{{sample}}/chunk_data/{region}.pkl", region=DCS_REGION_WILDCARDS)
    output:
        # The final, non-temporary file the rest of the workflow expects
        global_data = "variants/{sample}/dcs_global_data.pkl"
    threads:
        1
    resources:
        mem_gb = 16 
    run:
        import pickle
        from collections import defaultdict

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
        
        with open(output.global_data, "wb") as f:
            pickle.dump(final_global_data, f)

rule make_vcf_scatter:
    input:
        global_data  = "variants/{sample}/dcs_global_data.pkl", 
        depth_lookup = "variants/{sample}/depth_percentile_lookup.pkl", 
        dcs_bam      = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        bed_chunk    = "bed_chunks/{region}.bed", 
        ref          = REF,
        encode       = config["REF_FILES"]["encode_blacklist"],
        repeatmasker = config["REF_FILES"]["repeatmasker_blacklist"],
        indel_vcf    = "variants/{sample}/{sample}_mutect2_indels_only.vcf.gz"
    output:
        feature_vcf = temp("variants/{sample}/feature_rich.{region}.vcf")
    threads:
        1
    resources:
        mem_gb = 5
    run:
        import pysam, pandas as pd, math, statistics, bisect, gzip, pickle
        from intervaltree import IntervalTree
        from collections import defaultdict, Counter

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
        with open(input.global_data, "rb") as f:
            global_data = pickle.load(f)
            chrom_variant_positions = global_data["chrom_variant_positions"]
            read_to_variants = global_data["read_to_variants"]
        
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
                for pileupcolumn in samfile.pileup(chrom, start, end):
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

#### PART 3: GATHER DCS ###
### The dcs info has been called in parallel, now we need to gather the scattered VCFs into a single VCF per sample ###

rule clean_dcs_scatter:
    input:
        feature_vcf = "variants/{sample}/feature_rich.{region}.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # Temporary output that has a clean header, ready for concatenation
        vcf_gz = temp("variants/{sample}/cleaned_scatter.{region}.vcf.gz"),
        tbi = temp("variants/{sample}/cleaned_scatter.{region}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
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
        """

rule gather_vcf:
    input:
        # Collects all feature VCFs for the current sample across all regions
        scattered_vcf = expand("variants/{{sample}}/cleaned_scatter.{region}.vcf.gz", region=DCS_REGION_WILDCARDS),
        scattered_tbi = expand("variants/{{sample}}/cleaned_scatter.{region}.vcf.gz.tbi", region=DCS_REGION_WILDCARDS)
    output:
        # The final, full-sized, feature-rich VCF for the sample
        gathered_vcf = temp("variants/{sample}/full_feature_rich.vcf")
    threads:
        1
    resources:
        mem_gb = 8
    shell:
        """
        # Ensure the VCFs are concatenated correctly, preserving header order.
        # -a ensures that the header only comes from the first file.
        # -D suppresses duplicate INFO/FORMAT fields, useful if headers were slightly different.
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        """

# The final list of DCS candidates will guide the SSCS extraction.


rule clean_dcs_candidates:
    input:
        # Input VCF: Contains all features gathered from parallel runs (no filtering applied yet)
        gathered_vcf = "variants/{sample}/full_feature_rich.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # The output is the final, indexed guide VCF for SSCS
        vcf_gz = temp("final_results/{sample}/{sample}_dcs_candidates.vcf.gz"),
        idx = temp("final_results/{sample}/{sample}_dcs_candidates.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
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
        """

### STEP 4: COLLECTING SSCS DATA USING THE DCS CANDIDATES GUIDE VCF ###


rule split_candidate_vcf:
    input:
        dcs_vcf = "final_results/{sample}/{sample}_dcs_candidates.vcf.gz"
    output:
        # This output definition is correct and does not need to change
        chunks=temp(expand("sscs_chunks/{{sample}}/chunk_{i}.vcf.gz", i=range(1, SSCS_CHUNKS + 1)))
    params:
        n = SSCS_CHUNKS,
        outdir = "sscs_chunks/{sample}",
        split_suffix = 2 
    shell:
        r"""
        set -euo pipefail
        
        # 1. Prepare directory and extract header
        mkdir -p {params.outdir}
        bcftools view -h {input.dcs_vcf} > {params.outdir}/header.vcf
        
        # 2. Extract variant lines, count them, and determine lines per chunk
        DATA_LINES=$(bcftools view -H {input.dcs_vcf} | wc -l)
        # Add a check for 0 lines to avoid division by zero
        if [ "$DATA_LINES" -eq 0 ]; then
            LINES_PER_CHUNK=1
        else
            LINES_PER_CHUNK=$(awk "BEGIN {{ print int($DATA_LINES / {params.n}) + 1 }}")
        fi

        # --- THIS IS THE FIX ---
        # 3. Split the variant data lines using AWK (macOS/BSD compatible)
        bcftools view -H {input.dcs_vcf} | \
          awk -v lpc=$LINES_PER_CHUNK -v out_prefix="{params.outdir}/candidate_data_" '
            {{
              # Calculate file number (1-based)
              filenum = int((NR-1) / lpc) + 1;
              # Create zero-padded filename (e.g., candidate_data_01.vcf)
              file = sprintf("%s%02d.vcf", out_prefix, filenum);
              # Print the line to that file
              print $0 >> file;
            }}'
        
        # 4. Reconstruct, bgzip, and index each chunk
        # This 'for' loop and 'sed' logic will now work perfectly
        # because the awk command creates the 'candidate_data_01.vcf' files it expects.
        for datafile in {params.outdir}/candidate_data_*.vcf; do
            chunk_suffix=$(basename $datafile .vcf | sed 's/candidate_data_//')
            chunk_suffix_unpadded=$(echo $chunk_suffix | sed 's/^0*//')
            
            # This line correctly uses doubled-braces to avoid a Snakemake NameError
            output_file="{params.outdir}/chunk_${{chunk_suffix_unpadded}}.vcf.gz"
            
            # Combine header and data, compress, and index
            cat {params.outdir}/header.vcf $datafile | \
              bgzip -c > $output_file
              
            tabix -p vcf $output_file
            rm $datafile
        done
        
        rm {params.outdir}/header.vcf
        """

rule format_sscs_counts_scatter:
    input:
        # The small chunk of DCS candidates (the positions we need SSCS data for)
        candidate_vcf_chunk = "sscs_chunks/{sample}/{sscs_chunk}.vcf.gz",
        # The SSCS BAM file
        sscs_bam = lambda wildcards: config["samples"][wildcards.sample]["sscs_bam"],
        sscs_bai = lambda wildcards: config["samples"][wildcards.sample]["sscs_bam"] + ".bai",
    output:
        # Output VCF chunk containing SSCS counts and new metrics (e.g., RP_MED_SSCS)
        sscs_chunk_vcf = temp("sscs_scatter/{sample}/metrics.{sscs_chunk}.vcf")
    threads: 
        1
    resources:
        mem_gb = 2
    run:
        import pysam, statistics
        from collections import defaultdict

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
        
        sscs_bam = pysam.AlignmentFile(input.sscs_bam, "rb")
        vcf_in = pysam.VariantFile(input.candidate_vcf_chunk)
        
        try:
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

rule clean_sscs_scatter:
    input:
        feature_vcf = "sscs_scatter/{sample}/metrics.{sscs_chunk}.vcf",
        contigs_hdr = CONTIGS_HDR, # The file containing the final ##contig lines
        ref = REF
    output:
        # Temporary output that has a clean header, ready for concatenation
        vcf_gz = temp("sscs_scatter/{sample}/cleaned_scatter.{sscs_chunk}.vcf.gz"),
        tbi = temp("sscs_scatter/{sample}/cleaned_scatter.{sscs_chunk}.vcf.gz.tbi")
    
    threads: 1
    resources:
        mem_gb = 1
    shell:
        r"""
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
        """


rule gather_sscs_counts:
    input:
        # Collects all SSCS metric VCF chunks for the current sample
        scattered_vcf = expand("sscs_scatter/{{sample}}/cleaned_scatter.{sscs_chunk}.vcf.gz", sscs_chunk=SSCS_REGION_WILDCARDS),
        idx = expand("sscs_scatter/{{sample}}/cleaned_scatter.{sscs_chunk}.vcf.gz.tbi", sscs_chunk=SSCS_REGION_WILDCARDS)
    output:
        # The final, full VCF with all SSCS metrics
        gathered_vcf = temp("sscs_scatter/{sample}/full_sscs_metrics.vcf.gz")
    threads:
        1
    resources:
        mem_gb = 8
    shell:
        """
        # Concatenate VCFs, ensuring headers are merged correctly.
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        """

rule clean_sscs_vcf:
    input:
        gathered_vcf = "sscs_scatter/{sample}/full_sscs_metrics.vcf.gz",
        contigs_hdr = CONTIGS_HDR,
    output:
        # Final clean, compressed, and indexed SSCS VCF
        vcf_gz = temp("final_results/{sample}/{sample}_sscs_metrics_clean.vcf.gz"),
        idx = temp("final_results/{sample}/{sample}_sscs_metrics_clean.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf_gz})
        # Annotate header (ensures correct contig lines) and sort
        bcftools annotate -h {input.contigs_hdr} {input.gathered_vcf} | \
        bcftools sort -Oz -o {output.vcf_gz}
        
        # Index the VCF
        tabix -p vcf {output.vcf_gz}
        """

### STEP 5: COMBINE USE THE SSCS DATA TO ANNOTATE THE DCS VCF, CREATING THE FINAL VARIANT CALLSET ###


rule final_annotate:
    input:
        dcs_vcf = "final_results/{sample}/{sample}_dcs_candidates.vcf.gz", 
        sscs_vcf = "final_results/{sample}/{sample}_sscs_metrics_clean.vcf.gz",
        dcs_vcf_idx = "final_results/{sample}/{sample}_dcs_candidates.vcf.gz.tbi",
        sscs_vcf_idx = "final_results/{sample}/{sample}_sscs_metrics_clean.vcf.gz.tbi"
    output:
      vcf_gz = "final_results/{sample}/{sample}.final_combined.vcf.gz",
        idx = "final_results/{sample}/{sample}.final_combined.vcf.gz.tbi"
    params:
        annotation_fields='INFO/N_ALT_SSCS,INFO/N_TOTAL_SSCS,INFO/STRAND_BIAS_SSCS,INFO/RP_MED_SSCS,INFO/RP_SD_SSCS'
    shell:
        r"""
        bcftools annotate \
          -a {input.sscs_vcf} \
          -c "{params.annotation_fields}" \
          -Oz -o {output.vcf_gz} \
          {input.dcs_vcf}
          
        tabix -p vcf {output.vcf_gz}
        """

