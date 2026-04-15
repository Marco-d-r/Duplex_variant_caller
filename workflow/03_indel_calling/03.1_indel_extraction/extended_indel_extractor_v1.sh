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
TARGETS_COMBINED_VCF = expand("final_results/{sample}/{sample}.indel_candidates_annotated.vcf.gz", sample=SAMPLES)



### RULE ALL ###

rule all:
    input:
        TARGETS_COMBINED_VCF

#### PART 1: SCATTER ####
### the idea here is to use the WES .bed file to create small regions to parallelize the variant calling ###
### the key metrics rule can run in parallel with this step ###

rule split_bed_regions:
    input:
        BED
    output:
        # This expands to 'bed_chunks/chunk_1.bed' ... 'bed_chunks/chunk_N.bed'
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
        
        # 2. Calculate lines per chunk (rounding up) using {params.n}
        # FIX: Replaced hardcoded '60' with '{params.n}'
        lines_per_chunk=$(( (total_lines + {params.n} - 1) / {params.n} ))
        
        # 3. Clean existing chunks
        rm -f bed_chunks/chunk_*.bed
        rm -f bed_chunks/split_tmp_*
        
        # 4. Split the file
        split -a 4 -l $lines_per_chunk {input} bed_chunks/split_tmp_
        
        # 5. Rename and add .bed extension
        i=1
        for f in bed_chunks/split_tmp_*; do
            mv "$f" "bed_chunks/chunk_$i.bed"
            i=$((i + 1))
        done
        
        echo "Split {input} into $i chunks (target was {params.n})."
        """

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
        # -a ensures we still output '0' values for target bases that were missed (crucial for accurate stats)
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


rule scan_raw_signals:
    """
    (SCATTER) Scans the BAM chunk to build a 'Map of the Terrain'.
    Output: A pickle file containing:
      1. Indel Targets: Positions where we will attempt to call variants.
      2. Noise Maps: Sorted lists of messy positions (SNVs/Indels) for context metrics.
    """
    input:
        dcs_bam = lambda wildcards: config["samples"][wildcards.sample]["dcs_bam"],
        bed_chunk = "bed_chunks/{region}.bed",
        ref = REF
    output:
        scan_data = temp("variants/{sample}/raw_signals/{region}.pkl")
    threads: 1
    resources:
        mem_gb = 4
    run:
        import pysam, pandas as pd, pickle
        from intervaltree import IntervalTree
        from collections import defaultdict

        # 1. Load BED Regions
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

rule make_raw_indel_vcf:
    """
    (SCATTER) Generates a raw VCF for Indels with exhaustive feature extraction.
    Replaces Mutect2 to provide a 'rawest possible' pileup view.
    """
    input:
        scan_data = "variants/{sample}/raw_signals/{region}.pkl",
        dcs_bam   = lambda wc: config["samples"][wc.sample]["dcs_bam"],
        ref       = REF,
        depth_lookup = "variants/{sample}/depth_percentile_lookup.pkl",
        # Optional: Add specific blacklists if you have them, otherwise use a generic one
        blacklist = config["REF_FILES"]["blacklist_bed"]
    output:
        vcf = temp("variants/{sample}/raw_indels.{region}.vcf")
    threads: 1
    resources:
        mem_gb  = 4
    run:
        import pysam, pickle, statistics, math, bisect, pandas as pd
        from collections import defaultdict, Counter
        from intervaltree import IntervalTree

        # --- 1. HELPER FUNCTIONS ---

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

#### PART 3: GATHER DCS ###
### The dcs info has been called in parallel, now we need to gather the scattered VCFs into a single VCF per sample ###

rule clean_raw_indel_scatter:
    input:
        # The output from the Python rule we just wrote
        raw_vcf = "variants/{sample}/raw_indels.{region}.vcf",
        contigs_hdr = CONTIGS_HDR, # Essential: defines chromosome lengths for bcftools
        ref = REF
    output:
        vcf_gz = temp("variants/{sample}/clean_scatter.{region}.vcf.gz"),
        tbi = temp("variants/{sample}/clean_scatter.{region}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
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

        """

rule gather_raw_indels:
    input:
        # Collects the cleaned chunks
        scattered_vcf = expand("variants/{{sample}}/clean_scatter.{region}.vcf.gz", region=DCS_REGION_WILDCARDS),
        scattered_tbi = expand("variants/{{sample}}/clean_scatter.{region}.vcf.gz.tbi", region=DCS_REGION_WILDCARDS)
    output:
        # The merged file containing ALL raw candidates for the sample
        gathered_vcf = temp("variants/{sample}/raw_indels.gathered.vcf")
    threads: 1
    resources:
        mem_gb = 8
    shell:
        """
        # -a: Reset header attributes (safe for gathered files)
        # -D: Remove duplicate entries if chunks slightly overlapped (safety net)
        bcftools concat -a -D -o {output.gathered_vcf} {input.scattered_vcf}
        """

# The final list of DCS candidates will guide the SSCS extraction.

rule finalize_dcs_candidates:
    input:
        gathered_vcf = "variants/{sample}/raw_indels.gathered.vcf",
        # We don't strictly need contigs_hdr/ref here as they were applied in step 1,
        # but keeping ref for a sanity check sort is fine.
    output:
        # The final input for your Machine Learning model
        vcf_gz = "final_results/{sample}/{sample}_raw_indel_candidates.vcf.gz",
        idx    = "final_results/{sample}/{sample}_raw_indel_candidates.vcf.gz.tbi"
    threads: 1
    resources:
        mem_gb = 4
    shell:
        r"""
        set -euo pipefail
        
        # We just Sort and Compress.
        # The data was already normalized and cleaned in the scatter step.
        bcftools sort {input.gathered_vcf} -Oz -o {output.vcf_gz}
        
        tabix -p vcf {output.vcf_gz}
        """


### STEP 4: COLLECTING SSCS DATA USING THE DCS CANDIDATES GUIDE VCF ###


rule split_candidate_vcf:
    input:
        dcs_vcf = "final_results/{sample}/{sample}_raw_indel_candidates.vcf.gz"
    output:
        chunks = temp(expand("sscs_chunks/{{sample}}/chunk_{i}.vcf.gz", i=range(1, SSCS_CHUNKS + 1)))
    params:
        n = SSCS_CHUNKS,
        outdir = "sscs_chunks/{sample}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        
        # 1. Extract Header
        bcftools view -h {input.dcs_vcf} > {params.outdir}/header.vcf
        
        # 2. Count Lines
        DATA_LINES=$(bcftools view -H {input.dcs_vcf} | wc -l)
        if [ "$DATA_LINES" -eq 0 ]; then
            LINES_PER_CHUNK=1
        else
            # Double braces for awk
            LINES_PER_CHUNK=$(awk -v lines="$DATA_LINES" -v n="{params.n}" 'BEGIN {{ print int(lines/n) + 1 }}')
        fi

        # 3. Split Data
        bcftools view -H {input.dcs_vcf} | \
          awk -v lpc=$LINES_PER_CHUNK -v out_prefix="{params.outdir}/candidate_data_" '
            {{
              filenum = int((NR-1) / lpc) + 1;
              file = sprintf("%s%02d.vcf", out_prefix, filenum);
              print $0 >> file;
            }}'
        
        # 4. Rebuild Chunks
        for i in $(seq -w 1 {params.n}); do
            # BASH FIX: No braces needed for simple vars, or use double braces if needed
            clean_i=$(echo $i | sed 's/^0*//')
            
            # CRITICAL FIX HERE: Use ${{ }} for bash variables inside strings
            outfile="{params.outdir}/chunk_${{clean_i}}.vcf.gz"
            infile="{params.outdir}/candidate_data_${{i}}.vcf"
            
            if [ -f "$infile" ]; then
                cat {params.outdir}/header.vcf $infile | bgzip -c > $outfile
                rm $infile
            else
                cat {params.outdir}/header.vcf | bgzip -c > $outfile
            fi
            tabix -p vcf $outfile
        done
        
        rm {params.outdir}/header.vcf
        """

rule format_sscs_counts_scatter:
    input:
        candidate_vcf_chunk = "sscs_chunks/{sample}/chunk_{sscs_chunk}.vcf.gz",
        sscs_bam = lambda wc: config["samples"][wc.sample]["sscs_bam"],
        sscs_bai = lambda wc: config["samples"][wc.sample]["sscs_bam"] + ".bai",
    output:
        sscs_chunk_vcf = temp("sscs_scatter/{sample}/metrics.{sscs_chunk}.vcf")
    threads: 1
    resources:
        mem_gb = 2
    run:
        import pysam, statistics
        
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

rule clean_sscs_scatter:
    input:
        feature_vcf = "sscs_scatter/{sample}/metrics.{sscs_chunk}.vcf",
        contigs_hdr = CONTIGS_HDR,
        ref = REF
    output:
        vcf_gz = temp("sscs_scatter/{sample}/clean.{sscs_chunk}.vcf.gz"),
        tbi = temp("sscs_scatter/{sample}/clean.{sscs_chunk}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_gb = 3
    shell:
        r"""
        bcftools annotate -h {input.contigs_hdr} {input.feature_vcf} | \
        bcftools norm -m - -f {input.ref} - | \
        bcftools sort -Oz -o {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """

rule gather_sscs_counts:
    input:
        scattered_vcf = expand("sscs_scatter/{{sample}}/clean.{sscs_chunk}.vcf.gz", sscs_chunk=range(1, SSCS_CHUNKS+1)),
        idx = expand("sscs_scatter/{{sample}}/clean.{sscs_chunk}.vcf.gz.tbi", sscs_chunk=range(1, SSCS_CHUNKS+1))
    output:
        gathered_vcf = "final_results/{sample}/{sample}_sscs_metrics.vcf.gz",
        idx = "final_results/{sample}/{sample}_sscs_metrics.vcf.gz.tbi"
    shell:
        """
        bcftools concat -a -D -O z -o {output.gathered_vcf} {input.scattered_vcf}
        tabix -p vcf {output.gathered_vcf}
        """

### STEP 5: COMBINE USE THE SSCS DATA TO ANNOTATE THE DCS VCF, CREATING THE FINAL VARIANT CALLSET ###


rule final_annotate:
    input:
        dcs_vcf = "final_results/{sample}/{sample}_raw_indel_candidates.vcf.gz",
        sscs_vcf = "final_results/{sample}/{sample}_sscs_metrics.vcf.gz",
        dcs_idx = "final_results/{sample}/{sample}_raw_indel_candidates.vcf.gz.tbi",
        sscs_idx = "final_results/{sample}/{sample}_sscs_metrics.vcf.gz.tbi"
    output:
        vcf_gz = "final_results/{sample}/{sample}.indel_candidates_annotated.vcf.gz",
        idx = "final_results/{sample}/{sample}.indel_candidates_annotated.vcf.gz.tbi"
    params:
        # Include the slippage count!
        fields='INFO/N_ALT_SSCS,INFO/N_TOTAL_SSCS,INFO/N_SLIPPAGE_SSCS,INFO/STRAND_BIAS_SSCS,INFO/RP_MED_SSCS,INFO/RP_SD_SSCS'
    shell:
        r"""
        bcftools annotate \
          -a {input.sscs_vcf} \
          -c "{params.fields}" \
          -Oz -o {output.vcf_gz} \
          {input.dcs_vcf}
          
        tabix -p vcf {output.vcf_gz}
        """

