import pysam
import json
import os
import pybedtools

# Get inputs from Snakemake
bam_path = snakemake.input.bam
bed_path = snakemake.input.bed
chrom = snakemake.wildcards.chrom
ts = snakemake.config["thresholds"]

def evaluate_filters(read, query_pos, ts):
    """
    Conserves the EXACT logic of your check_read function.
    """
    # Softclipping: cigar_stats[0][4] is the count of S bases
    cigar_stats = read.get_cigar_stats()
    soft_clip_count = cigar_stats[0][4] if cigar_stats else 0

    results = {
        "MQ": read.mapping_quality >= ts["MIN_MQ"],
        "INSERT": abs(read.template_length) <= ts["MAX_INSERT"],
        "ASXS": False, # Default to fail, updated below
        "NM": True,    # Default to pass, updated below
        "NCOUNT": read.query_sequence.upper().count('N') <= ts["MAX_N_COUNT"],
        "SOFTCLIP": soft_clip_count <= ts["MAX_SOFTCLIP"],
        "READPOS": (query_pos >= ts["MIN_READ_POS"] and 
                    (read.query_length - query_pos) >= (150 - ts["MAX_READ_POS"]))
    }
    
    # NM check (matches your original try/except KeyError)
    try:
        if read.has_tag("NM"):
            if read.get_tag("NM") > ts["MAX_NM"]:
                results["NM"] = False
    except KeyError: pass

    # ASXS check (matches your original AS - XS logic)
    try:
        as_val = read.get_tag("AS")
        xs_val = read.get_tag("XS") if read.has_tag("XS") else 0
        if (as_val - xs_val) >= ts["MIN_ASXS"]:
            results["ASXS"] = True
    except KeyError: pass

    return results

metrics = {
    "total_observations": 0,
    "callable_bases_sum": 0,
    "marginal_gain_MQ": 0,
    "marginal_gain_BQ": 0,
    "marginal_gain_ASXS": 0,
    "marginal_gain_NM": 0,
    "marginal_gain_SOFTCLIP": 0
}

bam = pysam.AlignmentFile(bam_path, "rb")
bed = pybedtools.BedTool(bed_path).filter(lambda x: x.chrom == chrom)

for interval in bed:
    # NOTE: min_mapping_quality is REMOVED from pileup to allow marginal MQ calculation
    # but MQ is checked inside evaluate_filters to preserve callable_bases_sum logic.
    for pileupcolumn in bam.pileup(chrom, interval.start, interval.end, truncate=True):
        
        observations = []
        
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip: continue
            
            metrics["total_observations"] += 1 
            
            read = pileupread.alignment
            q_pos = pileupread.query_position
            
            # 1. Base Quality Check (Identical to your 'continue' logic)
            bq_pass = read.query_qualities[q_pos] >= ts["MIN_BQ"]
            
            # 2. Run all other health filters
            filter_results = evaluate_filters(read, q_pos, ts)
            
            # 3. An observation 'passes' only if BQ and ALL health filters pass
            all_pass = bq_pass and all(filter_results.values())
            
            observations.append({
                "is_reverse": read.is_reverse,
                "bq_pass": bq_pass,
                "filters": filter_results,
                "all_pass": all_pass
            })

        # --- SITE CALLABILITY (Preserves your original clean_count logic) ---
        passing_reads = [o for o in observations if o["all_pass"]]
        fwd_seen = any(not o["is_reverse"] for o in passing_reads)
        rev_seen = any(o["is_reverse"] for o in passing_reads)
        
        if len(passing_reads) >= ts["MIN_DEPTH"] and fwd_seen and rev_seen:
            # This is the 'final_sum' from your original script
            metrics["callable_bases_sum"] += len(passing_reads)

        # --- MARGINAL COSTS ---
        # Logic: If I removed ONLY this filter, how many more reads would have passed?
        for f_name in ["MQ", "ASXS", "NM", "SOFTCLIP"]:
            gain = 0
            for o in observations:
                if not o["all_pass"]:
                    # Check if it failed ONLY f_name (and passed BQ)
                    other_filters_pass = all(v for k, v in o["filters"].items() if k != f_name)
                    if other_filters_pass and o["bq_pass"]:
                        gain += 1
            metrics[f"marginal_gain_{f_name}"] += gain
            
        # Marginal Gain for BQ
        bq_gain = 0
        for o in observations:
            if not o["all_pass"] and not o["bq_pass"]:
                if all(o["filters"].values()):
                    bq_gain += 1
        metrics["marginal_gain_BQ"] += bq_gain

bam.close()

os.makedirs(os.path.dirname(snakemake.output.json), exist_ok=True)
with open(snakemake.output.json, 'w') as f:
    json.dump(metrics, f)