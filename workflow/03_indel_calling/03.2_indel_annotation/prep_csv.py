import pandas as pd
import numpy as np
import os
import gzip
import sys
import time

# -----------------------
# User Configuration
# -----------------------
INPUT_CSV = "../Sample_csv/pre_processed_lynch_samples.csv"      # Input: Must have 'mean_dcs_depth' & 'vcf_path'
OUTPUT_CSV = "../Sample_csv/processed_lynch_samples.csv"     # Output: Will include calculated 'bases'
REPORT_FILE = "lynch_samples_preprocessing_report.txt"

# Geography
EXOME_BED = "/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/reference_genome/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed"
BLACKLIST_FILES = ["master_blacklist.bed"]

# Quality Heuristic Settings
# We check the first N variants to estimate the "messiness" of the sample
VCF_SAMPLE_LIMIT = 50000 
MIN_MAPQ_HEURISTIC = 40.0   # Reads below this are considered "uncallable" noise
MIN_ASXS_HEURISTIC = 50.0  # Alignment score threshold

# -----------------------
# 1. BED Interval Logic
# -----------------------
def load_and_merge_intervals(bed_paths):
    """Loads one or multiple BED files and merges overlapping intervals."""
    intervals = []
    paths = [bed_paths] if isinstance(bed_paths, str) else bed_paths
    
    for path in paths:
        if not os.path.exists(path):
            print(f"[WARN] BED file not found: {path}. Skipping.", file=sys.stderr)
            continue
            
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split()
                if len(parts) < 3: continue
                # 0-based start, 1-based end
                intervals.append((parts[0], int(parts[1]), int(parts[2])))
    
    intervals.sort()
    
    merged = []
    if not intervals: return merged
        
    curr_chrom, curr_start, curr_end = intervals[0]
    for i in range(1, len(intervals)):
        next_chrom, next_start, next_end = intervals[i]
        if next_chrom == curr_chrom and next_start < curr_end:
            curr_end = max(curr_end, next_end)
        else:
            merged.append((curr_chrom, curr_start, curr_end))
            curr_chrom, curr_start, curr_end = next_chrom, next_start, next_end
            
    merged.append((curr_chrom, curr_start, curr_end))
    return merged

def calculate_bed_length(intervals):
    return sum(end - start for _, start, end in intervals)

def subtract_intervals(target_intervals, subtract_intervals):
    """Subtracts 'subtract_intervals' from 'target_intervals'."""
    final_intervals = []
    targets_by_chrom = {}
    for chrom, start, end in target_intervals:
        targets_by_chrom.setdefault(chrom, []).append((start, end))
        
    sub_by_chrom = {}
    for chrom, start, end in subtract_intervals:
        sub_by_chrom.setdefault(chrom, []).append((start, end))
        
    for chrom, t_regions in targets_by_chrom.items():
        if chrom not in sub_by_chrom:
            final_intervals.extend([(chrom, s, e) for s, e in t_regions])
            continue
            
        s_regions = sub_by_chrom[chrom]
        for t_start, t_end in t_regions:
            current_pos = t_start
            for s_start, s_end in s_regions:
                if s_end <= current_pos: continue 
                if s_start >= t_end: break        
                
                if s_start > current_pos:
                    final_intervals.append((chrom, current_pos, s_start))
                current_pos = max(current_pos, s_end)
                
            if current_pos < t_end:
                final_intervals.append((chrom, current_pos, t_end))
    return final_intervals

# -----------------------
# 2. VCF Quality Sampling
# -----------------------
def estimate_sample_quality(vcf_path):
    """
    Parses the VCF to calculate the fraction of variants that PASS 
    basic mapping/alignment thresholds.
    Returns: float (0.0 to 1.0)
    """
    if not os.path.exists(vcf_path):
        return 0.90 # Safe default if file missing
        
    total_checked = 0
    passed_checked = 0
    
    openf = gzip.open if vcf_path.endswith(".gz") else open
    
    try:
        with openf(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"): continue
                
                # Info extraction
                info_part = line.split("\t")[7]
                info_dict = {}
                # Quick parse of just the fields we need to save time
                for item in info_part.split(";"):
                    if item.startswith("MEAN_ALT_MAPQ=") or item.startswith("AVG_ALT_ASXS="):
                        k, v = item.split("=")
                        info_dict[k] = float(v)
                
                # Check metrics (defaults if missing in VCF line)
                mapq = info_dict.get("MEAN_ALT_MAPQ", info_dict.get("MEAN_MAPQ", 60.0))
                asxs = info_dict.get("AVG_ALT_ASXS", 150.0)
                
                total_checked += 1
                if mapq >= MIN_MAPQ_HEURISTIC and asxs >= MIN_ASXS_HEURISTIC:
                    passed_checked += 1
                
                if total_checked >= VCF_SAMPLE_LIMIT:
                    break
    except Exception as e:
        print(f"[WARN] Error parsing {vcf_path}: {e}")
        return 0.90

    if total_checked == 0:
        return 0.95 # Clean file or empty, assume high quality
        
    return passed_checked / total_checked

# -----------------------
# Main Execution
# -----------------------
if __name__ == "__main__":
    start_time = time.time()
    report_lines = []
    
    def log(msg):
        print(msg)
        report_lines.append(msg)

    log("=== PRE-PROCESSING & CALIBRATION REPORT ===")
    
    # 1. Load Data
    if not os.path.exists(INPUT_CSV):
        print(f"[ERROR] Input file {INPUT_CSV} not found.")
        sys.exit(1)
        
    df = pd.read_csv(INPUT_CSV)
    required_cols = ['sample_id', 'vcf_path', 'mean_dcs_depth']
    if not all(col in df.columns for col in required_cols):
        print(f"[ERROR] CSV missing required columns: {required_cols}")
        sys.exit(1)

    # 2. Geography Math
    log("\n[1] CALCULATING EFFECTIVE EXOME SIZE")
    exome_intervals = load_and_merge_intervals(EXOME_BED)
    raw_len = calculate_bed_length(exome_intervals)
    log(f"   Raw Exome Size:        {raw_len:,.0f} bp")

    blacklist_intervals = load_and_merge_intervals(BLACKLIST_FILES)
    bl_len = calculate_bed_length(blacklist_intervals)
    log(f"   Total Blacklist Area:  {bl_len:,.0f} bp")

    effective_intervals = subtract_intervals(exome_intervals, blacklist_intervals)
    effective_len = calculate_bed_length(effective_intervals)
    retention = (effective_len / raw_len * 100) if raw_len > 0 else 0
    
    log(f"   Effective Exome Size:  {effective_len:,.0f} bp")
    log(f"   Retention Rate:        {retention:.2f}%")

    # 3. Per-Sample Calibration
    log("\n[2] CALIBRATING PER-SAMPLE BASES")
    log(f"   Heuristic: Fraction of variants with MAPQ >= {MIN_MAPQ_HEURISTIC} & ASXS >= {MIN_ASXS_HEURISTIC}")
    
    new_rows = []
    
    print("   Processing samples...", end="", flush=True)
    
    for idx, row in df.iterrows():
        if idx % 10 == 0: print(".", end="", flush=True)
        
        # Calculate Quality Factor
        q_factor = estimate_sample_quality(row['vcf_path'])
        
        # The Core Formula
        # Bases = Effective_Bp * Mean_Depth * Quality_Factor
        # We round to int because 'bases' implies discrete count
        calc_bases = int(effective_len * row['mean_dcs_depth'] * q_factor)
        
        # Update row
        row_data = row.to_dict()
        row_data['bases'] = calc_bases
        row_data['quality_factor'] = round(q_factor, 4) # Keeping this for reference
        row_data['effective_territory'] = effective_len
        new_rows.append(row_data)

    print(" Done.")
    
    result_df = pd.DataFrame(new_rows)
    
    # 4. Save and Report
    cols_order = ['sample_id', 'sample_type', 'age', 'mean_dcs_depth', 'quality_factor', 'bases', 'vcf_path']
    # Add any extra columns that were in the original but not in our explicit list
    cols_order += [c for c in result_df.columns if c not in cols_order]
    
    result_df = result_df[cols_order]
    result_df.to_csv(OUTPUT_CSV, index=False)
    
    log(f"\n[3] SUMMARY")
    log(f"   Processed {len(result_df)} samples.")
    log(f"   Average Quality Factor: {result_df['quality_factor'].mean():.4f}")
    log(f"   Lowest Quality Factor:  {result_df['quality_factor'].min():.4f} (Sample: {result_df.sort_values('quality_factor').iloc[0]['sample_id']})")
    log(f"   Results saved to:       {OUTPUT_CSV}")
    log(f"   Full report saved to:   {REPORT_FILE}")

    with open(REPORT_FILE, "w") as f:
        f.write("\n".join(report_lines))

    print(f"\n[SUCCESS] Pre-processing complete. Ready for tuning script.")