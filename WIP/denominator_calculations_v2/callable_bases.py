#!/usr/bin/env python3
"""
callable_bases.py — Compute callable base counts for mutation burden analysis.

Appends three columns to an existing sample CSV:
  wes_territory              — bp of exome target after blacklist subtraction
  total_bases_in_wes         — total sequenced bases (all primary reads) landing in WES territory
  total_callable_bases_in_wes — genomic positions where a variant *could* have been called
                                 (filtered depth ∈ [20, 94th-percentile of that BAM])

Read-level discard filters (any one failing → read not counted toward callable depth):
  MAPQ        < 24
  AS − XS     < 22  (alignment score minus sub-optimal alignment score)
  N count     > 2.5 (bases in read equal to 'N')
  NM          > 4.5 (edit distance to reference)
  |insert|    > 410 bp
  soft clips  > 1.5 bases

Usage:
    python callable_bases.py \\
        --input   samples.csv \\
        --wes     exome_targets.bed \\
        --blacklist blacklist.bed \\
        --output  samples_callable.csv \\
        --workers 4          # samples processed in parallel

Dependencies:
    pip install pysam numpy
    (BAM files must be indexed: samtools index *.bam)
"""

import argparse
import csv
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from typing import List, Tuple, Dict, Any

import numpy as np
import pysam

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("callable_bases")


# ===========================================================================
# BED / interval helpers  (no pybedtools dependency)
# ===========================================================================

def parse_bed(path: str) -> List[Tuple[str, int, int]]:
    """Return list of (chrom, start, end) from a BED file (0-based, half-open)."""
    intervals: List[Tuple[str, int, int]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                intervals.append((parts[0], int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    return intervals


def merge_intervals(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Merge overlapping / adjacent intervals within the same chromosome."""
    if not intervals:
        return []
    merged = []
    for chrom, start, end in sorted(intervals, key=lambda x: (x[0], x[1])):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append([chrom, start, end])
    return [tuple(iv) for iv in merged]  # type: ignore[return-value]


def subtract_intervals(
    base: List[Tuple[str, int, int]],
    subtract: List[Tuple[str, int, int]],
) -> List[Tuple[str, int, int]]:
    """
    Efficiently return base intervals with subtract intervals removed.
    Inputs MUST be sorted and merged. Complexity: O(N + M)
    """
    sub_by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in subtract:
        sub_by_chrom[chrom].append((start, end))

    result: List[Tuple[str, int, int]] = []
    
    # Process one chromosome at a time
    base_by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in base:
        base_by_chrom[chrom].append((start, end))

    for chrom in sorted(base_by_chrom.keys()):
        b_intervals = base_by_chrom[chrom]
        s_intervals = sub_by_chrom.get(chrom, [])
        
        b_idx = 0
        s_idx = 0
        
        while b_idx < len(b_intervals):
            b_start, b_end = b_intervals[b_idx]
            
            # If no more blacklist intervals, keep the rest of the base
            if s_idx >= len(s_intervals):
                result.append((chrom, b_start, b_end))
                b_idx += 1
                continue
                
            s_start, s_end = s_intervals[s_idx]
            
            # Case 1: Blacklist is entirely before Base
            if s_end <= b_start:
                s_idx += 1
            # Case 2: Blacklist is entirely after Base
            elif s_start >= b_end:
                result.append((chrom, b_start, b_end))
                b_idx += 1
            # Case 3: Overlap
            else:
                # Part of base before blacklist
                if b_start < s_start:
                    result.append((chrom, b_start, s_start))
                
                # Update base interval to start after the current blacklist interval
                if b_end > s_end:
                    b_intervals[b_idx] = (s_end, b_end)
                    s_idx += 1 # Move to next blacklist to see if it clips more
                else:
                    b_idx += 1 # Base is exhausted
                    
    return result

def territory_size(intervals: List[Tuple[str, int, int]]) -> int:
    return sum(end - start for _, start, end in intervals)


# ---------------------------------------------------------------------------
# Read-level filter  (unchanged interface, but get all tags in one call)
# ---------------------------------------------------------------------------

_FETCH_FLAG_EXCLUDE = 0x4 | 0x100 | 0x400 | 0x800   # unmapped|not-primary|dup|supplementary

def _passes_filters(read: pysam.AlignedSegment) -> bool:
    """
    Return True iff the read passes all quality filters.
    Called ONCE per read (not once per covered position).
    """
    if read.mapping_quality < 24:
        return False

    # Fetch all tags in a single C-level scan → dict
    tags = dict(read.get_tags(with_value_type=False))

    AS = tags.get("AS")
    if AS is not None:
        XS = tags.get("XS", -10_000)
        if AS - XS < 22:
            return False

    seq = read.query_sequence
    if seq and seq.count("N") > 2:          # >2.5 ≡ ≥3 for ints
        return False

    NM = tags.get("NM")
    if NM is not None and NM > 4:           # >4.5 ≡ ≥5 for ints
        return False

    if abs(read.template_length) > 410:
        return False

    if read.cigartuples:
        if sum(l for op, l in read.cigartuples if op == 4) > 1:   # >1.5 ≡ ≥2
            return False

    return True


# ===========================================================================
# Per-sample worker  (fetch-based; each read processed exactly once)
# ===========================================================================

def _process_sample(args: Tuple) -> Dict[str, Any]:
    row, wes_territory_intervals, ter_size = args
    sample_id = row.get("sample_id", "<unknown>")
    bam_path  = row.get("bam_path", "")

    worker_log = logging.getLogger(sample_id)
    worker_log.info("Starting …")

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as exc:
        worker_log.error("Cannot open BAM %s: %s", bam_path, exc)
        return {**row, "wes_territory": ter_size,
                "total_bases_in_wes": "ERROR", "total_callable_bases_in_wes": "ERROR"}

    if not (os.path.exists(bam_path + ".bai") or
            os.path.exists(bam_path.replace(".bam", ".bai"))):
        worker_log.warning("No BAM index for %s — pileup may fail", bam_path)

    # One uint32 slot per bp in the WES territory (for the 94th-pct calculation)
    filtered_depths  = np.zeros(ter_size, dtype=np.uint32)
    total_bases_in_wes: int = 0
    offset: int = 0

    for chrom, reg_start, reg_end in wes_territory_intervals:
        reg_len = reg_end - reg_start

        # Diff arrays: +1 at read start, -1 at read end → cumsum = depth profile
        # Using int32 because depths won't exceed ~32 k and we need signed for diffs
        all_diff  = np.zeros(reg_len + 1, dtype=np.int32)
        filt_diff = np.zeros(reg_len + 1, dtype=np.int32)

        try:
            for read in bam.fetch(chrom, reg_start, reg_end):
                # Exclude unmapped / secondary / duplicate / supplementary at flag level
                if read.flag & _FETCH_FLAG_EXCLUDE:
                    continue
                ref_end = read.reference_end        # None for fully-soft-clipped reads
                if ref_end is None:
                    continue

                # Clamp read span to this region (fetch can return overhanging reads)
                s = max(read.reference_start, reg_start) - reg_start
                e = min(ref_end,              reg_end)   - reg_start
                if s >= e:
                    continue

                # Record this read's contribution to the all-reads depth
                all_diff[s] += 1
                all_diff[e] -= 1

                # Record to filtered depth only if it passes QC (called ONCE per read)
                if _passes_filters(read):
                    filt_diff[s] += 1
                    filt_diff[e] -= 1

        except (ValueError, KeyError) as exc:
            worker_log.warning("Skipping %s:%d-%d — %s", chrom, reg_start, reg_end, exc)

        # Convert diff arrays → depth profiles via cumulative sum (fast numpy op)
        all_depth  = np.cumsum(all_diff [:reg_len]).astype(np.uint32)
        filt_depth = np.cumsum(filt_diff[:reg_len]).astype(np.uint32)

        total_bases_in_wes += int(all_depth.sum())
        filtered_depths[offset : offset + reg_len] = filt_depth
        offset += reg_len

    bam.close()

    p94 = float(np.percentile(filtered_depths, 94))
    # 2. Create a boolean mask for positions that meet your criteria
    # We only care about sites where 20 <= depth <= p94
    mask = (filtered_depths >= 20) & (filtered_depths <= p94)

    # 3. Sum the actual depth values at those valid positions
    # This gives you the total count of "trusted" bases (nucleotides)
    callable_bases_sum = int(np.sum(filtered_depths[mask]))

    worker_log.info(
        "Done. total_bases=%d  p94_threshold=%.1f  callable_sum=%d",
        total_bases_in_wes, p94, callable_bases_sum,
    )
    
    return {
        **row,
        "wes_territory":              ter_size,
        "total_bases_in_wes":         total_bases_in_wes,
        "total_callable_bases_in_wes": callable_bases_sum,
    }

# ===========================================================================
# Main
# ===========================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True, help="Input CSV (must have a bam_path column)")
    p.add_argument("--wes", required=True, help="WES target BED file (0-based)")
    p.add_argument("--blacklist", required=True, help="Blacklist BED file (0-based)")
    p.add_argument("--output", required=True, help="Output CSV path")
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of samples to process in parallel (default: 1)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ------------------------------------------------------------------
    # Build WES territory (WES minus blacklist)
    # ------------------------------------------------------------------
    log.info("Loading WES BED: %s", args.wes)
    wes_raw = parse_bed(args.wes)
    wes_merged = merge_intervals(wes_raw)

    log.info("Loading blacklist BED: %s", args.blacklist)
    bl_raw = parse_bed(args.blacklist)
    bl_merged = merge_intervals(bl_raw)

    log.info("Subtracting blacklist from WES territory …")
    wes_final = subtract_intervals(wes_merged, bl_merged)
    ter_size = territory_size(wes_final)
    log.info(
        "WES territory: %d regions, %d bp after blacklist removal",
        len(wes_final),
        ter_size,
    )

    # ------------------------------------------------------------------
    # Read input CSV
    # ------------------------------------------------------------------
    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if not rows:
        log.error("Input CSV is empty.")
        sys.exit(1)

    required = {"bam_path"}
    missing = required - set(fieldnames)
    if missing:
        log.error("Input CSV is missing required columns: %s", missing)
        sys.exit(1)

    new_cols = ["wes_territory", "total_bases_in_wes", "total_callable_bases_in_wes"]
    out_fieldnames = fieldnames + [c for c in new_cols if c not in fieldnames]

    # ------------------------------------------------------------------
    # Process samples (parallel or serial)
    # ------------------------------------------------------------------
    worker_args = [(row, wes_final, ter_size) for row in rows]
    results: List[Dict[str, Any]] = [{}] * len(rows)

    # Map sample_id → original index so we can preserve CSV order
    id_to_idx = {row.get("sample_id", str(i)): i for i, row in enumerate(rows)}

    if args.workers > 1:
        log.info("Processing %d samples with %d parallel workers …", len(rows), args.workers)
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_to_idx = {
                executor.submit(_process_sample, wa): i
                for i, wa in enumerate(worker_args)
            }
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    results[idx] = future.result()
                except Exception as exc:
                    log.error("Sample at index %d failed: %s", idx, exc)
                    results[idx] = {
                        **rows[idx],
                        "wes_territory": ter_size,
                        "total_bases_in_wes": "ERROR",
                        "total_callable_bases_in_wes": "ERROR",
                    }
    else:
        log.info("Processing %d samples serially …", len(rows))
        for i, wa in enumerate(worker_args):
            results[i] = _process_sample(wa)

    # ------------------------------------------------------------------
    # Write output CSV
    # ------------------------------------------------------------------
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    log.info("Results written to: %s", args.output)


if __name__ == "__main__":
    main()
