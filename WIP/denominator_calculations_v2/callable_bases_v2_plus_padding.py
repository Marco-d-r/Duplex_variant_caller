#!/usr/bin/env python3
"""
callable_bases.py — Compute callable base counts for mutation burden analysis.

Appends six columns to an existing sample CSV:
  wes_territory                      — bp of exome target after blacklist subtraction
  total_bases_in_wes                 — total sequenced bases landing in WES territory
  total_callable_bases_in_wes        — genomic positions where a variant *could* have been called
  padded_wes_territory               — bp of exome target + padding after blacklist subtraction
  total_bases_in_padded_wes          — total sequenced bases in padded WES territory
  total_callable_bases_in_padded_wes — genomic positions where a variant *could* have been called in padded WES

Read-level discard filters (any one failing → read not counted toward callable depth):
  MAPQ        < 24
  AS − XS     < 22  (alignment score minus sub-optimal alignment score)
  N count     > 2.5 (bases in read equal to 'N')
  NM          > 4.5 (edit distance to reference)
  |insert|    > 410 bp
  soft clips  > 1.5 bases

Additional Constraints:
  The first 4 and last 4 bases of any given read are considered uncallable 
  and are trimmed from the callable depth calculations (but kept for total bases).

Usage:
    python callable_bases.py \
        --input   samples.csv \
        --wes     exome_targets.bed \
        --blacklist blacklist.bed \
        --output  samples_callable.csv \
        --padding 135 \
        --workers 4          

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


def pad_intervals(intervals: List[Tuple[str, int, int]], pad_bp: int) -> List[Tuple[str, int, int]]:
    """Add padding to both sides of each interval, clamping the start to 0."""
    return [(chrom, max(0, start - pad_bp), end + pad_bp) for chrom, start, end in intervals]


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
            
            if s_idx >= len(s_intervals):
                result.append((chrom, b_start, b_end))
                b_idx += 1
                continue
                
            s_start, s_end = s_intervals[s_idx]
            
            if s_end <= b_start:
                s_idx += 1
            elif s_start >= b_end:
                result.append((chrom, b_start, b_end))
                b_idx += 1
            else:
                if b_start < s_start:
                    result.append((chrom, b_start, s_start))
                
                if b_end > s_end:
                    b_intervals[b_idx] = (s_end, b_end)
                    s_idx += 1
                else:
                    b_idx += 1
                    
    return result

def territory_size(intervals: List[Tuple[str, int, int]]) -> int:
    return sum(end - start for _, start, end in intervals)


# ---------------------------------------------------------------------------
# Read-level filter  
# ---------------------------------------------------------------------------

_FETCH_FLAG_EXCLUDE = 0x4 | 0x100 | 0x400 | 0x800   # unmapped|not-primary|dup|supplementary

def _passes_filters(read: pysam.AlignedSegment) -> bool:
    if read.mapping_quality < 24:
        return False

    tags = dict(read.get_tags(with_value_type=False))

    AS = tags.get("AS")
    if AS is not None:
        XS = tags.get("XS", -10_000)
        if AS - XS < 22:
            return False

    seq = read.query_sequence
    if seq and seq.count("N") > 2:
        return False

    NM = tags.get("NM")
    if NM is not None and NM > 4:
        return False

    if abs(read.template_length) > 410:
        return False

    if read.cigartuples:
        if sum(l for op, l in read.cigartuples if op == 4) > 1:
            return False

    return True


# ===========================================================================
# Per-sample worker 
# ===========================================================================

def _process_sample(args: Tuple) -> Dict[str, Any]:
    row, wes_final, ter_size, padded_final, padded_ter_size = args
    sample_id = row.get("sample_id", "<unknown>")
    bam_path  = row.get("bam_path", "")

    worker_log = logging.getLogger(sample_id)
    worker_log.info("Starting …")

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as exc:
        worker_log.error("Cannot open BAM %s: %s", bam_path, exc)
        return {
            **row, "wes_territory": ter_size, "total_bases_in_wes": "ERROR", "total_callable_bases_in_wes": "ERROR",
            "padded_wes_territory": padded_ter_size, "total_bases_in_padded_wes": "ERROR", "total_callable_bases_in_padded_wes": "ERROR"
        }

    if not (os.path.exists(bam_path + ".bai") or os.path.exists(bam_path.replace(".bam", ".bai"))):
        worker_log.warning("No BAM index for %s — pileup may fail", bam_path)

    # -----------------------------------------------------------------------
    # Reusable core metric function applied to any interval list
    # -----------------------------------------------------------------------
    def _compute_metrics(intervals: List[Tuple[str, int, int]], size: int) -> Tuple[int, int, float]:
        filtered_depths = np.zeros(size, dtype=np.uint32)
        total_bases = 0
        offset = 0
        TRIM_BP = 4

        for chrom, reg_start, reg_end in intervals:
            reg_len = reg_end - reg_start
            all_diff  = np.zeros(reg_len + 1, dtype=np.int32)
            filt_diff = np.zeros(reg_len + 1, dtype=np.int32)

            try:
                for read in bam.fetch(chrom, reg_start, reg_end):
                    if read.flag & _FETCH_FLAG_EXCLUDE:
                        continue
                    ref_end = read.reference_end
                    if ref_end is None:
                        continue

                    # TOTAL BASES
                    s_all = max(read.reference_start, reg_start) - reg_start
                    e_all = min(ref_end,              reg_end)   - reg_start
                    if s_all < e_all:
                        all_diff[s_all] += 1
                        all_diff[e_all] -= 1

                    # CALLABLE BASES 
                    adj_start = read.reference_start + TRIM_BP
                    adj_end   = ref_end - TRIM_BP

                    s_filt = max(adj_start, reg_start) - reg_start
                    e_filt = min(adj_end,   reg_end)   - reg_start
                    
                    if s_filt < e_filt:
                        if _passes_filters(read):
                            filt_diff[s_filt] += 1
                            filt_diff[e_filt] -= 1

            except (ValueError, KeyError) as exc:
                worker_log.debug("Skipping %s:%d-%d — %s", chrom, reg_start, reg_end, exc)

            all_depth  = np.cumsum(all_diff [:reg_len]).astype(np.uint32)
            filt_depth = np.cumsum(filt_diff[:reg_len]).astype(np.uint32)

            total_bases += int(all_depth.sum())
            filtered_depths[offset : offset + reg_len] = filt_depth
            offset += reg_len
        
        # Calculate percentiles and sum for this specific territory
        if size > 0:
            p94 = float(np.percentile(filtered_depths, 94))
            mask = (filtered_depths >= 20) & (filtered_depths <= p94)
            callable_sum = int(np.sum(filtered_depths[mask]))
        else:
            p94 = 0.0
            callable_sum = 0
            
        return total_bases, callable_sum, p94

    # 1. Compute metrics for strict WES
    wes_tot, wes_call, wes_p94 = _compute_metrics(wes_final, ter_size)
    
    # 2. Compute metrics for padded WES
    pad_tot, pad_call, pad_p94 = _compute_metrics(padded_final, padded_ter_size)

    bam.close()

    worker_log.info(
        "Done. WES bases=%d (call=%d, p94=%.1f) | Padded WES bases=%d (call=%d, p94=%.1f)",
        wes_tot, wes_call, wes_p94, pad_tot, pad_call, pad_p94
    )
    
    return {
        **row,
        "wes_territory":                      ter_size,
        "total_bases_in_wes":                 wes_tot,
        "total_callable_bases_in_wes":        wes_call,
        "padded_wes_territory":               padded_ter_size,
        "total_bases_in_padded_wes":          pad_tot,
        "total_callable_bases_in_padded_wes": pad_call,
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
    p.add_argument("--padding", type=int, default=135, help="Padding in bp applied to both sides of WES targets (default: 135)")
    p.add_argument("--workers", type=int, default=1, help="Number of samples to process in parallel (default: 1)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ------------------------------------------------------------------
    # Build WES territories (Standard and Padded)
    # ------------------------------------------------------------------
    log.info("Loading blacklist BED: %s", args.blacklist)
    bl_raw = parse_bed(args.blacklist)
    bl_merged = merge_intervals(bl_raw)

    log.info("Loading WES BED: %s", args.wes)
    wes_raw = parse_bed(args.wes)
    
    # 1. Standard WES Territory
    wes_merged = merge_intervals(wes_raw)
    wes_final = subtract_intervals(wes_merged, bl_merged)
    ter_size = territory_size(wes_final)
    
    # 2. Padded WES Territory
    log.info(f"Generating padded territory (+{args.padding} bp) …")
    padded_raw = pad_intervals(wes_raw, args.padding)
    padded_merged = merge_intervals(padded_raw)
    padded_final = subtract_intervals(padded_merged, bl_merged)
    padded_ter_size = territory_size(padded_final)
    
    log.info(
        "WES territory: %d regions, %d bp after blacklist removal",
        len(wes_final), ter_size,
    )
    log.info(
        "Padded WES territory: %d regions, %d bp after blacklist removal",
        len(padded_final), padded_ter_size,
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

    new_cols = [
        "wes_territory", "total_bases_in_wes", "total_callable_bases_in_wes",
        "padded_wes_territory", "total_bases_in_padded_wes", "total_callable_bases_in_padded_wes"
    ]
    out_fieldnames = fieldnames + [c for c in new_cols if c not in fieldnames]

    # ------------------------------------------------------------------
    # Process samples 
    # ------------------------------------------------------------------
    worker_args = [(row, wes_final, ter_size, padded_final, padded_ter_size) for row in rows]
    results: List[Dict[str, Any]] = [{}] * len(rows)

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
                        "wes_territory": ter_size, "total_bases_in_wes": "ERROR", "total_callable_bases_in_wes": "ERROR",
                        "padded_wes_territory": padded_ter_size, "total_bases_in_padded_wes": "ERROR", "total_callable_bases_in_padded_wes": "ERROR"
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