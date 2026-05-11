#!/usr/bin/env python3
"""
callable_bases.py — Compute callable base counts for mutation burden analysis.

Appends five columns to an existing sample CSV:
  wes_territory                      — bp of exome target after blacklist subtraction
  total_bases_in_wes                 — total sequenced bases (all primary reads)
  total_callable_bases_in_wes        — standard filters (depth 20 to p95)
  total_callable_bases_in_wes_strict — 10% stricter filters (depth 22 to p94)
  total_callable_bases_in_wes_lax    — 10% more lax filters (depth 18 to p97)

Usage:
    python callable_bases.py \
        --input   samples.csv \
        --wes     exome_targets.bed \
        --blacklist blacklist.bed \
        --output  samples_callable.csv \
        --workers 4
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
    if not intervals:
        return []
    merged = []
    for chrom, start, end in sorted(intervals, key=lambda x: (x[0], x[1])):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append([chrom, start, end])
    return [tuple(iv) for iv in merged]


def subtract_intervals(
    base: List[Tuple[str, int, int]],
    subtract: List[Tuple[str, int, int]],
) -> List[Tuple[str, int, int]]:
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
# Read-level filter logic (Sensitivity Analysis)
# ---------------------------------------------------------------------------

_FETCH_FLAG_EXCLUDE = 0x4 | 0x100 | 0x400 | 0x800   # unmapped|not-primary|dup|supplementary

def _get_filter_status(read: pysam.AlignedSegment) -> Tuple[bool, bool, bool]:
    """
    Returns (passes_strict, passes_std, passes_lax) based on ~10% adjustments.
    """
    strict = std = lax = True
    
    # 1. MAPQ (Standard: 24)
    mapq = read.mapping_quality
    if mapq < 21: lax = False
    if mapq < 24: std = False
    if mapq < 26: strict = False
    
    # Fast exit if it fails the most permissive tier
    if not lax and not std and not strict:
        return False, False, False

    # 2. Alignment Tags
    tags = dict(read.get_tags(with_value_type=False))
    
    AS = tags.get("AS")
    if AS is not None:
        XS = tags.get("XS", -10_000)
        diff = AS - XS
        # Standard: 22
        if diff < 19: lax = False
        if diff < 22: std = False
        if diff < 24: strict = False

    # 3. N Count (Standard: <= 2)
    seq = read.query_sequence
    if seq:
        n_count = seq.count("N")
        if n_count > 3: lax = False
        if n_count > 2: std = False
        if n_count > 1: strict = False

    # 4. Edit distance (Standard: <= 4)
    NM = tags.get("NM")
    if NM is not None:
        if NM > 5: lax = False
        if NM > 4: std = False
        if NM > 3: strict = False

    # 5. Insert size (Standard: <= 410)
    ins_len = abs(read.template_length)
    if ins_len > 451: lax = False
    if ins_len > 410: std = False
    if ins_len > 369: strict = False

    # 6. Soft clips (Standard: <= 1)
    if read.cigartuples:
        sc = sum(l for op, l in read.cigartuples if op == 4)
        if sc > 2: lax = False
        if sc > 1: std = False
        if sc > 0: strict = False

    # Safety enforcement: A read cannot pass strict if it fails standard/lax
    std = std and lax
    strict = strict and std
    
    return strict, std, lax


# ===========================================================================
# Per-sample worker
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
        return {
            **row, "wes_territory": ter_size,
            "total_bases_in_wes": "ERROR", 
            "total_callable_bases_in_wes": "ERROR",
            "total_callable_bases_in_wes_strict": "ERROR",
            "total_callable_bases_in_wes_lax": "ERROR"
        }

    if not (os.path.exists(bam_path + ".bai") or
            os.path.exists(bam_path.replace(".bam", ".bai"))):
        worker_log.warning("No BAM index for %s — pileup may fail", bam_path)

    # 3 depth profiles for sensitivity analysis
    depths_strict = np.zeros(ter_size, dtype=np.uint32)
    depths_std    = np.zeros(ter_size, dtype=np.uint32)
    depths_lax    = np.zeros(ter_size, dtype=np.uint32)
    
    total_bases_in_wes: int = 0
    offset: int = 0
    TRIM_BP = 4

    for chrom, reg_start, reg_end in wes_territory_intervals:
        reg_len = reg_end - reg_start

        all_diff  = np.zeros(reg_len + 1, dtype=np.int32)
        diff_strict = np.zeros(reg_len + 1, dtype=np.int32)
        diff_std    = np.zeros(reg_len + 1, dtype=np.int32)
        diff_lax    = np.zeros(reg_len + 1, dtype=np.int32)

        try:
            for read in bam.fetch(chrom, reg_start, reg_end):
                if read.flag & _FETCH_FLAG_EXCLUDE:
                    continue
                ref_end = read.reference_end
                if ref_end is None:
                    continue

                # 1. TOTAL BASES
                s_all = max(read.reference_start, reg_start) - reg_start
                e_all = min(ref_end,              reg_end)   - reg_start
                if s_all < e_all:
                    all_diff[s_all] += 1
                    all_diff[e_all] -= 1

                # 2. CALLABLE BASES (With trimming & sensitivity checks)
                adj_start = read.reference_start + TRIM_BP
                adj_end   = ref_end - TRIM_BP

                s_filt = max(adj_start, reg_start) - reg_start
                e_filt = min(adj_end,   reg_end)   - reg_start
                
                if s_filt < e_filt:
                    p_strict, p_std, p_lax = _get_filter_status(read)
                    if p_strict:
                        diff_strict[s_filt] += 1
                        diff_strict[e_filt] -= 1
                    if p_std:
                        diff_std[s_filt] += 1
                        diff_std[e_filt] -= 1
                    if p_lax:
                        diff_lax[s_filt] += 1
                        diff_lax[e_filt] -= 1

        except (ValueError, KeyError) as exc:
            worker_log.warning("Skipping %s:%d-%d — %s", chrom, reg_start, reg_end, exc)

        # Apply cumulative sum
        all_depth  = np.cumsum(all_diff [:reg_len]).astype(np.uint32)
        total_bases_in_wes += int(all_depth.sum())

        depths_strict[offset : offset + reg_len] = np.cumsum(diff_strict[:reg_len]).astype(np.uint32)
        depths_std[offset : offset + reg_len]    = np.cumsum(diff_std[:reg_len]).astype(np.uint32)
        depths_lax[offset : offset + reg_len]    = np.cumsum(diff_lax[:reg_len]).astype(np.uint32)
        
        offset += reg_len

    bam.close()

    # Calculate different percentiles dynamically for each array
    p_strict = float(np.percentile(depths_strict, 94))
    p_std    = float(np.percentile(depths_std, 95))
    p_lax    = float(np.percentile(depths_lax, 97))

    # Masking thresholds
    mask_strict = (depths_strict >= 22) & (depths_strict <= p_strict)
    mask_std    = (depths_std >= 20) & (depths_std <= p_std)
    mask_lax    = (depths_lax >= 18) & (depths_lax <= p_lax)

    callable_strict = int(np.sum(depths_strict[mask_strict]))
    callable_std    = int(np.sum(depths_std[mask_std]))
    callable_lax    = int(np.sum(depths_lax[mask_lax]))

    worker_log.info(
        "Done. Total=%d  | Std Callable=%d (p95=%.1f) | Strict Callable=%d (p94=%.1f) | Lax Callable=%d (p97=%.1f)",
        total_bases_in_wes, 
        callable_std, p_std, 
        callable_strict, p_strict, 
        callable_lax, p_lax
    )
    
    return {
        **row,
        "wes_territory":                      ter_size,
        "total_bases_in_wes":                 total_bases_in_wes,
        "total_callable_bases_in_wes":        callable_std,
        "total_callable_bases_in_wes_strict": callable_strict,
        "total_callable_bases_in_wes_lax":    callable_lax,
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
        len(wes_final), ter_size,
    )

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
        "wes_territory", 
        "total_bases_in_wes", 
        "total_callable_bases_in_wes",
        "total_callable_bases_in_wes_strict",
        "total_callable_bases_in_wes_lax"
    ]
    out_fieldnames = fieldnames + [c for c in new_cols if c not in fieldnames]

    worker_args = [(row, wes_final, ter_size) for row in rows]
    results: List[Dict[str, Any]] = [{}] * len(rows)

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
                        "total_callable_bases_in_wes_strict": "ERROR",
                        "total_callable_bases_in_wes_lax": "ERROR",
                    }
    else:
        log.info("Processing %d samples serially …", len(rows))
        for i, wa in enumerate(worker_args):
            results[i] = _process_sample(wa)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    log.info("Results written to: %s", args.output)


if __name__ == "__main__":
    main()