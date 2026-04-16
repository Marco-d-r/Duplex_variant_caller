#!/usr/bin/env python3
"""
callable_bases_calculator.py
─────────────────────────────────────────────────────────────────────────────
Computes the number of *callable sequenced bases* for each duplex sequencing
sample, so that a mutation rate can be expressed as:

    mutation_rate = variant_count / callable_sequenced_bases

"Callable sequenced bases" is the sum of per-position DCS read depth across
all positions that (a) lie inside the target BED, (b) are not blacklisted,
and (c) reach at least --min_depth coverage after all read-level filters have
been applied.  In other words, if 100 reads each 150 bp long cover a region
and all pass every filter, those reads contribute 15,000 bases to the
denominator — not 150 genomic positions.

Read-level filters (applied before depth is computed, matching the variant
caller exactly):
  • Unmapped, non-primary, QC-failed, duplicate, and supplementary reads
    excluded via samtools FLAG filter -F 3844
  • Mapping quality ≥ --min_mapq  (default 40)
  • AS − XS ≥ --min_asxs          (default 50; uniquely mapping reads that
    lack an XS tag are treated as XS = 0 and therefore always pass)

Blacklisted positions are subtracted from both BEDs before any depth
calculation, so they never contribute reads to the denominator.

Output CSV preserves every column from the input CSV and appends:

  raw_total_bases                overall sequenced bases in the input BAM
                                 (no filters, no blacklist)
  total_wes_bp                   genomic positions in WES BED after blacklist
                                 subtraction (reference only — not used as the
                                 mutation-rate denominator)
  total_sequenced_bases_wes      Σ depth across all non-blacklisted WES
                                 positions (all reads that pass filters)
  callable_sequenced_bases_wes   Σ depth at positions where depth ≥ min_depth
                                 ← use this as the mutation-rate denominator
  pct_callable_bases_wes         callable / total × 100
  mean_depth_wes                 mean per-position depth over the WES BED

  (identical set of columns repeated for the padded BED)

Usage
─────
python callable_bases_calculator.py \
    --input_csv   samples.csv \
    --bam_col     bam_path \
    --wes_bed     capture_panel.bed \
    --blacklist   hg38_blacklist.bed \
    --genome_fai  hg38.fa.fai \
    --output_csv  callable_bases_results.csv \
    [--padded_bed already_padded.bed] \
    [--padding    100] \
    [--min_mapq   40] \
    [--min_asxs   50] \
    [--min_depth  1] \
    [--threads    8] \
    [--verbose]

Dependencies
────────────
  samtools ≥ 1.14   (required for -e expression filter support)
  bedtools ≥ 2.29
  Python   ≥ 3.8
"""

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile

# ─── logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ─── utilities ────────────────────────────────────────────────────────────────

def check_tool(name: str) -> None:
    """Abort early if a required command-line tool is not on PATH."""
    if shutil.which(name) is None:
        log.error("Required tool '%s' not found on PATH. Please install it.", name)
        sys.exit(1)


def run(cmd: str, check: bool = True, capture: bool = False) -> subprocess.CompletedProcess:
    """Run a shell command; forward stderr to the debug logger."""
    log.debug("CMD: %s", cmd)
    result = subprocess.run(
        cmd,
        shell=True,
        check=check,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            log.debug("  STDERR: %s", line)
    return result


# ─── BED helpers ──────────────────────────────────────────────────────────────

def bed_size(bed_path: str) -> int:
    """Return the total number of base-pairs spanned by a BED file."""
    total = 0
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            total += int(parts[2]) - int(parts[1])
    return total


def clean_bed(wes_bed: str, blacklist: str, out_path: str) -> None:
    """
    Subtract blacklisted regions from wes_bed, then sort and merge intervals.
    Positions in the blacklist will not appear in the depth output and
    therefore cannot contribute any reads to the callable-bases denominator.
    """
    cmd = (
        f"bedtools subtract -a {wes_bed} -b {blacklist} "
        f"| bedtools sort -i - "
        f"| bedtools merge -i - "
        f"> {out_path}"
    )
    run(cmd)


def make_padded_bed(
    wes_bed: str,
    genome_fai: str,
    padding: int,
    blacklist: str,
    out_path: str,
) -> None:
    """
    Expand every interval in wes_bed by ±padding bp (clamped to chromosome
    ends), subtract blacklisted regions, then sort and merge.
    """
    cmd = (
        f"bedtools slop -i {wes_bed} -g {genome_fai} -b {padding} "
        f"| bedtools subtract -a - -b {blacklist} "
        f"| bedtools sort -i - "
        f"| bedtools merge -i - "
        f"> {out_path}"
    )
    run(cmd)


# ─── depth computation ────────────────────────────────────────────────────────

def get_raw_total_bases(bam: str, threads: int) -> int:
    """
    Compute the overall number of sequenced bases in the BAM without any filters.
    Utilizes samtools stats and parses the 'total length' metric.
    """
    cmd = f"samtools stats -@ {threads} {bam}"
    log.debug("Raw total bases pipeline: %s", cmd)
    proc = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        check=True
    )
    for line in proc.stdout.splitlines():
        if line.startswith("SN\ttotal length:"):
            # SN format: SN\ttotal length:\t<number>
            return int(line.split("\t")[2])
    return 0


def compute_depth_stats(
    bam: str,
    region_bed: str,
    min_mapq: int,
    min_asxs: int,
    min_depth: int,
    threads: int,
) -> dict:
    """
    Compute callable sequenced base statistics for one BAM over one target BED.
    """
    asxs_expr = f"(exists([XS]) && ([AS] - [XS] >= {min_asxs})) || (!exists([XS]) && ([AS] >= {min_asxs}))"

    depth_cmd = (
        f"samtools view -@ {threads} -b "
        f"-F 3844 "
        f"-q {min_mapq} "
        f"-e '{asxs_expr}' "
        f"{bam} "
        f"| samtools depth -a -b {region_bed} -"
    )

    log.debug("Depth pipeline: %s", depth_cmd)

    proc = subprocess.Popen(
        depth_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    callable_sequenced_bases = 0   # Σ depth at positions with depth ≥ min_depth
    total_sequenced_bases    = 0   # Σ depth at all positions
    total_positions          = 0   # number of genomic loci in region_bed

    for line in proc.stdout:
        parts = line.rstrip().split("\t")
        if len(parts) < 3:
            continue
        depth = int(parts[2])
        total_positions       += 1
        total_sequenced_bases += depth          # accumulate read-base count
        if depth >= min_depth:
            callable_sequenced_bases += depth   # only callable positions count

    proc.wait()
    if proc.returncode != 0:
        stderr_msg = proc.stderr.read()
        raise RuntimeError(
            f"samtools depth pipeline failed (exit {proc.returncode}):\n{stderr_msg}"
        )

    mean_depth = (
        total_sequenced_bases / total_positions if total_positions > 0 else 0.0
    )

    return {
        "callable_sequenced_bases": callable_sequenced_bases,
        "total_sequenced_bases":    total_sequenced_bases,
        "mean_depth":               round(mean_depth, 4),
        "total_positions":          total_positions,
    }


# ─── per-sample processing ────────────────────────────────────────────────────

# Column names appended to the output CSV, in order.
NEW_COLS = [
    "raw_total_bases",                 # overall sequenced bases in the input BAM
    "total_wes_bp",                    # genomic positions targeted (reference only)
    "total_sequenced_bases_wes",       # Σ depth — all filtered reads over WES BED
    "callable_sequenced_bases_wes",    # Σ depth at positions ≥ min_depth  ← denominator
    "pct_callable_bases_wes",          # callable / total × 100
    "mean_depth_wes",
    "total_padded_bp",
    "total_sequenced_bases_padded",
    "callable_sequenced_bases_padded",
    "pct_callable_bases_padded",
    "mean_depth_padded",
]

ERROR_ROW = {col: "ERROR" for col in NEW_COLS}


def process_sample(
    row: dict,
    bam_col: str,
    clean_wes_bed: str,
    clean_padded_bed: str,
    min_mapq: int,
    min_asxs: int,
    min_depth: int,
    threads: int,
) -> dict:
    """
    Compute callable sequenced base statistics for one sample.
    Returns the original CSV row with new columns appended.
    """
    sample_id = row.get("sample_id", "UNKNOWN")
    bam = row[bam_col].strip()

    if not os.path.isfile(bam):
        raise FileNotFoundError(f"BAM not found: {bam}")

    # Index the BAM if no .bai is present.
    if not (
        os.path.isfile(bam + ".bai")
        or os.path.isfile(bam.replace(".bam", ".bai"))
    ):
        log.warning("  No .bai index found — indexing %s …", bam)
        run(f"samtools index -@ {threads} {bam}")

    log.info("  [%s] Computing raw total bases in BAM …", sample_id)
    raw_total_bases = get_raw_total_bases(bam, threads)

    log.info("  [%s] Computing depth over WES BED …", sample_id)
    wes = compute_depth_stats(
        bam=bam,
        region_bed=clean_wes_bed,
        min_mapq=min_mapq,
        min_asxs=min_asxs,
        min_depth=min_depth,
        threads=threads,
    )

    log.info("  [%s] Computing depth over padded BED …", sample_id)
    pad = compute_depth_stats(
        bam=bam,
        region_bed=clean_padded_bed,
        min_mapq=min_mapq,
        min_asxs=min_asxs,
        min_depth=min_depth,
        threads=threads,
    )

    wes_total = wes["total_sequenced_bases"]
    pad_total = pad["total_sequenced_bases"]

    new_fields = {
        "raw_total_bases":                  raw_total_bases,
        "total_wes_bp":                     wes["total_positions"],
        "total_sequenced_bases_wes":        wes_total,
        "callable_sequenced_bases_wes":     wes["callable_sequenced_bases"],
        "pct_callable_bases_wes":           (
            round(100 * wes["callable_sequenced_bases"] / wes_total, 4)
            if wes_total else 0
        ),
        "mean_depth_wes":                   wes["mean_depth"],
        "total_padded_bp":                  pad["total_positions"],
        "total_sequenced_bases_padded":     pad_total,
        "callable_sequenced_bases_padded":  pad["callable_sequenced_bases"],
        "pct_callable_bases_padded":        (
            round(100 * pad["callable_sequenced_bases"] / pad_total, 4)
            if pad_total else 0
        ),
        "mean_depth_padded":                pad["mean_depth"],
    }

    return {**row, **new_fields}


# ─── argument parsing ─────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input_csv", required=True,
        help="Input CSV with sample metadata and a BAM path column.",
    )
    p.add_argument(
        "--bam_col", default="bam_path",
        help="Name of the column holding BAM file paths (default: bam_path).",
    )
    p.add_argument(
        "--wes_bed", required=True,
        help="Strict WES capture panel BED file.",
    )
    p.add_argument(
        "--blacklist", required=True,
        help="BED file of blacklisted regions to subtract from all targets.",
    )
    p.add_argument(
        "--genome_fai", default=None,
        help=(
            "Reference genome .fai index — required when --padded_bed is not "
            "supplied, so bedtools slop can clamp intervals to chromosome ends."
        ),
    )
    p.add_argument(
        "--padded_bed", default=None,
        help=(
            "Pre-computed padded BED (WES ± padding).  "
            "Auto-generated from --wes_bed if omitted."
        ),
    )
    p.add_argument(
        "--padding", type=int, default=100,
        help="Base-pairs to add on each side of WES intervals (default: 100).",
    )
    p.add_argument(
        "--min_mapq", type=int, default=40,
        help="Minimum mapping quality; reads below this are excluded (default: 40).",
    )
    p.add_argument(
        "--min_asxs", type=int, default=50,
        help=(
            "Minimum AS−XS score delta; reads below this are excluded (default: 50).  "
            "Reads lacking an XS tag are treated as XS = 0 and always pass."
        ),
    )
    p.add_argument(
        "--min_depth", type=int, default=1,
        help=(
            "Minimum DCS depth for a position to be callable (default: 1).  "
            "Raise to 2 for strict duplex-strand callability."
        ),
    )
    p.add_argument(
        "--threads", type=int, default=4,
        help="Number of samtools threads (default: 4).",
    )
    p.add_argument(
        "--output_csv", default="callable_bases_results.csv",
        help="Output CSV path (default: callable_bases_results.csv).",
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    return p.parse_args()


# ─── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    for tool in ("samtools", "bedtools"):
        check_tool(tool)

    ver = run("samtools --version", capture=True).stdout.splitlines()[0]
    log.info("Using %s", ver)

    # ── read input CSV ────────────────────────────────────────────────────────
    with open(args.input_csv) as fh:
        samples = list(csv.DictReader(fh))

    if not samples:
        log.error("Input CSV is empty.")
        sys.exit(1)

    if args.bam_col not in samples[0]:
        log.error(
            "Column '%s' not found in CSV.  Available columns: %s",
            args.bam_col,
            list(samples[0].keys()),
        )
        sys.exit(1)

    original_cols = list(samples[0].keys())

    with tempfile.TemporaryDirectory(prefix="callable_bases_") as tmp_dir:

        # ── prepare BEDs once; shared across all samples ──────────────────────

        clean_wes = os.path.join(tmp_dir, "wes_clean.bed")
        log.info("Preparing WES BED (blacklist subtraction + merge) …")
        clean_bed(args.wes_bed, args.blacklist, clean_wes)
        log.info("  WES BED size after cleaning: %d bp", bed_size(clean_wes))

        clean_pad = os.path.join(tmp_dir, "padded_clean.bed")
        if args.padded_bed:
            log.info("Preparing supplied padded BED (blacklist subtraction + merge) …")
            clean_bed(args.padded_bed, args.blacklist, clean_pad)
        else:
            if not args.genome_fai:
                log.error(
                    "--genome_fai is required when --padded_bed is not supplied."
                )
                sys.exit(1)
            log.info(
                "Generating padded BED (WES ± %d bp, blacklist subtracted) …",
                args.padding,
            )
            make_padded_bed(
                args.wes_bed,
                args.genome_fai,
                args.padding,
                args.blacklist,
                clean_pad,
            )
        log.info("  Padded BED size after cleaning: %d bp", bed_size(clean_pad))

        # ── process each sample ───────────────────────────────────────────────

        results = []
        for i, row in enumerate(samples, 1):
            sid = row.get("sample_id", f"sample_{i}")
            log.info("─── [%d/%d] %s ───", i, len(samples), sid)
            try:
                result = process_sample(
                    row=row,
                    bam_col=args.bam_col,
                    clean_wes_bed=clean_wes,
                    clean_padded_bed=clean_pad,
                    min_mapq=args.min_mapq,
                    min_asxs=args.min_asxs,
                    min_depth=args.min_depth,
                    threads=args.threads,
                )
                log.info(
                    "  callable_bases_wes=%d (%.1f%%)  "
                    "callable_bases_padded=%d (%.1f%%)  "
                    "mean_depth_wes=%.2f  mean_depth_padded=%.2f",
                    result["callable_sequenced_bases_wes"],
                    result["pct_callable_bases_wes"],
                    result["callable_sequenced_bases_padded"],
                    result["pct_callable_bases_padded"],
                    result["mean_depth_wes"],
                    result["mean_depth_padded"],
                )
            except Exception as exc:
                log.error("  FAILED for %s: %s", sid, exc)
                result = {**row, **ERROR_ROW}

            results.append(result)

    # ── write output ──────────────────────────────────────────────────────────
    if not results:
        log.error("No results to write.")
        sys.exit(1)

    fieldnames = original_cols + NEW_COLS
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    log.info("Done. Results written → %s", args.output_csv)


if __name__ == "__main__":
    main()