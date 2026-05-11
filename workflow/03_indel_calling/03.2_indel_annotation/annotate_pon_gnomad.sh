#!/bin/bash
set -e

# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

# 1. Path to gnomAD (Original)
GNOMAD_VCF="/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/reference_genome/af-only-gnomad.hg38_common_variants.vcf.gz"

# 2. Path to Reference Genome (REQUIRED for Indel Normalization)
#    Updated to match the file you listed in your directory
REF_FASTA="/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/reference_genome/GRCh38.primary_assembly.genome.fa"

# 3. Path to CSV
PON_CSV="pon_list.csv"

# 4. Output directory
OUTDIR="../Sample_files"

# 5. Input Directory
INPUT_DIR="/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/variant_caller/indel_callers/final_results/UDI54/"

# ==============================================================================
# HELPER FUNCTIONS & SETUP
# ==============================================================================

get_abs_path() {
    if [ -f "$1" ]; then echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; else echo "$1"; fi
}

check_index() {
    local vcf="$1"
    if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.csi" ]; then
        echo "      [!] Index missing for $(basename "$vcf"). Indexing now..."
        bcftools index -t "$vcf"
    fi
}

GNOMAD_VCF_ABS=$(get_abs_path "$GNOMAD_VCF")
REF_FASTA_ABS=$(get_abs_path "$REF_FASTA")
PON_CSV_ABS=$(get_abs_path "$PON_CSV")
OUTDIR_ABS=$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)
WORK_DIR="$OUTDIR_ABS/temp_work"

# Validate Reference Fasta
if [ ! -f "$REF_FASTA_ABS" ]; then
    echo "CRITICAL ERROR: Reference Fasta not found at $REF_FASTA_ABS"
    echo "Please update the REF_FASTA path in the configuration section."
    exit 1
fi

mkdir -p "$WORK_DIR"
echo "--- Starting Annotation Pipeline (Indel Optimized + Canonical Fix) ---"

# ------------------------------------------------------------------------------
# STEP 0: PREPARE GNOMAD (FIXED)
# ------------------------------------------------------------------------------
GNOMAD_SPLIT="$WORK_DIR/gnomad_split_norm.vcf.gz"

# We explicitly list canonical chromosomes to filter out "random" contigs
# that cause normalization crashes.
CANONICAL_REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

if [ ! -f "$GNOMAD_SPLIT" ]; then
    echo "[0/4] Preparing gnomAD (Filtering & Normalizing)..."
    
    # 1. view --regions: Restrict to main chromosomes (fixes the crash)
    # 2. norm -m -: Split multiallelic sites (e.g. A -> T,G becomes two lines)
    # 3. norm -f: Left-align indels using the reference genome
    bcftools view --regions "$CANONICAL_REGIONS" "$GNOMAD_VCF_ABS" -Ou \
    | bcftools norm -m - -f "$REF_FASTA_ABS" -Oz -o "$GNOMAD_SPLIT"
    
    bcftools index -t "$GNOMAD_SPLIT"
else
    echo "[0/4] Using existing prepared gnomAD."
fi

# ------------------------------------------------------------------------------
# STEP 1: PREPARE PON
# ------------------------------------------------------------------------------
MERGED_PON="$WORK_DIR/merged_pon_split.vcf.gz"

if [ ! -f "$MERGED_PON" ]; then
    echo "[1/4] Merging and Normalizing PoN files..."
    if [ ! -f "$PON_CSV_ABS" ]; then echo "Error: CSV file not found"; exit 1; fi
    
    grep -v '^$' "$PON_CSV_ABS" | cut -d',' -f1 > "$WORK_DIR/pon_file_list.txt"
    
    # Merge -> Normalize (Split multiallelics & Left-align) -> Output
    bcftools merge --file-list "$WORK_DIR/pon_file_list.txt" --force-samples -Ou \
    | bcftools norm -m - -f "$REF_FASTA_ABS" \
    | bcftools view -G -Oz -o "$MERGED_PON"
    
    bcftools index -t "$MERGED_PON"
else
    echo "[1/4] Using existing PoN."
fi

# Create Header Definitions
cat <<EOF > "$WORK_DIR/new_headers.txt"
##INFO=<ID=IN_GNOMAD_COMMON,Number=1,Type=String,Description="Variant present in gnomAD common variants (Yes/No)">
##INFO=<ID=IN_CB_PON,Number=1,Type=String,Description="Variant present in Panel of Normals (Yes/No)">
EOF

# ==============================================================================
# MAIN LOOP
# ==============================================================================

for VCF_IN in "$INPUT_DIR"/*indel_candidates_annotated.vcf.gz; do
    
    if [ ! -f "$VCF_IN" ]; then
        echo "Error: No files matching *final_combined.vcf.gz found in $INPUT_DIR"
        exit 1
    fi

    VCF_IN_ABS=$(get_abs_path "$VCF_IN")
    BASENAME=$(basename "$VCF_IN_ABS")
    SAFE_NAME="${BASENAME//./_}"
    FINAL_VCF="$OUTDIR_ABS/$BASENAME"

    echo "----------------------------------------------------------------"
    echo "[2/4] Processing: $BASENAME"
    
    # Normalize Input (Safety check for Indels)
    VCF_IN_NORM="$WORK_DIR/${SAFE_NAME}_norm.vcf.gz"
    bcftools norm -m - -f "$REF_FASTA_ABS" "$VCF_IN_ABS" -Oz -o "$VCF_IN_NORM"
    bcftools index -t "$VCF_IN_NORM"

    # --------------------------------------------------------------------------
    # STEP 2: GNOMAD CHECK
    # --------------------------------------------------------------------------
    echo "      Checking against gnomAD..."
    DIR_ISEC="$WORK_DIR/isec_gnomad_${SAFE_NAME}"
    mkdir -p "$DIR_ISEC"
    
    # Compare Normalized Input vs Normalized gnomAD
    bcftools isec -p "$DIR_ISEC" -n~10 "$VCF_IN_NORM" "$GNOMAD_SPLIT"

    if [ ! -f "$DIR_ISEC/0000.vcf" ]; then echo "Error: isec failed."; exit 1; fi

    # 1. Annotate NO (Private)
    bcftools view "$DIR_ISEC/0000.vcf" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} { if($8==".") $8="IN_GNOMAD_COMMON=No"; else $8=$8";IN_GNOMAD_COMMON=No"; print }' \
    | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC/tagged_no.vcf.gz"
    bcftools index "$DIR_ISEC/tagged_no.vcf.gz"

    # 2. Annotate YES (Shared)
    if [ -f "$DIR_ISEC/0002.vcf" ]; then
        bcftools view "$DIR_ISEC/0002.vcf" \
        | awk 'BEGIN{OFS="\t"} /^#/ {print; next} { if($8==".") $8="IN_GNOMAD_COMMON=Yes"; else $8=$8";IN_GNOMAD_COMMON=Yes"; print }' \
        | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC/tagged_yes.vcf.gz"
    else
        bcftools view -h "$VCF_IN_NORM" \
        | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC/tagged_yes.vcf.gz"
    fi
    bcftools index "$DIR_ISEC/tagged_yes.vcf.gz"

    # Combine
    bcftools concat -a "$DIR_ISEC/tagged_no.vcf.gz" "$DIR_ISEC/tagged_yes.vcf.gz" \
    | bcftools sort -Oz -o "$WORK_DIR/step1.vcf.gz"
    bcftools index "$WORK_DIR/step1.vcf.gz"

    # --------------------------------------------------------------------------
    # STEP 3: PON CHECK
    # --------------------------------------------------------------------------
    echo "      Checking against PoN..."
    DIR_ISEC_PON="$WORK_DIR/isec_pon_${SAFE_NAME}"
    mkdir -p "$DIR_ISEC_PON"

    # Compare against Normalized PoN
    bcftools isec -p "$DIR_ISEC_PON" "$WORK_DIR/step1.vcf.gz" "$MERGED_PON"

    # 1. Annotate NO
    bcftools view "$DIR_ISEC_PON/0000.vcf" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} { if($8==".") $8="IN_CB_PON=No"; else $8=$8";IN_CB_PON=No"; print }' \
    | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC_PON/tagged_no.vcf.gz"
    bcftools index "$DIR_ISEC_PON/tagged_no.vcf.gz"

    # 2. Annotate YES
    if [ -f "$DIR_ISEC_PON/0002.vcf" ]; then
        bcftools view "$DIR_ISEC_PON/0002.vcf" \
        | awk 'BEGIN{OFS="\t"} /^#/ {print; next} { if($8==".") $8="IN_CB_PON=Yes"; else $8=$8";IN_CB_PON=Yes"; print }' \
        | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC_PON/tagged_yes.vcf.gz"
    else
        bcftools view -h "$WORK_DIR/step1.vcf.gz" \
        | bcftools annotate -h "$WORK_DIR/new_headers.txt" --no-version -Oz -o "$DIR_ISEC_PON/tagged_yes.vcf.gz"
    fi
    bcftools index "$DIR_ISEC_PON/tagged_yes.vcf.gz"

    # Final Combine
    bcftools concat -a "$DIR_ISEC_PON/tagged_no.vcf.gz" "$DIR_ISEC_PON/tagged_yes.vcf.gz" \
    | bcftools sort -Oz -o "$FINAL_VCF"
    
    bcftools index -t "$FINAL_VCF"
    echo "      Done! Saved to $FINAL_VCF"

    # Cleanup intermediate files
    rm -rf "$DIR_ISEC" "$DIR_ISEC_PON" "$WORK_DIR"/step1* "$VCF_IN_NORM"*

done

echo "--- All processing complete ---"