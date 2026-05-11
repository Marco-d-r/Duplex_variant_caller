#!/bin/bash
set -e

# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

# 1. Path to gnomAD
GNOMAD_VCF="/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/reference_genome/af-only-gnomad.hg38_common_variants.vcf.gz"

# 2. Path to CSV
PON_CSV="pon_list.csv"

# 3. Output directory
OUTDIR="../Sample_files"

# 4. Input Directory (No wildcards here, just the folder path)
INPUT_DIR="/Users/marcodupuisrodriguez/Documents/PhD/VCS/Marco_sequencing_data_and_variant_processing/Extended_variant_extracting/final_results/UDI56"

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
PON_CSV_ABS=$(get_abs_path "$PON_CSV")
OUTDIR_ABS=$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)
WORK_DIR="$OUTDIR_ABS/temp_work"

mkdir -p "$WORK_DIR"
echo "--- Starting Annotation Pipeline (Fixed Header) ---"

check_index "$GNOMAD_VCF_ABS"

# ------------------------------------------------------------------------------
# STEP 1: PREPARE PON
# ------------------------------------------------------------------------------
MERGED_PON="$WORK_DIR/merged_pon.vcf.gz"
if [ ! -f "$MERGED_PON" ]; then
    echo "[1/4] Merging PoN files..."
    if [ ! -f "$PON_CSV_ABS" ]; then echo "Error: CSV file not found"; exit 1; fi
    grep -v '^$' "$PON_CSV_ABS" | cut -d',' -f1 > "$WORK_DIR/pon_file_list.txt"
    bcftools merge --file-list "$WORK_DIR/pon_file_list.txt" --force-samples | bcftools view -G -Oz -o "$MERGED_PON"
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

for VCF_IN in "$INPUT_DIR"/*final_combined.vcf.gz; do
    
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
    
    check_index "$VCF_IN_ABS"

    # --------------------------------------------------------------------------
    # STEP 2: GNOMAD CHECK
    # --------------------------------------------------------------------------
    echo "      Checking against gnomAD..."
    DIR_ISEC="$WORK_DIR/isec_gnomad_${SAFE_NAME}"
    mkdir -p "$DIR_ISEC"
    
    # Run isec (generates 0000.vcf [private] and 0002.vcf [shared])
    bcftools isec -p "$DIR_ISEC" -n~10 "$VCF_IN_ABS" "$GNOMAD_VCF_ABS"

    if [ ! -f "$DIR_ISEC/0000.vcf" ]; then echo "Error: isec failed."; exit 1; fi

    # 1. Annotate NO (Private)
    #    We pipe awk -> annotate to fix the header immediately
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
        # Create empty dummy with correct header if no shared variants
        bcftools view -h "$VCF_IN_ABS" \
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

    # Cleanup
    rm -rf "$DIR_ISEC" "$DIR_ISEC_PON" "$WORK_DIR"/step1*

done

echo "--- All processing complete ---"