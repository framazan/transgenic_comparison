#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
ROOT_DIR="./"
QUERY_DIR="$ROOT_DIR/standardized_results"   
REF_DIR="$ROOT_DIR/reference_annotations"   
OUT_DIR="$ROOT_DIR/gffcompare_results"
TRANSGENIC_ONLY=false

# --- Parse command line arguments ---
while [[ $# -gt 0 ]]; do
  case $1 in
    --transgenic-only)
      TRANSGENIC_ONLY=true
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--transgenic-only]"
      echo "  --transgenic-only  Only process transgenic GFF files"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

mkdir -p "$OUT_DIR"

# --- 1. Define Known Species ---
# We use this list to scan filenames. Order matters slightly (longer matches first usually better, but these are distinct enough).
KNOWN_SPECIES="A_thaliana B_distachyon B_rapa G_max L_sativa O_sativa P_patens P_trichocarpa S_bicolor S_italica S_lycopersicum V_vinifera Z_mays"

# --- 2. Helper: Map Species to Reference Filename ---
get_ref_gff() {
  local species="$1"
  case "$species" in
    A_thaliana)    echo "Athaliana_167_TAIR10.gene.clean.gff3" ;;
    B_distachyon)  echo "Bdistachyon_314_v3.1.gene_exons.clean.gff3" ;;
    B_rapa)        echo "BrapaO_302V_711_v1.1.gene.gff3" ;;
    G_max)          echo "Gmax_880_Wm82.a6.v1.gene_exons.clean.gff3" ;;
    L_sativa)       echo "Lsativa_467_v5.gene_exons.gff3" ;;
    O_sativa)      echo "Osativa_323_v7.0.gene_exons.exon.gff3" ;;
    P_patens)      echo "Ppatens_318_v3.3.gene_exons.clean.gff3" ;;
    P_trichocarpa) echo "Ptrichocarpa_533_v4.1.gene_exons.clean.gff3" ;;
    S_bicolor)     echo "Sbicolor_730_v5.1.gene_exons.clean.gff3" ;;
    S_italica)     echo "Sitalica_312_v2.2.gene_exons.clean.gff3" ;;
    S_lycopersicum) echo "Slycopersicum_796_ITAG5.0.gene.gff3" ;;
    V_vinifera)    echo "Vvinifera_PN40024_5.1_on_T2T_ref.exon.gff3" ;;
    Z_mays)        echo "Zmays_493_RefGen_V4.gene_exons.exon.gff3" ;;
    *)             echo "" ;;
  esac
}

# --- 3. Main Loop ---
shopt -s nullglob
for query_gff in "$QUERY_DIR"/*.gff3; do
    
    filename=$(basename "$query_gff")
    
    # Skip non-transgenic files if --transgenic-only is set
    if [[ "$TRANSGENIC_ONLY" == true && ! "$filename" == *transgenic* ]]; then
        continue
    fi      
    base="${filename%.gff3}"               
    
    # --- FIX: Robust Species Detection ---
    # Instead of splitting by underscores, we loop through known species 
    # and see if the filename STARTS with one of them.
    species=""
    for s in $KNOWN_SPECIES; do
        if [[ "$base" == "$s"* ]]; then
            species="$s"
            break
        fi
    done

    # Handle unmatched files
    if [[ -z "$species" ]]; then
        echo "[SKIP] Could not detect species in filename: $filename"
        continue
    fi

    # Find Reference
    ref_filename=$(get_ref_gff "$species")
    ref_path="$REF_DIR/$ref_filename"

    if [[ -z "$ref_filename" ]]; then
        echo "[SKIP] No reference mapping configured for: $species"
        continue
    fi

    if [[ ! -f "$ref_path" ]]; then
        echo "[WARN] Reference file missing: $ref_path"
        continue
    fi

    echo "Comparing: $base vs $ref_filename"

    # Run gffcompare
    gffcompare -r "$ref_path" -o "$OUT_DIR/$base" "$query_gff"
    
done