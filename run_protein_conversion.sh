#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
ROOT_DIR="/Users/filipr/Desktop_comparison"
GENOME_DIR="$ROOT_DIR/genomes"
STD_DIR="$ROOT_DIR/standardized_results"
OUT_DIR="$ROOT_DIR/protein_results"
PYTHON_CLEANER="$ROOT_DIR/match_headers.py"  # <--- Make sure this points to the python script above

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# --- 1. Auto-Detect PASA ---
# (Same detection logic as before)
PASA_MAIN=$(which Launch_PASA_pipeline.pl)
if [[ -z "$PASA_MAIN" ]]; then echo "Error: PASA not found. Activate your conda env."; exit 1; fi
PASA_SCRIPT=$(find "$(dirname "$(dirname "$PASA_MAIN")")" -name "gff3_file_to_proteins.pl" | head -n 1)

# --- 2. Helper Function (Species -> Genome) ---
get_genome_fasta() {
  local species="$1"
  local type="${2:-standard}"

  if [[ "$type" == "softmasked" ]]; then
    case "$species" in
      A_thaliana)    echo "softmasked/Athaliana_167_TAIR10.softmasked.fa" ;;
      B_distachyon)  echo "softmasked/Bdistachyon_314_v3.0.softmasked.fa" ;;
      B_rapa)        echo "softmasked/BrapaO_302V_711_v1.0.softmasked.fa" ;;
      Gmax|G_max)    echo "softmasked/Gmax_880_v6.0.softmasked.fa" ;;
      Lsativa|L_sativa) echo "softmasked/Lsativa_467_v8.softmasked.fa" ;;
      O_sativa)      echo "softmasked/Osativa_323_v7.0.softmasked.fa" ;;
      P_patens)      echo "softmasked/Ppatens_318_v3.softmasked.fa" ;;
      P_trichocarpa) echo "softmasked/Ptrichocarpa_533_v4.0.softmasked.fa" ;;
      S_bicolor)     echo "softmasked/Sbicolor_730_v5.0.softmasked.fa" ;;
      S_italica)     echo "softmasked/Sitalica_312_v2.softmasked.fa" ;;
      Slycopersicum|S_lycopersicum) echo "softmasked/Slycopersicum_796_SL5.0.softmasked.fa" ;;
      Z_mays)        echo "softmasked/Zmays_493_APGv4.softmasked.fa" ;;
      *)             echo "" ;;
    esac
  else
    case "$species" in
      A_thaliana)    echo "Athaliana_167_TAIR10.fa" ;;
      B_distachyon)  echo "Bdistachyon_314_v3.0.fa" ;;
      B_rapa)        echo "BrapaO_302V_711_v1.0.fa" ;;
      Gmax|G_max)    echo "Gmax_880_v6.0.fa" ;;
      Lsativa|L_sativa) echo "Lsativa_467_v8.fa" ;;
      O_sativa)      echo "Osativa_323_v7.0.fa" ;;
      P_patens)      echo "Ppatens_318_v3.fa" ;;
      P_trichocarpa) echo "Ptrichocarpa_533_v4.0.fa" ;;
      S_bicolor)     echo "Sbicolor_730_v5.0.fa" ;;
      S_italica)     echo "Sitalica_312_v2.fa" ;;
      Slycopersicum|S_lycopersicum) echo "Slycopersicum_796_ITAG5.0.fa" ;;
      V_vinifera)    echo "Vvinifera_T2T_ref.fa" ;;
      Z_mays)        echo "Zmays_493_APGv4.fa" ;;
      *)             echo "" ;;
    esac
  fi
}

# --- 3. Main Execution ---
cd "$ROOT_DIR"
shopt -s nullglob

for gff in "$STD_DIR"/*.gff3; do
  
  base="$(basename "$gff")"           
  name_no_ext="${base%.gff3}"
  
  # Check for softmasked
  if [[ "$name_no_ext" == *"_softmasked" ]]; then
      is_softmasked="true"
      name_no_ext="${name_no_ext%_softmasked}"
  else
      is_softmasked="false"
  fi
  
  species="${name_no_ext%_*}"        
  tool="${name_no_ext#${species}_}"      

  if [[ "$is_softmasked" == "true" ]]; then
      genome_filename="$(get_genome_fasta "$species" "softmasked")"
  else
      genome_filename="$(get_genome_fasta "$species")"
  fi

  if [[ -z "$genome_filename" ]]; then
    echo "[WARN] No genome map for $species, skipping..."
    continue
  fi

  raw_genome_path="$GENOME_DIR/$genome_filename"
  
  suffix=""
  if [[ "$is_softmasked" == "true" ]]; then suffix="_softmasked"; fi

  # Define paths for outputs and temp files
  out_fa="$OUT_DIR/${species}_${tool}${suffix}.prot.fasta"
  clean_genome_path="$OUT_DIR/.temp_clean_genome_${species}_${tool}${suffix}.fa" # Temp location
  
  tmp_prot="$OUT_DIR/.tmp_prot_${species}_${tool}${suffix}.fasta"

  echo "---------------------------------------------------"
  echo "[INFO] Processing: $base"

  # --- STEP A: Run Python Cleaner ---
  # This creates a temporary genome file where headers match the GFF3 exactly
  python3 "$PYTHON_CLEANER" "$gff" "$raw_genome_path" "$clean_genome_path"

  # --- STEP B: Run PASA ---
  # We use the $clean_genome_path instead of the raw genome
  if ! perl "$PASA_SCRIPT" --gff3 "$gff" --fasta "$clean_genome_path" > "$tmp_prot"; then
      echo "[ERROR] PASA failed for $base even after header cleaning."
      # If it fails, we keep the temp files for debugging
  else
      # --- STEP C: Finalize ---
      # Add metadata header and save to final location
      if [[ ! -s "$tmp_prot" ]]; then
         echo "[WARN] Generated protein file is empty (excluding headers) for $base"
      fi

      {
        echo "# genome_original=$genome_filename"
        echo "# gff3_file=$base"
        cat "$tmp_prot"
      } > "$out_fa"
      
      echo "[SUCCESS] Generated: $(basename "$out_fa")"
      
      # Cleanup temp files only on success
      rm -f "$tmp_prot" "$clean_genome_path"
  fi

done