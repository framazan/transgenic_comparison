#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
ROOT_DIR="./"
GENOME_DIR="$ROOT_DIR/genomes"
STD_DIR="$ROOT_DIR/standardized_results"
OUT_DIR="$ROOT_DIR/protein_results"
LOG_DIR="$ROOT_DIR/protein_results/logs"
PYTHON_CLEANER="$ROOT_DIR/match_headers.py"  # <--- Make sure this points to the python script above
MAX_PARALLEL=${MAX_PARALLEL:-6}  # Maximum number of parallel tmux sessions
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
      echo "  --transgenic-only  Only process transgenic files"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Ensure output directories exist
mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- 1. Auto-Detect PASA ---
PASA_MAIN=$(which Launch_PASA_pipeline.pl)
if [[ -z "$PASA_MAIN" ]]; then echo "Error: PASA not found. Activate your conda env."; exit 1; fi
PASA_SCRIPT=$(find "$(dirname "$(dirname "$PASA_MAIN")")" -name "gff3_file_to_proteins.pl" | head -n 1)

# Check for tmux
if ! command -v tmux &> /dev/null; then
    echo "Error: tmux not found. Please install tmux."
    exit 1
fi

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

# Session prefix for tmux
SESSION_PREFIX="pasa_job_$$"
declare -a PASA_JOBS=()

# --- PHASE 1: Clean all files sequentially and prepare jobs ---
echo "==========================================="
echo "PHASE 1: Preparing files (sequential)"
echo "==========================================="

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
  
  # Filter for transgenic-only if flag is set
  if [[ "$TRANSGENIC_ONLY" == "true" ]] && [[ ! "$tool" =~ transgenic ]]; then
      echo "[SKIP] $base (not transgenic)"
      continue
  fi      

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
  clean_genome_path="$OUT_DIR/.temp_clean_genome_${species}_${tool}${suffix}.fa"
  tmp_prot="$OUT_DIR/.tmp_prot_${species}_${tool}${suffix}.fasta"

  echo "---------------------------------------------------"
  echo "[INFO] Cleaning headers for: $base"

  # --- STEP A: Run Python Cleaner (sequential) ---
  python3 "$PYTHON_CLEANER" "$gff" "$raw_genome_path" "$clean_genome_path"

  # Store job info for parallel execution
  PASA_JOBS+=("$gff|$clean_genome_path|$tmp_prot|$out_fa|$base|$genome_filename")

done

# --- PHASE 2: Run PASA in parallel using tmux ---
echo ""
echo "==========================================="
echo "PHASE 2: Running PASA jobs in parallel (max $MAX_PARALLEL)"
echo "==========================================="

declare -a ACTIVE_SESSIONS=()
job_index=0

for job_info in "${PASA_JOBS[@]}"; do
  IFS='|' read -r gff clean_genome_path tmp_prot out_fa base genome_filename <<< "$job_info"
  
  session_name="${SESSION_PREFIX}_${job_index}"
  
  # Wait if we've reached max parallel jobs
  while [[ ${#ACTIVE_SESSIONS[@]} -ge $MAX_PARALLEL ]]; do
    sleep 2
    # Check which sessions are still running
    declare -a STILL_RUNNING=()
    for sess in "${ACTIVE_SESSIONS[@]}"; do
      if tmux has-session -t "$sess" 2>/dev/null; then
        STILL_RUNNING+=("$sess")
      fi
    done
    ACTIVE_SESSIONS=("${STILL_RUNNING[@]}")
  done
  
  echo "[INFO] Launching PASA for: $base (session: $session_name)"
  
  # Create a wrapper script for this job
  wrapper_script="$OUT_DIR/.pasa_wrapper_${job_index}.sh"
  log_file="$LOG_DIR/pasa_${base%.gff3}.log"
  cat > "$wrapper_script" << EOF
#!/usr/bin/env bash
set -euo pipefail

gff="$gff"
clean_genome_path="$clean_genome_path"
tmp_prot="$tmp_prot"
out_fa="$out_fa"
base="$base"
genome_filename="$genome_filename"
PASA_SCRIPT="$PASA_SCRIPT"
log_file="$log_file"

# Redirect stderr to log file to prevent contamination of protein output
if perl "\$PASA_SCRIPT" --gff3 "\$gff" --fasta "\$clean_genome_path" > "\$tmp_prot" 2>"\$log_file"; then
    if [[ ! -s "\$tmp_prot" ]]; then
        echo "[WARN] Generated protein file is empty for \$base"
    fi
    {
      echo "# genome_original=\$genome_filename"
      echo "# gff3_file=\$base"
      cat "\$tmp_prot"
    } > "\$out_fa"
    echo "[SUCCESS] Generated: \$(basename "\$out_fa")"
    rm -f "\$tmp_prot" "\$clean_genome_path"
    touch "$OUT_DIR/.done_${job_index}"
else
    echo "[ERROR] PASA failed for \$base (see \$log_file)"
    touch "$OUT_DIR/.failed_${job_index}"
fi
EOF
  chmod +x "$wrapper_script"
  
  # Launch in tmux (detached)
  tmux new-session -d -s "$session_name" "bash '$wrapper_script'; exit"
  
  ACTIVE_SESSIONS+=("$session_name")
  ((job_index++))
done

# --- PHASE 3: Wait for all jobs to complete ---
echo ""
echo "==========================================="
echo "PHASE 3: Waiting for all PASA jobs to complete..."
echo "==========================================="

while [[ ${#ACTIVE_SESSIONS[@]} -gt 0 ]]; do
  sleep 5
  declare -a STILL_RUNNING=()
  for sess in "${ACTIVE_SESSIONS[@]}"; do
    if tmux has-session -t "$sess" 2>/dev/null; then
      STILL_RUNNING+=("$sess")
    fi
  done
  ACTIVE_SESSIONS=("${STILL_RUNNING[@]}")
  echo "[INFO] ${#ACTIVE_SESSIONS[@]} jobs still running..."
done

# --- PHASE 4: Report results ---
echo ""
echo "==========================================="
echo "PHASE 4: Summary"
echo "==========================================="

success_count=$(find "$OUT_DIR" -name ".done_*" 2>/dev/null | wc -l | tr -d ' ')
fail_count=$(find "$OUT_DIR" -name ".failed_*" 2>/dev/null | wc -l | tr -d ' ')

echo "Completed: $success_count"
echo "Failed: $fail_count"

# Cleanup marker and wrapper files
rm -f "$OUT_DIR"/.done_* "$OUT_DIR"/.failed_* "$OUT_DIR"/.pasa_wrapper_*.sh

echo ""
echo "All done!"