#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
export ROOT_DIR="./"
export INPUT_DIR="$ROOT_DIR/protein_results"
export OUTPUT_DIR="$ROOT_DIR/busco_results"
export LINEAGE="viridiplantae_odb10"

# --- 1. Choose your Weapon (Conda vs Docker) ---
# Set to "conda" if you fixed the install. Set to "docker" if Conda is still broken.
MODE="conda"  

mkdir -p "$OUTPUT_DIR"

# --- 2. Setup Function ---
setup_busco() {
    if [[ "$MODE" == "conda" ]]; then
        # Try to find BUSCO
        if ! command -v busco &> /dev/null; then
            echo "[ERROR] BUSCO not found. Did you 'mamba activate busco_env'?"
            exit 1
        fi
        # Trigger download of lineage in serial (prevents parallel race conditions)
        busco --list-datasets > /dev/null
    elif [[ "$MODE" == "docker" ]]; then
        # Pull the image once before starting
        echo "[INFO] Checking Docker image..."
        docker pull ezlabgva/busco:v5.7.1_cv1
    fi
}

# --- 3. The Worker Function ---
run_sample() {
    local fasta="$1"
    local base=$(basename "$fasta" .prot.fasta) # Extract "A_thaliana_Tiberius"
    local out_path="$OUTPUT_DIR/$base"

    # Skip if done
    if [[ -d "$out_path" && -f "$out_path/short_summary.json" ]]; then
        echo "[SKIP] $base already finished."
        return
    fi

    echo "[START] $base ($MODE mode)..."

    if [[ "$MODE" == "conda" ]]; then
        # --- METHOD A: CONDA ---
        busco \
            -i "$fasta" \
            -o "$base" \
            --out_path "$OUTPUT_DIR" \
            -m prot \
            -l viridiplantae_odb10 \
            --download_path "$ROOT_DIR/busco_downloads" \
            -c 2 \
            --quiet --force --offline

    elif [[ "$MODE" == "docker" ]]; then
        # --- METHOD B: DOCKER ---
        # We mount ROOT_DIR to /data inside container
        # Paths must be relative to /data
        local rel_input="/data/protein_results/$(basename "$fasta")"
        local rel_out_path="/data/busco_results"
        
        docker run --rm \
            -v "$ROOT_DIR":/data \
            -w /data \
            --user $(id -u):$(id -g) \
            ezlabgva/busco:v5.7.1_cv1 \
            busco \
            -i "$rel_input" \
            -o "$base" \
            --out_path "$rel_out_path" \
            -m prot \
            -l "$LINEAGE" \
            -c 2 \
            --quiet --force --offline
    fi

    echo "[DONE] $base"
}

export -f run_sample
export MODE
export ROOT_DIR
export OUTPUT_DIR
export LINEAGE

# --- 4. Execution ---
setup_busco

echo "[INFO] Found $(find "$INPUT_DIR" -name "*.prot.fasta" | wc -l) protein files."
echo "[INFO] Starting Parallel BUSCO (4 jobs at once)..."

# Use GNU Parallel to run 4 jobs, 2 CPUs each = 8 cores total usage
find "$INPUT_DIR" -name "*.prot.fasta" | parallel -j 4 --eta run_sample {}