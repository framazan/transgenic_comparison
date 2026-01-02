#!/bin/bash
#SBATCH --job-name=masking
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --output=mask_%j.log
#SBATCH --error=mask_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32     # ADJUST THIS to your cluster limits
#SBATCH --mem=64G             # ADJUST THIS to your cluster limits
#SBATCH --time=48:00:00        # Time limit

# --- CONFIGURATION ---
# UPDATE THESE PATHS - you can add more genomes here, and then add a line at the end of the file (scroll down there)
GENOME1="/data/gpfs/assoc/pgl/data/Transgenic/genomes/Athaliana_167_TAIR10.fa"
GENOME2="/data/gpfs/assoc/pgl/data/Transgenic/genomes/Vvinifera_T2T_ref.fa"
OUTDIR="./softmasked"

# Path to the container and tools (Assumes they are in current directory)
CONTAINER="./tetools.sif"
BIN_DIR="$PWD/bin"

# --- SETUP & CHECKS ---
THREADS=$SLURM_CPUS_PER_TASK

# Ensure we have absolute paths for tools to avoid "not found" errors
CONTAINER=$(realpath "$CONTAINER")
SPLIT_SCRIPT=$(realpath "$BIN_DIR/splitMfasta.pl")
PARSE_SCRIPT=$(realpath "$BIN_DIR/parseTrfOutput.py")
BEDTOOLS=$(realpath "$BIN_DIR/bedtools")

# Add bin to path for this session
export PATH=$BIN_DIR:$PATH

if [ ! -f "$CONTAINER" ]; then
    echo "Error: tetools.sif not found at $CONTAINER"
    exit 1
fi

# Load GNU Parallel (Most HPCs use the module system)
module load parallel 2>/dev/null || echo "Module parallel load skipped"

# Helper function to run commands INSIDE the container
# We bind (-B) the current working directory tree so the container sees files
sing_exec() {
    singularity exec -B "$PWD":"$PWD" -B "$OUTDIR":"$OUTDIR" "$CONTAINER" "$@"
}

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
TEMP_BASE="$OUTDIR/temp_work_files"
mkdir -p "$TEMP_BASE"

# --- MAIN PROCESS ---
process_genome() {
    local INPUT_FASTA=$1
    INPUT_FASTA=$(realpath "$INPUT_FASTA")
    local BASENAME=$(basename "$INPUT_FASTA" .fasta)
    BASENAME=$(basename "$BASENAME" .fa)
    
    local FINAL_FASTA="$OUTDIR/$BASENAME.final_softmasked.fasta"
    local FINAL_GFF="$OUTDIR/$BASENAME.trf.gff3"

    echo "=========================================="
    echo "Processing $BASENAME"
    echo "=========================================="

    if [ -f "$FINAL_FASTA" ]; then
        echo "Output exists. Skipping."
        return
    fi

    # Create temp dir
    WORKDIR=$(mktemp -d "${TEMP_BASE}/${BASENAME}_working_XXXXXX")
    
    # Copy genome to working dir
    cp "$INPUT_FASTA" "$WORKDIR/$BASENAME.fa"
    cd "$WORKDIR" || exit 1

    # 1. RepeatModeler (Inside Container)
    echo "[$BASENAME] Building Database..."
    sing_exec BuildDatabase -name "$BASENAME" "$BASENAME.fa" > build_db.log 2>&1

    echo "[$BASENAME] Running RepeatModeler (LTRStruct enabled)..."
    sing_exec RepeatModeler -database "$BASENAME" -threads "$THREADS" -LTRStruct > repmod.log 2>&1

    # Find library
    if [ -f "$BASENAME-families.fa" ]; then
        LIB_FILE="$BASENAME-families.fa"
    else
        echo "Error: RepeatModeler failed. Check logs in $WORKDIR"
        exit 1
    fi

    # 2. RepeatMasker (Inside Container)
    echo "[$BASENAME] Running RepeatMasker..."
    sing_exec RepeatMasker -pa "$THREADS" -lib "$LIB_FILE" -xsmall "$BASENAME.fa" > repmask.log 2>&1
    
    if [ ! -f "$BASENAME.fa.masked" ]; then
        echo "Error: RepeatMasker failed."
        exit 1
    fi

    # 3. TRF Refinement
    echo "[$BASENAME] Running TRF refinement..."
    mkdir -p trf
    cd trf || exit 1
    ln -s ../"$BASENAME.fa.masked" genome.fa

    # FIX: Run perl script INSIDE container, passing the absolute path to the script
    echo "Splitting Fasta..."
    sing_exec perl "$SPLIT_SCRIPT" --minsize=25000000 genome.fa

    # Parallel Steps
    # We use 'parallel' to spawn many singularity jobs
    # Note: We use the container for TRF and parsing, but local static bedtools for merging
    
    echo "Running TRF..."
    ls genome.split.*.fa | parallel "singularity exec -B $PWD:$PWD -B $OUTDIR:$OUTDIR $CONTAINER trf {} 2 7 7 80 10 50 500 -d -m -h > {}.log 2>&1"

    echo "Parsing TRF..."
    ls genome.split.*.fa.2.7.7.80.10.50.500.dat | parallel "singularity exec -B $PWD:$PWD -B $OUTDIR:$OUTDIR $CONTAINER python3 $PARSE_SCRIPT {} --minCopies 1 --statistics {}.STATS > {}.raw.gff 2> {}.parsedLog"

    echo "Sorting GFFs..."
    ls genome.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff | parallel "sort -k1,1 -k4,4n -k5,5n {} > {}.sorted 2> {}.sortLog"

    echo "Merging GFFs..."
    # Using local static bedtools
    find . -name "*.sorted" | while read f; do
         "$BEDTOOLS" merge -i "$f" | awk 'BEGIN{OFS="\t"} {print $1,"trf","repeat",$2+1,$3,".",".",".","."}' > "$f.merged.gff" 2> "$f.bedtools_merge.log"
    done

    echo "Final Masking..."
    # Using local static bedtools
    ls genome.split.*.fa | parallel "$BEDTOOLS maskfasta -fi {} -bed {}.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff -fo {}.combined.masked -soft > {}.bedtools_mask.log 2>&1"

    # 4. Consolidate
    echo "[$BASENAME] Saving..."
    cat genome.split.*.fa.combined.masked > "$FINAL_FASTA"
    
    if ls *.merged.gff 1> /dev/null 2>&1; then
        cat *.merged.gff > "$FINAL_GFF"
    else
        touch "$FINAL_GFF"
    fi

    # 5. Cleanup
    cd "$OUTDIR"
    rm -rf "$WORKDIR"
    echo "[$BASENAME] Done."
}

# Run
process_genome "$GENOME1"
process_genome "$GENOME2"
# copy above line if you have added more genomes...