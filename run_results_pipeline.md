# 🧬 Postprocessing Environment Setup and Analysis

This documentation guides you through how to set up mamba/conda environments, intall PASA and gffcompare, and run accuracy/precision analysis on annotation outputs.

**Prerequisites:**
* You must have **Miniforge** (or Mambaforge) installed.
* You need to have a few directories set up:
    - `./results`, with subdirs with '_results' appended that contain gff3/gff files with annotations from each software. For example, for Helixer, you will have `./results/helixer_results`. If you would like to indicate a particular difference between two software, you can add an extra underscores (e.g. `results/tiberius_results` vs. `results/tiberius_softmasked_results`).
    - `./genomes` and `./genomes/<particular difference>`. Store fasta files for genome here. For `./genomes/<particular difference>`, a good example of use is `./genomes/softmasked` to indicate the genomes used for the `./results/tiberius_softmasked_results`. The text for the "particular difference" must match.
    - `./reference_annotations`, where you store the gffs with reference annotations for each and every genome. These must be indicated in the scripts below too, before running.

## **Git Repository Note**
Due to GitHub file size limits, the large genome files (`*.fa`) in `genomes/` and `genomes/softmasked/` are **not tracked** in this repository.
You must download these files separately and place them in the `genomes/` directory before running the pipeline.
See `genomes/Download_1043824_Publications.csv` for download information.

---

## **1. PASA Environment (`pasa_env`)**
**Purpose:** Pre-processing, repairing GTF/GFF formats, and extracting protein sequences.
**NOTE**: You must already have PASA installed. For these results, pasa-2.5.3 was used.

### **Installation**
Copy and paste this entire block into your terminal:

```bash
# 1. Create environment
mamba create -n pasa_env python=3.10 -y
# Run this for MacOS instead:
# CONDA_SUBDIR=osx-64 mamba create -n pasa_env python=3.10 -y

# 2. Activate the environment
mamba activate pasa_env

# For MacOS only:
# 3. PERMANENTLY lock this environment to Intel architecture
# This ensures any future 'mamba install' uses the correct x86 packages
# conda config --env --set subdir osx-64

# 4. Install Dependencies
# - gffread: Standard tool for converting GTF -> GFF3
# - perl: Required for PASA
# - perl-bioperl: Required for gff3_file_to_proteins.pl
# - perl-db-file: Common database requirement for PASA scripts
mamba install -c bioconda -c conda-forge gffread perl perl-bioperl perl-db-file -y

# For MacOS:
# 5. Verify Architecture (Should say 'x86_64', NOT 'arm64')
# python -c "import platform; print(platform.machine())"

```

### **Scripts for this Env**

* `standardize_gff.py` (runs on all GFFs under `./results`, organizing them in a common format under `./standardized_results` to be used in the rest of pipeline)
* `fix_tiberius.py` (convert GTFs to GFFs, particularly for Tiberius)
* `match_header.py` (run to match GFF output and genome FASTA input headers, to not confuse PASA)
* `run_protein_conversion.sh` (Protein extraction with PASA) to `./protein_results`)

---

## **2. BUSCO Environment (`busco_env`)**

**Purpose:** Quality Control (BUSCO) and Accuracy Benchmarking (gffcompare).

### **Installation**

Copy and paste this entire block into your terminal:

```bash
# 1. Creates environment
mamba create -n busco_env -c conda-forge -c bioconda busco -y
# For MacOS:
# CONDA_SUBDIR=osx-64 mamba create -n busco_env -c conda-forge -c bioconda busco -y

# 2. Activate
mamba activate busco_env

# For MacOS:
# 3. Lock architecture
# conda config --env --set subdir osx-64

# 4. Install gffcompare
mamba install -c bioconda gffcompare -y

# 5. Verify Installation
echo "--- BUSCO Version ---"
busco --version
echo "--- gffcompare Version ---"
gffcompare --version
# For MacOS:
echo "--- Architecture (Must be x86_64) ---"
python -c "import platform; print(platform.machine())"

```

### **Scripts for this Env**

* `setup_for_busco.sh` (Download busco lineage model, removes bad lines from GFFs, like comments)
* `run_busco_annotations.sh` (Runs the assessment)
* `summarize_busco_v3.py` (Creates the stats table in a csv)
* `run_gffcompare.sh` (Runs the accuracy comparison against reference. Edit reference paths in script)
* `summarize_gffcompare.py` (Creates the sensitivity/precision table)

---
Note: for some versions of Busco, you might need to run `sed -i '' 's/\*//g' *.fa` to clean all the stop codons ('*') from the output prot.fastas; otherwise Busco will give an error.