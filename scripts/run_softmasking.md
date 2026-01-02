This softmasking pipeline was used to annotate **Arabidopsis thaliana** and **Vitis vinifera**, which Phytozome did not have. This pipeline was taken from what Tiberius recommends here.

### **1. Installation & Setup**

We use **Singularity (Apptainer)** to encapsulate the complex dependencies (RepeatModeler, RepeatMasker, TRF, RECON, RMBlast) into a single, portable file.

#### **A. Pull the Container**

Run these commands on the login node to download the Dfam TE Tools container.

```bash
# 1. Create a clean workspace
mkdir -p $HOME/softmasked/bin
cd $HOME/softmasked

# 2. Pull the official Dfam TE Tools container
singularity pull tetools.sif docker://dfam/tetools:1.89

```

#### **B. Install External Dependencies**

The pipeline requires a few standalone tools for the final steps (TRF refinement). We install them locally since the HPC lacks them or uses outdated versions.

**1. GNU Parallel (For multicore speed)**
Make sure that HPC has GNU Parallel. Script searches and uses the installation, no matter where it may be. Run this to install it if needed:
```bash
cd bin
wget https://ftp.gnu.org/gnu/parallel/parallel-20230722.tar.bz2
tar -xjf parallel-20230722.tar.bz2
cd parallel-20230722
./configure --prefix=$PWD/..
make && make install
cd ../..

```

**2. Helper Scripts**
We use two helper scripts to handle the splitting and parsing of TRF data. These are taken from Tiberius' documentation.

* **Create `bin/splitMfasta.pl**`
* **Create `bin/parseTrfOutput.py**`

* **Permissions:**
```bash
chmod +x bin/splitMfasta.pl bin/parseTrfOutput.py

```

---

### **2. Execution**

We run the pipeline using a **SLURM Job Array**. This allows us to process multiple genomes simultaneously using a single script.

#### **A. The Directory Structure**

Ensure your working directory looks like this:

```text
softmasked/
├── tetools.sif                # The container
├── bin/                       # Local tools (parallel, scripts)
├── genomes/                   # Input FASTA files
│   ├── Athaliana_167_TAIR10.fa
│   └── Vvinifera_T2T_ref.fa
└── submit_masking.sh          # The script below

```

#### **B. The Submission Script**

Save as `repeat_masking.sh`.  It dynamically adjusts resources and handles the entire workflow (Modeling  Masking  TRF). You can view it in this file's working directory.

### **3. Final Outputs**

Your results will appear in `./softmasked`, as specified in `repeat_masking.slurm`:

* `*.final_softmasked.fasta`: The genome with repeats lower-cased (e.g., `atgc`).
* `*.trf.gff3`: The location of simple tandem repeats.
* `*-families.fa`: The library of TE families discovered in your specific genome.

Use this with Tiberius, Augustus, or whatever other annotation software that requires it.