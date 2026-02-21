busco --download viridiplantae_odb10 \
  --download_path "/Users/filipr/Desktop/transgenic_comparison/busco_downloads"
cd /Users/filipr/Desktop/transgenic_comparison/protein_results

# Loop through files and clean for BUSCO compatibility
for f in *.prot.fasta; do
    # Remove comment lines starting with '#'
    sed -i '' '/^#/d' "$f"
    # Remove '/' characters from sequences (illegal in hmmsearch)
    # Only modifies non-header lines
    sed -i '' '/^[^>]/s/\///g' "$f"
    echo "Cleaned $f for BUSCO"
done