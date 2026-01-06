busco --download viridiplantae_odb10 \
  --download_path "/Users/filipr/Desktop/transgenic_comparison/busco_downloads"
cd /Users/filipr/Desktop/transgenic_comparison/protein_results

# Loop through files and remove lines starting with '#'
# sed -i '' '/^#/d' deletes any line that starts with #
for f in *.prot.fasta; do
    sed -i '' '/^#/d' "$f"
    echo "Fixed header in $f"
done