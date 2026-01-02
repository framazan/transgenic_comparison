import os
import csv
import glob

# --- Configuration ---
ROOT_DIR = "/Users/filipr/Desktop/transgenic"
STATS_DIR = os.path.join(ROOT_DIR, "gffcompare_results")
OUTPUT_CSV = os.path.join(ROOT_DIR, "gffcompare_summary.csv")

# List of known species to correctly split filenames
# (Longer names first to avoid partial matches, e.g., finding 'P_patens' inside 'P_patens_complex')
KNOWN_SPECIES = [
    "A_thaliana", "B_distachyon", "B_rapa", "G_max", "L_sativa", 
    "O_sativa", "P_patens", "P_trichocarpa", "S_bicolor", "S_italica", 
    "S_lycopersicum", "V_vinifera", "Z_mays"
]

def parse_stats_file(filepath):
    """Parses .stats file to extract Accuracy metrics, ignoring pipes."""
    stats = {
        "Base_Sens": 0.0, "Base_Prec": 0.0,
        "Exon_Sens": 0.0, "Exon_Prec": 0.0,
        "Trans_Sens": 0.0, "Trans_Prec": 0.0
    }
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Identify the row
            if line.startswith("Base level:"):
                key_sens, key_prec = "Base_Sens", "Base_Prec"
            elif line.startswith("Exon level:"):
                key_sens, key_prec = "Exon_Sens", "Exon_Prec"
            elif line.startswith("Transcript level:"):
                key_sens, key_prec = "Trans_Sens", "Trans_Prec"
            else:
                continue

            # Parsing Logic:
            # Original Line: "Base level:    76.6     |    97.7    |"
            # 1. Remove "Base level:" label -> "76.6     |    97.7    |"
            content = line.split(":", 1)[1]
            
            # 2. Replace pipes with spaces -> "76.6          97.7     "
            content_clean = content.replace("|", " ")
            
            # 3. Split by whitespace -> ['76.6', '97.7']
            parts = content_clean.split()
            
            if len(parts) >= 2:
                stats[key_sens] = parts[0]
                stats[key_prec] = parts[1]

    return stats

def get_species_and_tool(filename_base):
    """Intelligently splits filename into Species and Tool."""
    for species in KNOWN_SPECIES:
        if filename_base.startswith(species):
            # The tool is everything AFTER the species name (and the underscore)
            # Example: "A_thaliana_tiberius_softmasked"
            # Species: "A_thaliana"
            # Tool: "tiberius_softmasked"
            
            # Check if there's actually a tool suffix
            if len(filename_base) > len(species):
                # Remove species + 1 char (the underscore)
                tool = filename_base[len(species)+1:]
                return species, tool
            else:
                return species, "Unknown"
    
    # Fallback if species isn't in list
    if "_" in filename_base:
        return filename_base.split("_", 1)
    return filename_base, "Unknown"

def main():
    print(f"Scanning {STATS_DIR}...")
    stat_files = glob.glob(os.path.join(STATS_DIR, "*.stats"))

    if not stat_files:
        print("[ERROR] No .stats files found.")
        return

    headers = [
        "Species", "Tool", 
        "Base_Sensitivity", "Base_Precision",
        "Exon_Sensitivity", "Exon_Precision",
        "Transcript_Sensitivity", "Transcript_Precision"
    ]

    rows = []
    for sf in stat_files:
        filename = os.path.basename(sf) # e.g. A_thaliana_tiberius_softmasked.stats
        base = filename.replace(".stats", "")
        
        # 1. Correct Species/Tool splitting
        species, tool = get_species_and_tool(base)

        # 2. Correct Stats Parsing
        metrics = parse_stats_file(sf)
        
        rows.append({
            "Species": species,
            "Tool": tool,
            "Base_Sensitivity": metrics["Base_Sens"],
            "Base_Precision": metrics["Base_Prec"],
            "Exon_Sensitivity": metrics["Exon_Sens"],
            "Exon_Precision": metrics["Exon_Prec"],
            "Transcript_Sensitivity": metrics["Trans_Sens"],
            "Transcript_Precision": metrics["Trans_Prec"]
        })
        print(f"  -> Parsed: {species} ({tool})")

    # Write CSV
    with open(OUTPUT_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        
        # Sort by Species then Tool
        rows.sort(key=lambda x: (x['Species'], x['Tool']))
        
        for r in rows:
            writer.writerow(r)

    print("------------------------------------------------")
    print(f"[SUCCESS] Table saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()