import os
import json
import csv
import glob

# --- Configuration ---
ROOT_DIR = "./"
BUSCO_RESULTS_DIR = os.path.join(ROOT_DIR, "busco_results")
OUTPUT_CSV = os.path.join(ROOT_DIR, "busco_summary_final.csv")

def get_stats_from_json(json_path):
    """Safely extracts counts from BUSCO v6 JSON."""
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            
        res = data.get('results', {})
        
        # --- BUSCO v6 SPECIFIC KEYS ---
        # We use .get(key, 0) to prevent crashes if a key is missing
        s = res.get('Single copy BUSCOs', 0)
        d = res.get('Multi copy BUSCOs', 0)
        f = res.get('Fragmented BUSCOs', 0)
        m = res.get('Missing BUSCOs', 0)
        total = res.get('n_markers', 0)
        c_pct = res.get('Complete percentage', 0.0)

        # Fallback: If 'Complete percentage' is missing, calculate it
        if c_pct == 0.0 and total > 0:
            c_pct = round(((s + d) / total) * 100, 1)

        return {
            "C_pct": c_pct,
            "S": s,
            "D": d,
            "F": f,
            "M": m,
            "Total": total
        }
    except Exception as e:
        print(f"[WARN] Failed to parse {os.path.basename(json_path)}: {e}")
        return None

def main():
    print(f"Scanning {BUSCO_RESULTS_DIR} for BUSCO v6 results...")

    # Look for 'short_summary.specific*.json' files recursively
    search_pattern = os.path.join(BUSCO_RESULTS_DIR, "**", "short_summary.specific*.json")
    json_files = glob.glob(search_pattern, recursive=True)

    if not json_files:
        print("[ERROR] No JSON files found.")
        return

    # CSV Headers
    headers = ['Species', 'Tool', 'Complete (%)', 'Single (S)', 'Duplicated (D)', 'Fragmented (F)', 'Missing (M)', 'Total BUSCOs']
    
    unique_rows = {}

    for jf in json_files:
        # --- 1. Infer Name from Folder Structure ---
        # Path: busco_results/B_distachyon_tiberius/run_viridiplantae.../short_summary.json
        
        # Get the parent folder name (e.g., "run_viridiplantae_odb10")
        parent_dir = os.path.basename(os.path.dirname(jf))
        
        # If we are in the 'run_' folder, the sample name is one level up
        if parent_dir.startswith("run_"):
            sample_folder = os.path.basename(os.path.dirname(os.path.dirname(jf)))
        else:
            sample_folder = os.path.basename(os.path.dirname(jf))

        # Split "Species_Tool"
        if "_" in sample_folder:
            species, tool = sample_folder.rsplit("_", 1)
        else:
            species, tool = sample_folder, "Unknown"

        # --- 2. Extract Stats ---
        stats = get_stats_from_json(jf)
        
        if stats and stats['Total'] > 0:
            # Save to dictionary to deduplicate (in case of multiple summary files per run)
            key = f"{species}_{tool}"
            unique_rows[key] = {
                'Species': species,
                'Tool': tool,
                'Complete (%)': stats['C_pct'],
                'Single (S)': stats['S'],
                'Duplicated (D)': stats['D'],
                'Fragmented (F)': stats['F'],
                'Missing (M)': stats['M'],
                'Total BUSCOs': stats['Total']
            }
            print(f"  -> Found: {species} ({tool})")

    # --- 3. Write CSV ---
    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        
        # Sort alphabetically by species
        for key in sorted(unique_rows.keys()):
            writer.writerow(unique_rows[key])

    print("------------------------------------------------")
    print(f"[SUCCESS] Table saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()