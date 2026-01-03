import sys
import os
import glob
import re

def parse_tiberius_line(line):
    parts = line.strip().split('\t')
    if len(parts) < 9:
        return None

    seqid, source, feature, start, end, score, strand, frame = parts[:8]
    attributes_raw = parts[8].strip()

    # --- LOGIC 1: Handle "Naked ID" lines (Gene / Transcript) ---
    # Example: "g1" or "g1.t1"
    if feature in ['gene', 'transcript']:
        
        # 1. Rename 'transcript' to 'mRNA' for PASA
        if feature == 'transcript':
            feature = 'mRNA'
            
        # 2. Construct ID and Parent
        # We assume the ID is just the raw text in column 9
        curr_id = attributes_raw
        
        if feature == 'gene':
            # Genes have no parent
            new_attrs = f"ID={curr_id}"
        else:
            # mRNA needs a Parent. 
            # Logic: "g1.t1" -> Parent is "g1" (strip after last period)
            if '.' in curr_id:
                parent_id = curr_id.rsplit('.', 1)[0]
            else:
                # Fallback if no dot exists
                parent_id = curr_id
            new_attrs = f"ID={curr_id};Parent={parent_id}"

    # --- LOGIC 2: Handle "GTF Style" lines (CDS / Exon / Intron) ---
    # Example: transcript_id "g1.t1"; gene_id "g1"; ...
    else:
        # We only strictly need CDS and exon for PASA.
        # Introns are usually inferred, but we can keep them if valid.
        
        # Extract transcript_id using Regex to be safe
        # Looks for: transcript_id "MATCH";
        match = re.search(r'transcript_id "([^"]+)"', attributes_raw)
        
        if match:
            parent_id = match.group(1)
            # CDS/Exons link to the mRNA (the transcript_id)
            new_attrs = f"Parent={parent_id}"
        else:
            # If we can't find a transcript_id, skip this line (it's orphaned)
            return None

    # Return the reconstructed GFF3 line
    return "\t".join([seqid, source, feature, start, end, score, strand, frame, new_attrs])

def convert_file(input_path, output_path):
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            # Skip comments or empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            try:
                new_line = parse_tiberius_line(line)
                if new_line:
                    fout.write(new_line + "\n")
            except Exception as e:
                pass # Silently skip malformed lines to prevent crashing

if __name__ == "__main__":
    # --- CONFIGURATION ---
    # Adjust these paths to match your folders
    OUTPUT_DIR = "/Users/filipr/Desktop/transgenic_comparison/tiberius_softmasked_gff3s"
    INPUT_DIR = "/Users/filipr/Desktop/transgenic_comparison/gtfs/tiberius_softmasked_results"

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Find all files (assuming they end in .gtf or .gff)
    files = glob.glob(os.path.join(INPUT_DIR, "*"))
    
    print(f"Found {len(files)} files. Starting conversion...")

    count = 0
    for f in files:
        if f.endswith(".DS_Store"): continue # Skip Mac system files
        
        # Output filename: same name but forced .gff3 extension
        base_name = os.path.splitext(os.path.basename(f))[0]
        out_file = os.path.join(OUTPUT_DIR, base_name + ".gff3")
        
        convert_file(f, out_file)
        count += 1
        print(f"  [{count}] Converted: {base_name} -> .gff3")

    print("Done.")