import sys
import os

def get_gff_ids(gff_path):
    """Extracts unique sequence IDs (column 1) from a GFF3 file."""
    ids = set()
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) > 0:
                ids.add(parts[0]) # Column 1 is the SeqID
    return ids

def clean_fasta(gff_ids, fasta_input, fasta_output):
    """
    Reads FASTA. If the first word of the header matches a GFF ID, 
    it truncates the rest of the header.
    """
    found_count = 0
    missing_count = 0
    
    with open(fasta_input, 'r') as fin, open(fasta_output, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                # remove '>' and newlines
                full_header = line.strip()[1:] 
                
                # Logic: Take the first word (everything before the first space)
                # This fixes ">Chr1 CHROMOSOME..." -> "Chr1"
                potential_id = full_header.split()[0]
                
                if potential_id in gff_ids:
                    fout.write(f">{potential_id}\n")
                    found_count += 1
                else:
                    # If the GFF doesn't have this ID, we keep the original 
                    # (or you could choose to skip it)
                    fout.write(line)
                    missing_count += 1
            else:
                fout.write(line)
    
    return found_count, missing_count

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python match_headers.py <gff_file> <input_fasta> <output_fasta>")
        sys.exit(1)
    
    gff_file = sys.argv[1]
    in_fasta = sys.argv[2]
    out_fasta = sys.argv[3]
    
    print(f"[Py] Reading GFF IDs from {os.path.basename(gff_file)}...")
    ids = get_gff_ids(gff_file)
    print(f"[Py] Found {len(ids)} unique sequence IDs in GFF.")
    
    print(f"[Py] Cleaning FASTA headers in {os.path.basename(in_fasta)}...")
    matched, missed = clean_fasta(ids, in_fasta, out_fasta)
    print(f"[Py] Done. {matched} headers matched and cleaned. {missed} unmatched.")