import os
import sys
import re


# Canonical species prefixes so annevo/helixer (and others) match
SPECIES_CANONICAL = {
    "Athaliana": "A_thaliana",
    "A_thaliana": "A_thaliana",
    "Bdistachyon": "B_distachyon",
    "B_distachyon": "B_distachyon",
    "BrapaO": "B_rapa",
    "B_rapa": "B_rapa",
    "Brapa": "B_rapa",
    "Gmax": "G_max",  # keep common shorthand without underscore
    "Lsativa": "L_sativa",
    "Osativa": "O_sativa",
    "O_sativa": "O_sativa",
    "Ppatens": "P_patens",
    "P_patens": "P_patens",
    "Ptrichocarpa": "P_trichocarpa",
    "P_trichocarpa": "P_trichocarpa",
    "Sbicolor": "S_bicolor",
    "S_bicolor": "S_bicolor",
    "Sitalica": "S_italica",
    "S_italica": "S_italica",
    "Slycopersicum": "S_lycopersicum",
    "Vvinifera": "V_vinifera",
    "V_vinifera": "V_vinifera",
    "Zmays": "Z_mays",
    "Z_mays": "Z_mays",
}


def parse_attributes(attr_string):
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes
    
    # Check if it's GTF style (key "value";)
    if ' "' in attr_string or ' ";' in attr_string:
        parts = attr_string.strip().split(';')
        for part in parts:
            part = part.strip()
            if not part: continue
            if ' ' in part:
                try:
                    key, value = part.split(' ', 1)
                    value = value.strip('"')
                    attributes[key] = value
                except ValueError:
                    pass
    else:
        # GFF3 style (key=value;)
        parts = attr_string.strip().split(';')
        for part in parts:
            part = part.strip()
            if not part: continue
            if '=' in part:
                try:
                    key, value = part.split('=', 1)
                    attributes[key] = value
                except ValueError:
                    pass
    return attributes

def format_attributes(attributes):
    if not attributes:
        return "."
    # Ensure ID and Parent are first
    ordered = []
    if 'ID' in attributes:
        ordered.append(f"ID={attributes['ID']}")
    if 'Parent' in attributes:
        ordered.append(f"Parent={attributes['Parent']}")
    
    for k, v in attributes.items():
        if k not in ['ID', 'Parent']:
            ordered.append(f"{k}={v}")
            
    return ";".join(ordered)


def infer_species(tool_name, filename):
    """Infer a clean species name from tool + filename.

    The returned string is used as the sole prefix for all feature IDs
    and as the species part of the standardized output filename.
    """

    base = os.path.basename(filename)

    # Helixer patterns, e.g. A_thaliana_annotated_helixer.gff3, Z_mays_helixer.gff3
    if tool_name == "helixer":
        if base.endswith("_annotated_helixer.gff3"):
            species = base[:-len("_annotated_helixer.gff3")]
        elif base.endswith("_annotated_helixer.gff"):
            species = base[:-len("_annotated_helixer.gff")]
        elif base.endswith("_helixer.gff3"):
            species = base[:-len("_helixer.gff3")]
        elif base.endswith("_helixer.gff"):
            species = base[:-len("_helixer.gff")]
        else:
            species = os.path.splitext(base)[0]

    # ANNEVO patterns, e.g. Athaliana_167_TAIR10.fa_annotated.gff
    elif tool_name == "annevo":
        tmp = re.sub(r"\.gff3?$", "", base)
        # Handle .fa_annotated, .fa_annnotated (typo) and similar
        tmp = re.sub(r"\.fa_a+n+otated$", "", tmp)
        # Now keep only the first token before an underscore as species
        species = tmp.split("_")[0]

    # Tiberius or other tools: usually Species.ext
    else:
        # Strip extension and take first token before underscore or dot
        tmp = os.path.splitext(base)[0]
        species = re.split(r"[_.]", tmp)[0]

    # Clean up characters: species name only (letters, digits, underscore)
    species = re.sub(r"[^A-Za-z0-9_]", "_", species)
    species = re.sub(r"_+", "_", species).strip("_")

    # Normalize to a canonical prefix so tools agree
    if species in SPECIES_CANONICAL:
        species = SPECIES_CANONICAL[species]

    return species or "species"


def standardize_gff(input_file, output_file, species_prefix, tool_name):
    print(f"Processing {input_file}...")
    features = []
    fixes = []
    
    # 1. Read all features
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            
            seqid, source, type_, start, end, score, strand, phase, attributes_str = parts
            attributes = parse_attributes(attributes_str)
            
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            if start > end:
                start, end = end, start
                fixes.append(f"Line {line_num}: Swapped coordinates {start}-{end}")
            
            # Normalize type
            original_type = type_
            if type_ == 'transcript':
                type_ = 'mRNA'
                fixes.append(f"Line {line_num}: Normalized type 'transcript' to 'mRNA'")
            
            # GTF to GFF3 ID/Parent logic for initial linking
            if 'ID' not in attributes:
                if 'gene_id' in attributes:
                    if type_ == 'gene':
                        attributes['ID'] = attributes['gene_id']
                    elif type_ == 'mRNA':
                        attributes['ID'] = attributes.get('transcript_id', attributes['gene_id'] + ".t1")
                        attributes['Parent'] = attributes['gene_id']
                    else:
                        # exon, CDS, UTR
                        attributes['Parent'] = attributes.get('transcript_id', attributes['gene_id'])
                        # Temp ID
                        attributes['ID'] = f"temp_{line_num}"
                else:
                    # No gene_id, generate temp ID
                    attributes['ID'] = f"temp_{line_num}"
            
            feature = {
                'seqid': seqid,
                'source': source,
                'type': type_,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attributes,
                'children': [],
                'original_line': line_num
            }
            features.append(feature)

    # 2. Build Hierarchy
    genes = {}
    mrnas = {}
    orphans = []
    
    # Index by ID
    feature_map = {f['attributes']['ID']: f for f in features}
    
    # Separate by type
    gene_features = [f for f in features if f['type'] == 'gene']
    mrna_features = [f for f in features if f['type'] == 'mRNA']
    other_features = [f for f in features if f['type'] not in ['gene', 'mRNA']]
    
    # Populate genes
    for g in gene_features:
        genes[g['attributes']['ID']] = g

    # Link mRNAs to genes
    for m in mrna_features:
        parent_id = m['attributes'].get('Parent')
        if parent_id and parent_id in genes:
            genes[parent_id]['children'].append(m)
        else:
            # Orphan mRNA - create a gene for it
            new_gene_id = f"gene_for_{m['attributes']['ID']}"
            new_gene = {
                'seqid': m['seqid'],
                'source': m['source'],
                'type': 'gene',
                'start': m['start'],
                'end': m['end'],
                'score': '.',
                'strand': m['strand'],
                'phase': '.',
                'attributes': {'ID': new_gene_id},
                'children': [m],
                'original_line': 'generated'
            }
            genes[new_gene_id] = new_gene
            m['attributes']['Parent'] = new_gene_id
            fixes.append(f"Created parent gene for orphan mRNA {m['attributes']['ID']}")

    # Link exons/CDS to mRNAs
    for o in other_features:
        parent_id = o['attributes'].get('Parent')
        if parent_id and parent_id in feature_map:
            parent = feature_map[parent_id]
            if parent['type'] == 'mRNA':
                parent['children'].append(o)
            elif parent['type'] == 'gene':
                # Direct child of gene? Move to a new mRNA
                # Check if gene already has an mRNA that covers this feature?
                # For simplicity, create a new mRNA for this gene
                # But wait, if a gene has exons directly, it's weird.
                # Let's look for an existing mRNA in this gene
                found_mrna = False
                for child in parent['children']:
                    if child['type'] == 'mRNA':
                        # Assign to first mRNA? Risky.
                        # Let's create a new mRNA
                        pass
                
                # Create synthetic mRNA
                new_mrna_id = f"mRNA_for_{o['attributes']['ID']}"
                new_mrna = {
                    'seqid': o['seqid'],
                    'source': o['source'],
                    'type': 'mRNA',
                    'start': o['start'],
                    'end': o['end'],
                    'score': '.',
                    'strand': o['strand'],
                    'phase': '.',
                    'attributes': {'ID': new_mrna_id, 'Parent': parent['attributes']['ID']},
                    'children': [o],
                    'original_line': 'generated'
                }
                parent['children'].append(new_mrna)
                o['attributes']['Parent'] = new_mrna_id
                fixes.append(f"Created intermediate mRNA for feature {o['type']} attached to gene {parent['attributes']['ID']}")
            else:
                # Parent is something else?
                orphans.append(o)
        else:
            # No parent or parent not found
            # Try to find an mRNA that overlaps?
            # For now, treat as orphan
            orphans.append(o)
            fixes.append(f"Orphan feature {o['type']} at {o['seqid']}:{o['start']} (Parent={parent_id}) dropped")

    # 3. Verify and Fix Coordinates
    sorted_genes = sorted(genes.values(), key=lambda x: (x['seqid'], x['start']))
    
    for gene in sorted_genes:
        # Update gene coordinates to cover all children
        min_start = gene['start']
        max_end = gene['end']
        
        for mrna in gene['children']:
            # Update mRNA coordinates to cover all exons
            m_min = mrna['start']
            m_max = mrna['end']
            
            if mrna['children']:
                child_starts = [c['start'] for c in mrna['children']]
                child_ends = [c['end'] for c in mrna['children']]
                m_min = min(child_starts)
                m_max = max(child_ends)
            
            if m_min < mrna['start'] or m_max > mrna['end']:
                fixes.append(f"Expanded mRNA {mrna['attributes']['ID']} coordinates from {mrna['start']}-{mrna['end']} to {m_min}-{m_max}")
                mrna['start'] = m_min
                mrna['end'] = m_max
            
            min_start = min(min_start, mrna['start'])
            max_end = max(max_end, mrna['end'])
            
        if min_start < gene['start'] or max_end > gene['end']:
            fixes.append(f"Expanded gene {gene['attributes']['ID']} coordinates from {gene['start']}-{gene['end']} to {min_start}-{max_end}")
            gene['start'] = min_start
            gene['end'] = max_end

    # 4. Uniform Renaming: IDs use species name ONLY as prefix
    #    (no chromosome, model, or file-specific suffixes).
    prefix = species_prefix
    
    gene_counter = 0
    
    for gene in sorted_genes:
        gene_counter += 1
        new_gene_id = f"{prefix}_g{gene_counter:06d}"
        gene['attributes']['ID'] = new_gene_id
        # Remove old ID/Parent/Name to avoid confusion, or keep as Alias?
        # User wants uniform IDs.
        
        mrna_counter = 0
        for mrna in gene['children']:
            mrna_counter += 1
            new_mrna_id = f"{new_gene_id}.t{mrna_counter}"
            mrna['attributes']['ID'] = new_mrna_id
            mrna['attributes']['Parent'] = new_gene_id
            
            # Sort children by start
            mrna['children'].sort(key=lambda x: x['start'])
            
            # Counters for different types
            type_counters = {}
            
            for child in mrna['children']:
                ctype = child['type']
                if ctype not in type_counters:
                    type_counters[ctype] = 0
                type_counters[ctype] += 1
                
                new_child_id = f"{new_mrna_id}.{ctype}.{type_counters[ctype]}"
                child['attributes']['ID'] = new_child_id
                child['attributes']['Parent'] = new_mrna_id

    # Record global normalization fixes (once per file)
    fixes.append(
        f"Standardized all gene/mRNA/exon/CDS IDs to use species prefix '{prefix}'"
    )
    fixes.append(
        "Rebuilt and verified Parent relationships to match gene→mRNA→subfeature hierarchy"
    )
    fixes.append("Sorted genes and all child features by genomic coordinates")

    # 5. Write Output
    with open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        if fixes:
            out.write("# Fixes applied:\n")
            # Deduplicate fixes
            unique_fixes = sorted(list(set(fixes)))
            for fix in unique_fixes[:20]:
                out.write(f"# - {fix}\n")
            if len(unique_fixes) > 20:
                out.write(f"# - ... and {len(unique_fixes) - 20} more\n")
        else:
            out.write("# No fixes applied.\n")
            
        for gene in sorted_genes:
            write_feature(out, gene)


def write_feature(out, feature):
    attr_str = format_attributes(feature['attributes'])
    line = f"{feature['seqid']}\t{feature['source']}\t{feature['type']}\t{feature['start']}\t{feature['end']}\t{feature['score']}\t{feature['strand']}\t{feature['phase']}\t{attr_str}\n"
    out.write(line)
    
    if 'children' in feature:
        # Sort children: mRNA first (though usually inside gene), but here children are mRNAs
        # Inside mRNA, children are exons/CDS.
        # Standard GFF3: Gene -> mRNA -> Exon/CDS
        for child in feature['children']:
            write_feature(out, child)


def main():
    root_dir = "/Users/filipr/Desktop/transgenic"
    results_root = os.path.join(root_dir, "results")
    output_dir = os.path.join(root_dir, "standardized_results")
    os.makedirs(output_dir, exist_ok=True)
    
    result_folders = [
        d
        for d in os.listdir(results_root)
        if os.path.isdir(os.path.join(results_root, d))
    ]

    for folder in result_folders:
        folder_path = os.path.join(results_root, folder)
        tool_name = folder.replace("_results", "")

        files = [
            f
            for f in os.listdir(folder_path)
            if f.endswith(".gff") or f.endswith(".gff3") or f.endswith(".gtf")
        ]

        for file_name in files:
            input_path = os.path.join(folder_path, file_name)

            # Infer species from tool + filename
            species = infer_species(tool_name, file_name)

            # Standardize output filename: <species>_<tool>.gff3 (no extra characters)
            clean_species = species
            clean_tool = re.sub(r"[^A-Za-z0-9_]", "_", tool_name)
            clean_tool = re.sub(r"_+", "_", clean_tool).strip("_")

            output_basename = f"{clean_species}_{clean_tool}.gff3"
            output_path = os.path.join(output_dir, output_basename)

            standardize_gff(input_path, output_path, species, tool_name)

if __name__ == "__main__":
    main()
