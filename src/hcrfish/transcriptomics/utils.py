import pandas as pd 

def get_sequence(seq_dict, chrom, start, end, strand):
    # Check that chromosome is in the sequence dictionary
    if chrom not in seq_dict:
        raise ValueError(f"Chromosome {chrom} not found in seq_dict.")
    
    # Check that start < end 
    if start >= end:
        raise ValueError(f"Start position {start} is greater than or equal to end position {end}.")
    
    # Get the sequence from the dictionary
    sequence = seq_dict[chrom].seq[start:end]
    
    # Reverse complement if necessary
    if strand == '-':
        sequence = sequence.reverse_complement()
    
    # Return the sequence as a string in uppercase
    return str(sequence).upper()


def gtf_to_dataframe(gtf_file):
    # Define column names as per GTF format plus additional attribute columns
    columns = [
        "seqname",   # Chromosome or scaffold
        "source",    # Source of the feature (e.g., Ensembl, HAVANA)
        "feature",   # Feature type (e.g., gene, transcript, exon)
        "start",     # Start position of the feature
        "end",       # End position of the feature
        "score",     # Score (can be '.')
        "strand",    # Strand (+ or -)
        "frame",     # Frame (0, 1, or 2, or '.')
        "gene_id",   # Extracted gene_id
        "gene_name", # Extracted gene_name
        "transcript_id", # Extracted transcript_id
        "attributes" # Full attributes dictionary
    ]

    # Read the GTF file line by line and parse into DataFrame rows
    data = []
    with open(gtf_file, 'r') as file:
        for line in file:
            # Skip comments
            if line.startswith("#"):
                continue
            
            # Split the line into fields and extract the attributes column
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip if there are not enough fields
            
            # Parse the attributes column into a dictionary
            attributes = {}
            for attribute in fields[8].strip().split(';'):
                if attribute.strip():
                    key_value = attribute.strip().split(' ', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        attributes[key] = value.replace('"', '')
            
            # Extract specific attributes if they exist, otherwise use None
            gene_id = attributes.get("gene_id")
            gene_name = attributes.get("gene_name")
            transcript_id = attributes.get("transcript_id")
            
            # Add each field, specific attributes, and the full attributes dictionary as a row
            row = fields[:8] + [gene_id, gene_name, transcript_id, attributes]
            data.append(row)

    # Create a DataFrame with the parsed data
    df = pd.DataFrame(data, columns=columns)

    # Convert start and end to integers
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    return df



    