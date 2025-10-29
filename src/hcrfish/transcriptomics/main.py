from tqdm import tqdm 
from Bio import SeqIO
from hcrfish.transcriptomics.classes import Transcriptome, Gene, Transcript, Exon, Intron, UTR, CDS 
from hcrfish.transcriptomics.utils import get_sequence, gtf_to_dataframe

def generate_transcriptome_object(transcriptome_gtf_path, genome_fasta_path): 
    
    # Initialize the genome sequence
    genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fasta_path, "fasta"))

    # Initialize the Transcriptome object
    transcriptome = Transcriptome()

    # Create a database from the GTF file
    db = gtf_to_dataframe(transcriptome_gtf_path)

    # Remove columns source, score, frame and attributes 
    db = db.drop(columns=['source', 'score', 'frame', 'attributes'])

    # Sort the DataFrame by `gene_id` for efficient processing
    db = db.sort_values(by='gene_id').reset_index(drop=True)

    # Group by `gene_id`, this avoids the need for repeated filtering
    grouped_genes = db.groupby('gene_id')

    # Collect all unique genes based on gene_id column 
    unique_gene_ids = db['gene_id'].unique()
    print(f"Found {len(unique_gene_ids)} unique genes.")

    # Convert all "features" to lower case 
    db['feature'] = db['feature'].str.lower()

    # Iterate through each group of gene_id
    for gene_id, gene_rows in tqdm(grouped_genes, total=len(unique_gene_ids)):

        # Get the gene name (if available)
        gene_name = gene_rows['gene_name'].dropna().iloc[0] if gene_rows['gene_name'].notna().any() else gene_id

        # Get chromosome and skip if not in genome 
        chromosome = gene_rows['seqname'].dropna().iloc[0]
        if chromosome not in genome_seq.keys():
            continue 

        # Get strand and skip if not valid
        strand = gene_rows['strand'].dropna().iloc[0]
        if strand not in ['+', '-']:   
                    continue

        # Create a Gene object
        gene = Gene(gene_id, gene_name)
        gene.chromosome = chromosome
        gene.strand = strand

        # Group rows by gene_id and transcript_id 
        gene_groups = gene_rows.groupby(['gene_id', 'transcript_id'])

        for (gene_id, transcript_id), transcript_rows in gene_groups:
            # If transcript_id is not present, skip
            if transcript_id == "":
                continue

            # Create a Transcript object
            transcript = Transcript(transcript_id)
            transcript.chromosome = chromosome
            transcript.strand = strand

            # Collect all exons, CDS, 5' UTR, and 3' UTR features for this transcript
            exons = transcript_rows[transcript_rows['feature'] == 'exon']
            cds = transcript_rows[transcript_rows['feature'] == 'cds']
            utr5 = transcript_rows[(transcript_rows['feature'] == '5utr') | (transcript_rows['feature'] == 'five_prime_utr')]
            utr3 = transcript_rows[(transcript_rows['feature'] == '3utr') | (transcript_rows['feature'] == 'three_prime_utr')]

            # Add exons to the transcript
            for _, exon in exons.iterrows():
                exon_sequence = get_sequence(genome_seq, chromosome, exon['start'] - 1, exon['end'], strand)
                exon_obj = Exon(transcript_id, exon_sequence, chromosome, (exon['start'], exon['end']), strand)
                transcript.exons.append(exon_obj)
            
            # Add CDS to the transcript
            for _, cd in cds.iterrows():
                cds_sequence = get_sequence(genome_seq, chromosome, cd['start'] - 1, cd['end'], strand)
                cds_obj = CDS(transcript_id, cds_sequence, chromosome, (cd['start'], cd['end']), strand)
                transcript.cds.append(cds_obj)
            
            # Add 5' UTR to the transcript
            for _, utr in utr5.iterrows():
                utr_sequence = get_sequence(genome_seq, chromosome, utr['start'] - 1, utr['end'], strand)
                utr_obj = UTR(transcript_id, utr_sequence, chromosome, (utr['start'], utr['end']), strand, '5UTR')
                transcript.utrs.append(utr_obj)

            # Add 3' UTR to the transcript
            for _, utr in utr3.iterrows():
                utr_sequence = get_sequence(genome_seq, chromosome, utr['start'] - 1, utr['end'], strand)
                utr_obj = UTR(transcript_id, utr_sequence, chromosome, (utr['start'], utr['end']), strand, '3UTR')
                transcript.utrs.append(utr_obj)

            # Define introns
            exons = sorted(transcript.exons, key=lambda x: x.position[0], reverse=False)
            for i in range(len(exons) - 1):
                intron_start = exons[i].position[1] + 1  # end of previous exon
                intron_end = exons[i+1].position[0] - 1  # start of next exon
                if intron_end - intron_start > 1:  # only add introns with length > 1
                    intron_sequence = get_sequence(genome_seq, chromosome, intron_start - 1, intron_end, strand)
                    intron_obj = Intron(transcript_id, intron_sequence, chromosome, (intron_start, intron_end), strand)
                    transcript.introns.append(intron_obj)

            # Sort exons, UTRs, CDS, and introns
            reverse_sort = True if strand == '-' else False # reverse sort if strand is - else sort normally if strand is +
            transcript.utrs = sorted(transcript.utrs, key=lambda x: x.position[0], reverse=reverse_sort)
            transcript.exons = sorted(transcript.exons, key=lambda x: x.position[0], reverse=reverse_sort)
            transcript.cds = sorted(transcript.cds, key=lambda x: x.position[0], reverse=reverse_sort)
            transcript.introns = sorted(transcript.introns, key=lambda x: x.position[0], reverse=reverse_sort)
            
            # Get CDS sequence
            transcript.cds_sequence = ''.join([cds.sequence for cds in transcript.cds])
            transcript.cds_length = len(transcript.cds_sequence)

            # Get mRNA sequence 
            transcript.mrna_sequence = ''.join([exon.sequence for exon in transcript.exons])

            # Get DNA sequence (exon and intron sequence)
            combined = transcript.exons[:] + transcript.introns[:] 
            combined = [x for x in combined if x is not None]
            combined = sorted(combined, key=lambda x: x.position[0], reverse=reverse_sort) 
            transcript.dna_sequence = ''.join([feature.sequence for feature in combined]) 

            # Add position to transcript
            transcript.position = transcript.get_bounds() 

            # Add the transcript to the gene
            gene.add_transcript(transcript)

        # Add the gene to the transcriptome
        transcriptome.add_gene(gene)

    # The transcriptome object is now populated with genes, transcripts, exons, etc.
    print(transcriptome)

    return transcriptome 