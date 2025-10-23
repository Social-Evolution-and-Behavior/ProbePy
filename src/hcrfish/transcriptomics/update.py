import pickle
from hcrfish.transcriptomics.main import generate_transcriptome_object

def update_transcriptome_object(genome_path, transcriptome_path, output_filename, species):
    """ Creates the transcriptome object and saves it to a file."""    

    # Generate the transcriptome object
    transcriptome_obj = generate_transcriptome_object(transcriptome_path, genome_path, species)
    
    # Add .pkl to output_filename if it's not already there
    if not output_filename.endswith('.pkl'):
        output_filename += '.pkl'
    
    # Serialize and save the object to a file
    with open(output_filename, 'wb') as f:
        pickle.dump(transcriptome_obj, f)
    
    print(f"Transcriptome object has been updated and saved to {output_filename}")


# Load the transcriptome object from a file
def load_transcriptome_object(filename):
    # Add .pkl to output_filename if it's not already there
    if not filename.endswith('.pkl'):
        filename += '.pkl'
    try:
        with open(filename, 'rb') as f:
            transcriptome_obj = pickle.load(f)
    except FileNotFoundError:
        print(f"File {filename} not found. Please run update(genome_path, transcriptome_path, output_filename) to generate the transcriptome object.")
        return None

    return transcriptome_obj


def check_exons_contain_all_features(transcriptome_obj): 
    """Check if all features (CDS and UTRs) are contained within exons."""
    def is_contained(exon, feature):
        """Check if a feature is contained within an exon."""
        return (
            exon.chromosome == feature.chromosome and
            exon.position[0] <= feature.position[0] and
            exon.position[1] >= feature.position[1]
        )

    def verify_exons_contain_features(exons, features):
        """Verify that all features are contained within exons."""
        for feature in features:
            if not any(is_contained(exon, feature) for exon in exons):
                return False  # If any feature isn't contained, return False
        return True  # All features are contained
    
    for gene_name in transcriptome_obj.genes.keys():
        gene = transcriptome_obj.get_gene(gene_name) 
        for transcript in gene.transcripts: 
            exons = transcript.exons
            features = transcript.cds + transcript.utrs
            contains_all_features = verify_exons_contain_features(exons, features)
            if not contains_all_features: 
                print(f"Gene {gene.name} transcript {transcript.name} does not contain all features in exons")
