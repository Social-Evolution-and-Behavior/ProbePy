import pickle
import os
from hcrfish.transcriptomics.main import generate_transcriptome_object

def update_transcriptome_object(
        genome_path, 
        transcriptome_path, 
        species_identifier, 
        base_dir="",
        overwrite=False
):
    """ Creates the transcriptome object and saves it to a file."""    

    # Check that both paths exist
    if not os.path.exists(genome_path):
        raise FileNotFoundError(f"Genome path {genome_path} does not exist.")
    if not os.path.exists(transcriptome_path):
        raise FileNotFoundError(f"Transcriptome path {transcriptome_path} does not exist.")

    # Create output path 
    output_path = os.path.join(base_dir, "input", species_identifier, f"{species_identifier}_transcriptome.pkl")

    # If the output file exists and overwrite is False, raise an error
    if os.path.exists(output_path) and not overwrite:
        print(f"File {output_path} already exists. Set overwrite=True to overwrite it.")
        return  # Stop the update if we don't want to overwrite

    # Generate the transcriptome object
    transcriptome_obj = generate_transcriptome_object(transcriptome_path, genome_path)

    # Serialize and save the object to a file
    with open(output_path, 'wb') as f:
        pickle.dump(transcriptome_obj, f)
    
    print(f"Transcriptome object has been updated and saved to {output_path}")


# Load the transcriptome object from a file
def load_transcriptome_object(species_identifier, base_dir=""):
    """Load the transcriptome object from a file."""
    input_path = os.path.join(base_dir, "input", species_identifier, f"{species_identifier}_transcriptome.pkl")

    try:
        with open(input_path, 'rb') as f:
            transcriptome_obj = pickle.load(f)
        print(f"Loaded transcriptome object from {input_path}")
        return transcriptome_obj
    except FileNotFoundError:
        print(f"File {input_path} not found.")
        print("Please run update_transcriptome_object() to generate the transcriptome object.")
        return None


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
