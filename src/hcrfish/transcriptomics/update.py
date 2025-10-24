import pickle
import os
from hcrfish.transcriptomics.main import generate_transcriptome_object

def update_transcriptome_object(genome_path, transcriptome_path, output_filename, species):
    """ Creates the transcriptome object and saves it to a file."""    

    # Generate the transcriptome object
    transcriptome_obj = generate_transcriptome_object(transcriptome_path, genome_path, species)
    
    # Add .pkl to output_filename if it's not already there
    if not output_filename.endswith('.pkl'):
        output_filename += '.pkl'
    
    # Create the species-specific directory structure
    species_dir = f"input/{species}"
    os.makedirs(species_dir, exist_ok=True)
    
    # Create the full path for the output file
    output_path = os.path.join(species_dir, output_filename)
    
    # Serialize and save the object to a file
    with open(output_path, 'wb') as f:
        pickle.dump(transcriptome_obj, f)
    
    print(f"Transcriptome object has been updated and saved to {output_path}")


# Load the transcriptome object from a file
def load_transcriptome_object(filename, species=None):
    # Add .pkl to output_filename if it's not already there
    if not filename.endswith('.pkl'):
        filename += '.pkl'
    
    # Try to find the file in species directory first, then fallback to current directory
    file_paths = []
    
    # If species is provided, look in the species directory first
    if species:
        species_path = os.path.join(f"input/{species}", filename)
        file_paths.append(species_path)
    
    # Also try to find in common species directories
    for common_species in ['dmel', 'dyak']:
        species_path = os.path.join(f"input/{common_species}", filename)
        if species_path not in file_paths:
            file_paths.append(species_path)
    
    # Finally, try the current directory and docs directory for backward compatibility
    file_paths.extend([filename, os.path.join("docs", filename)])
    
    # Try each path until we find the file
    for filepath in file_paths:
        try:
            with open(filepath, 'rb') as f:
                transcriptome_obj = pickle.load(f)
            print(f"Loaded transcriptome object from {filepath}")
            return transcriptome_obj
        except FileNotFoundError:
            continue
    
    # If no file found, print error message
    print(f"File {filename} not found in any of the expected locations:")
    for path in file_paths:
        print(f"  - {path}")
    print("Please run update_transcriptome_object(genome_path, transcriptome_path, output_filename, species) to generate the transcriptome object.")
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
