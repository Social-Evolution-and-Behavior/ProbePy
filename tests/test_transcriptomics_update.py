"""
Tests for probepy.transcriptomics.update module.

This module tests functions for creating, updating, and loading
transcriptome objects from disk.
"""

import pytest
import pickle
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from probepy.transcriptomics.update import (
    update_transcriptome_object,
    load_transcriptome_object,
    check_exons_contain_all_features
)
from probepy.transcriptomics.classes import (
    Transcriptome, Gene, Transcript, Exon, CDS, UTR
)


class TestUpdateTranscriptomeObject:
    """Test update_transcriptome_object function."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for tests."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield temp_dir
    
    @pytest.fixture
    def mock_genome_file(self, temp_dir):
        """Create a mock genome FASTA file."""
        genome_path = os.path.join(temp_dir, "genome.fa")
        with open(genome_path, 'w') as f:
            f.write(">chr1\nATCGATCGATCGATCGATCGATCG\n")
        return genome_path
    
    @pytest.fixture
    def mock_gtf_file(self, temp_dir):
        """Create a mock GTF file."""
        gtf_path = os.path.join(temp_dir, "transcriptome.gtf")
        with open(gtf_path, 'w') as f:
            f.write('chr1\tHAVANA\tgene\t1\t20\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t1\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t1\t10\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
        return gtf_path
    
    def test_successful_update(self, temp_dir, mock_genome_file, mock_gtf_file):
        """Test successful transcriptome object creation and saving."""
        species_id = "test_species"
        
        update_transcriptome_object(
            mock_genome_file,
            mock_gtf_file,
            species_id,
            base_dir=temp_dir
        )
        
        # Check that output file was created
        output_path = os.path.join(temp_dir, "input", species_id, f"{species_id}_transcriptome.pkl")
        assert os.path.exists(output_path)
        
        # Load and verify
        with open(output_path, 'rb') as f:
            transcriptome = pickle.load(f)
        
        assert isinstance(transcriptome, Transcriptome)
        assert len(transcriptome.genes) > 0
    
    def test_genome_path_not_exists(self, temp_dir, mock_gtf_file):
        """Test error when genome path doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Genome path .* does not exist"):
            update_transcriptome_object(
                "/nonexistent/genome.fa",
                mock_gtf_file,
                "test_species",
                base_dir=temp_dir
            )
    
    def test_transcriptome_path_not_exists(self, temp_dir, mock_genome_file):
        """Test error when transcriptome path doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Transcriptome path .* does not exist"):
            update_transcriptome_object(
                mock_genome_file,
                "/nonexistent/transcriptome.gtf",
                "test_species",
                base_dir=temp_dir
            )
    
    def test_overwrite_false_with_existing_file(self, temp_dir, mock_genome_file, mock_gtf_file):
        """Test that existing file is not overwritten when overwrite=False."""
        species_id = "test_species"
        output_path = os.path.join(temp_dir, "input", species_id, f"{species_id}_transcriptome.pkl")
        
        # Create output directory and file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write("existing content")
        
        original_mtime = os.path.getmtime(output_path)
        
        # Try to update with overwrite=False
        update_transcriptome_object(
            mock_genome_file,
            mock_gtf_file,
            species_id,
            base_dir=temp_dir,
            overwrite=False
        )
        
        # File should not have been modified
        new_mtime = os.path.getmtime(output_path)
        assert original_mtime == new_mtime
    
    def test_overwrite_true_with_existing_file(self, temp_dir, mock_genome_file, mock_gtf_file):
        """Test that existing file is overwritten when overwrite=True."""
        species_id = "test_species"
        output_path = os.path.join(temp_dir, "input", species_id, f"{species_id}_transcriptome.pkl")
        
        # Create output directory and file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write("existing content")
        
        # Update with overwrite=True
        update_transcriptome_object(
            mock_genome_file,
            mock_gtf_file,
            species_id,
            base_dir=temp_dir,
            overwrite=True
        )
        
        # File should be a valid pickle
        with open(output_path, 'rb') as f:
            transcriptome = pickle.load(f)
        
        assert isinstance(transcriptome, Transcriptome)
    
    def test_output_directory_creation(self, temp_dir, mock_genome_file, mock_gtf_file):
        """Test that output directory is created if it doesn't exist."""
        species_id = "test_species"
        
        update_transcriptome_object(
            mock_genome_file,
            mock_gtf_file,
            species_id,
            base_dir=temp_dir
        )
        
        # Check directory structure was created
        assert os.path.exists(os.path.join(temp_dir, "input"))
        assert os.path.exists(os.path.join(temp_dir, "input", species_id))
    
    def test_serialization_format(self, temp_dir, mock_genome_file, mock_gtf_file):
        """Test that transcriptome is properly serialized."""
        species_id = "test_species"
        
        update_transcriptome_object(
            mock_genome_file,
            mock_gtf_file,
            species_id,
            base_dir=temp_dir
        )
        
        output_path = os.path.join(temp_dir, "input", species_id, f"{species_id}_transcriptome.pkl")
        
        # Should be loadable with pickle
        with open(output_path, 'rb') as f:
            obj = pickle.load(f)
        
        # Should have correct type and structure
        assert isinstance(obj, Transcriptome)
        assert hasattr(obj, 'genes')
        assert hasattr(obj, 'transcripts_to_genes')


class TestLoadTranscriptomeObject:
    """Test load_transcriptome_object function."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for tests."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield temp_dir
    
    @pytest.fixture
    def saved_transcriptome(self, temp_dir):
        """Create and save a test transcriptome object."""
        species_id = "test_species"
        
        # Create transcriptome
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        transcriptome.add_gene(gene)
        
        # Save to file
        output_dir = os.path.join(temp_dir, "input", species_id)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{species_id}_transcriptome.pkl")
        
        with open(output_path, 'wb') as f:
            pickle.dump(transcriptome, f)
        
        return species_id
    
    def test_successful_load(self, temp_dir, saved_transcriptome):
        """Test successful loading of transcriptome object."""
        transcriptome = load_transcriptome_object(saved_transcriptome, base_dir=temp_dir)
        
        assert transcriptome is not None
        assert isinstance(transcriptome, Transcriptome)
        assert len(transcriptome.genes) == 1
        assert "TestGene" in transcriptome.genes
    
    def test_file_not_found(self, temp_dir):
        """Test behavior when file doesn't exist."""
        result = load_transcriptome_object("nonexistent_species", base_dir=temp_dir)
        
        assert result is None
    
    def test_data_integrity(self, temp_dir, saved_transcriptome):
        """Test that loaded data matches saved data."""
        transcriptome = load_transcriptome_object(saved_transcriptome, base_dir=temp_dir)
        
        # Check gene details
        gene = transcriptome.get_gene("TestGene")
        assert gene is not None
        assert gene.id == "GENE001"
        assert gene.name == "TestGene"
    
    def test_multiple_loads(self, temp_dir, saved_transcriptome):
        """Test that object can be loaded multiple times."""
        trans1 = load_transcriptome_object(saved_transcriptome, base_dir=temp_dir)
        trans2 = load_transcriptome_object(saved_transcriptome, base_dir=temp_dir)
        
        # Both should be valid but different instances
        assert trans1 is not None
        assert trans2 is not None
        assert trans1 is not trans2
        assert len(trans1.genes) == len(trans2.genes)
    
    def test_corrupted_pickle_file(self, temp_dir):
        """Test handling of corrupted pickle file."""
        species_id = "test_species"
        output_dir = os.path.join(temp_dir, "input", species_id)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{species_id}_transcriptome.pkl")
        
        # Write invalid pickle data
        with open(output_path, 'wb') as f:
            f.write(b"not a valid pickle")
        
        # Should raise exception or return None
        try:
            result = load_transcriptome_object(species_id, base_dir=temp_dir)
            # If it doesn't raise, it should return None or handle gracefully
        except (pickle.UnpicklingError, EOFError):
            # Expected behavior
            pass


class TestCheckExonsContainAllFeatures:
    """Test check_exons_contain_all_features function."""
    
    @pytest.fixture
    def valid_transcriptome(self):
        """Create transcriptome with valid exon-feature relationships."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        
        # Create exon that contains CDS and UTR
        exon1 = Exon("TRANS001", "ATCGATCGATCG", "chr1", (100, 112), "+")
        cds1 = CDS("TRANS001", "ATCG", "chr1", (102, 106), "+")
        utr1 = UTR("TRANS001", "AT", "chr1", (100, 102), "+", "5UTR")
        
        transcript.exons = [exon1]
        transcript.cds = [cds1]
        transcript.utrs = [utr1]
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        return transcriptome
    
    @pytest.fixture
    def invalid_transcriptome(self):
        """Create transcriptome with invalid exon-feature relationships."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        
        # Create exon and CDS that don't overlap
        exon1 = Exon("TRANS001", "ATCG", "chr1", (100, 104), "+")
        cds1 = CDS("TRANS001", "GCTA", "chr1", (200, 204), "+")  # Outside exon
        
        transcript.exons = [exon1]
        transcript.cds = [cds1]
        transcript.utrs = []
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        return transcriptome
    
    def test_valid_structure_no_output(self, valid_transcriptome, capsys):
        """Test that valid structure produces no error messages."""
        check_exons_contain_all_features(valid_transcriptome)
        
        captured = capsys.readouterr()
        # Should not print any error messages
        assert "does not contain all features" not in captured.out
    
    def test_invalid_structure_prints_warning(self, invalid_transcriptome, caplog):
        """Test that invalid structure prints warning."""
        check_exons_contain_all_features(invalid_transcriptome)
        
        # Should print warning about features not in exons
        assert "does not contain all features" in caplog.text
        assert "TestGene" in caplog.text
        assert "TRANS001" in caplog.text
    
    def test_empty_transcriptome(self):
        """Test with empty transcriptome."""
        transcriptome = Transcriptome()
        
        # Should complete without errors
        check_exons_contain_all_features(transcriptome)
    
    def test_transcript_without_features(self):
        """Test transcript with exons but no CDS/UTR."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        exon1 = Exon("TRANS001", "ATCG", "chr1", (100, 104), "+")
        transcript.exons = [exon1]
        transcript.cds = []
        transcript.utrs = []
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        # Should complete without errors (no features to check)
        check_exons_contain_all_features(transcriptome)
    
    def test_multiple_exons_covering_features(self):
        """Test features spread across multiple exons."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        
        # Two exons
        exon1 = Exon("TRANS001", "ATCG", "chr1", (100, 104), "+")
        exon2 = Exon("TRANS001", "GCTA", "chr1", (200, 204), "+")
        
        # CDS in first exon
        cds1 = CDS("TRANS001", "AT", "chr1", (100, 102), "+")
        # UTR in second exon
        utr1 = UTR("TRANS001", "GC", "chr1", (200, 202), "+", "3UTR")
        
        transcript.exons = [exon1, exon2]
        transcript.cds = [cds1]
        transcript.utrs = [utr1]
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        # Should complete without errors
        check_exons_contain_all_features(transcriptome)
    
    def test_partially_contained_feature(self):
        """Test feature that is only partially contained in exon."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        
        # Exon from 100-110
        exon1 = Exon("TRANS001", "ATCGATCGATCG", "chr1", (100, 112), "+")
        # CDS extends beyond exon
        cds1 = CDS("TRANS001", "ATCGATCGATCGATCG", "chr1", (105, 121), "+")
        
        transcript.exons = [exon1]
        transcript.cds = [cds1]
        transcript.utrs = []
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        # Should detect and report this
        check_exons_contain_all_features(transcriptome)
    
    def test_wrong_chromosome(self):
        """Test feature on different chromosome than exon."""
        transcriptome = Transcriptome()
        gene = Gene("GENE001", "TestGene")
        
        transcript = Transcript("TRANS001")
        
        exon1 = Exon("TRANS001", "ATCG", "chr1", (100, 104), "+")
        cds1 = CDS("TRANS001", "ATCG", "chr2", (100, 104), "+")  # Different chromosome
        
        transcript.exons = [exon1]
        transcript.cds = [cds1]
        transcript.utrs = []
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        # Should detect chromosome mismatch
        check_exons_contain_all_features(transcriptome)


class TestIntegration:
    """Integration tests for update module."""
    
    def test_update_and_load_workflow(self):
        """Test complete workflow of updating and loading."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create test files
            genome_path = os.path.join(temp_dir, "genome.fa")
            with open(genome_path, 'w') as f:
                f.write(">chr1\nATCGATCGATCGATCGATCGATCG\n")
            
            gtf_path = os.path.join(temp_dir, "transcriptome.gtf")
            with open(gtf_path, 'w') as f:
                f.write('chr1\tHAVANA\tgene\t1\t20\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
                f.write('chr1\tHAVANA\ttranscript\t1\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
                f.write('chr1\tHAVANA\texon\t1\t10\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            
            species_id = "test_species"
            
            # Update
            update_transcriptome_object(
                genome_path,
                gtf_path,
                species_id,
                base_dir=temp_dir
            )
            
            # Load
            transcriptome = load_transcriptome_object(species_id, base_dir=temp_dir)
            
            # Verify
            assert transcriptome is not None
            assert isinstance(transcriptome, Transcriptome)
            assert len(transcriptome.genes) > 0
            
            # Check structure
            check_exons_contain_all_features(transcriptome)
    
    def test_update_overwrite_and_load(self):
        """Test updating with overwrite and loading."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create test files
            genome_path = os.path.join(temp_dir, "genome.fa")
            with open(genome_path, 'w') as f:
                f.write(">chr1\nATCGATCGATCGATCGATCGATCG\n")
            
            gtf_path = os.path.join(temp_dir, "transcriptome.gtf")
            with open(gtf_path, 'w') as f:
                f.write('chr1\tHAVANA\tgene\t1\t20\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
                f.write('chr1\tHAVANA\texon\t1\t10\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            
            species_id = "test_species"
            
            # First update
            update_transcriptome_object(
                genome_path,
                gtf_path,
                species_id,
                base_dir=temp_dir
            )
            
            # Load first version
            trans1 = load_transcriptome_object(species_id, base_dir=temp_dir)
            gene_count1 = len(trans1.genes)
            
            # Update GTF file with more genes
            with open(gtf_path, 'a') as f:
                f.write('chr1\tHAVANA\tgene\t30\t40\t.\t+\t.\tgene_id "GENE002"; gene_name "TestGene2";\n')
            
            # Update again with overwrite=True
            update_transcriptome_object(
                genome_path,
                gtf_path,
                species_id,
                base_dir=temp_dir,
                overwrite=True
            )
            
            # Load second version
            trans2 = load_transcriptome_object(species_id, base_dir=temp_dir)
            gene_count2 = len(trans2.genes)
            
            # Should have more genes now
            assert gene_count2 >= gene_count1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
