"""
Tests for probepy.transcriptomics.main module.

This module tests the generate_transcriptome_object function which is the
main function for creating complete Transcriptome objects from GTF and genome files.
"""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, Mock, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO

from probepy.transcriptomics.main import generate_transcriptome_object
from probepy.transcriptomics.classes import Transcriptome, Gene, Transcript


class TestGenerateTranscriptomeObject:
    """Test generate_transcriptome_object function."""
    
    @pytest.fixture
    def simple_genome_file(self):
        """Create a simple genome FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
            f.write(">chr2\n")
            f.write("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n")
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except:
            pass
    
    @pytest.fixture
    def simple_gtf_file(self):
        """Create a simple GTF file with one gene."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            # Gene with one transcript and one exon
            f.write('chr1\tHAVANA\tgene\t10\t30\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t10\t30\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tCDS\t12\t18\t.\t+\t0\tgene_id "GENE001"; gene_name "TestGene"; transcript_id "TRANS001";\n')
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except:
            pass
    
    @pytest.fixture
    def multi_exon_gtf_file(self):
        """Create GTF file with multi-exon gene."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            # Gene with one transcript and multiple exons
            f.write('chr1\tHAVANA\tgene\t10\t50\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t10\t50\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t30\t40\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t45\t50\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tCDS\t12\t18\t.\t+\t0\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tCDS\t30\t38\t.\t+\t0\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except:
            pass
    
    @pytest.fixture
    def gtf_with_utrs(self):
        """Create GTF file with UTR features."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t50\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t10\t50\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t10\t50\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\t5UTR\t10\t15\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tCDS\t16\t40\t.\t+\t0\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\t3UTR\t41\t50\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except:
            pass
    
    @pytest.fixture
    def reverse_strand_gtf_file(self):
        """Create GTF file with gene on reverse strand."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t30\t.\t-\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t-\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tCDS\t12\t18\t.\t-\t0\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except:
            pass
    
    def test_basic_transcriptome_creation(self, simple_genome_file, simple_gtf_file):
        """Test basic transcriptome object creation."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        assert isinstance(transcriptome, Transcriptome)
        assert len(transcriptome.genes) > 0
    
    def test_gene_creation(self, simple_genome_file, simple_gtf_file):
        """Test that genes are created correctly."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        # Should have one gene
        assert len(transcriptome.genes) == 1
        
        # Check gene details
        gene = transcriptome.get_gene("TestGene")
        assert gene is not None
        assert gene.name == "TestGene"
        assert gene.id == "GENE001"
        assert gene.chromosome == "chr1"
        assert gene.strand == "+"
    
    def test_transcript_creation(self, simple_genome_file, simple_gtf_file):
        """Test that transcripts are created correctly."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        assert len(gene.transcripts) == 1
        
        transcript = gene.transcripts[0]
        assert transcript.name == "TRANS001"
        assert transcript.chromosome == "chr1"
        assert transcript.strand == "+"
    
    def test_exon_creation(self, simple_genome_file, simple_gtf_file):
        """Test that exons are created with sequences."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert len(transcript.exons) > 0
        
        exon = transcript.exons[0]
        assert exon.sequence is not None
        assert len(exon.sequence) > 0
        assert isinstance(exon.sequence, str)
        assert exon.chromosome == "chr1"
        assert exon.strand == "+"
    
    def test_cds_creation(self, simple_genome_file, simple_gtf_file):
        """Test that CDS features are created."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert len(transcript.cds) > 0
        
        cds = transcript.cds[0]
        assert cds.sequence is not None
        assert len(cds.sequence) > 0
        assert isinstance(cds.sequence, str)
    
    def test_mrna_sequence_generation(self, simple_genome_file, simple_gtf_file):
        """Test that mRNA sequence is generated from exons."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert transcript.mrna_sequence is not None
        assert len(transcript.mrna_sequence) > 0
        assert isinstance(transcript.mrna_sequence, str)
        
        # mRNA should be concatenation of exons
        exon_concat = ''.join([exon.sequence for exon in transcript.exons])
        assert transcript.mrna_sequence == exon_concat
    
    def test_dna_sequence_generation(self, simple_genome_file, multi_exon_gtf_file):
        """Test that DNA sequence includes introns."""
        transcriptome = generate_transcriptome_object(multi_exon_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert transcript.dna_sequence is not None
        assert len(transcript.dna_sequence) > 0
        
        # DNA should be longer than mRNA if introns exist
        if len(transcript.introns) > 0:
            assert len(transcript.dna_sequence) >= len(transcript.mrna_sequence)
    
    def test_intron_detection(self, simple_genome_file, multi_exon_gtf_file):
        """Test that introns are detected between exons."""
        transcriptome = generate_transcriptome_object(multi_exon_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        # Should detect introns between the 3 exons
        assert len(transcript.introns) == 2
        
        # Each intron should have sequence
        for intron in transcript.introns:
            assert intron.sequence is not None
            assert len(intron.sequence) > 0
    
    def test_utr_creation(self, simple_genome_file, gtf_with_utrs):
        """Test that UTR features are created."""
        transcriptome = generate_transcriptome_object(gtf_with_utrs, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert len(transcript.utrs) > 0
        
        # Should have both 5' and 3' UTRs
        utr_types = [utr.utr_type for utr in transcript.utrs]
        assert '5UTR' in utr_types or 'five_prime_utr' in utr_types
        assert '3UTR' in utr_types or 'three_prime_utr' in utr_types
    
    def test_reverse_strand_handling(self, simple_genome_file, reverse_strand_gtf_file):
        """Test correct handling of reverse strand genes."""
        transcriptome = generate_transcriptome_object(reverse_strand_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        assert gene.strand == "-"
        
        transcript = gene.transcripts[0]
        assert transcript.strand == "-"
        
        # Sequences should be reverse complemented
        assert transcript.mrna_sequence is not None
        assert len(transcript.mrna_sequence) > 0
    
    def test_feature_sorting(self, simple_genome_file, multi_exon_gtf_file):
        """Test that features are sorted correctly."""
        transcriptome = generate_transcriptome_object(multi_exon_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        # Exons should be sorted
        if len(transcript.exons) > 1:
            for i in range(len(transcript.exons) - 1):
                # For forward strand, positions should increase
                if transcript.strand == "+":
                    assert transcript.exons[i].position[0] <= transcript.exons[i+1].position[0]
    
    def test_cds_length_calculation(self, simple_genome_file, simple_gtf_file):
        """Test that CDS length is calculated correctly."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert transcript.cds_length > 0
        assert transcript.cds_length == len(transcript.cds_sequence)
    
    def test_transcript_position_assignment(self, simple_genome_file, simple_gtf_file):
        """Test that transcript positions are calculated."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        gene = transcriptome.get_gene("TestGene")
        transcript = gene.transcripts[0]
        
        assert transcript.position is not None
        assert len(transcript.position) == 2
        assert transcript.position[0] < transcript.position[1]
    
    def test_invalid_chromosome_skipped(self, simple_genome_file):
        """Test that genes on chromosomes not in genome are skipped."""
        # Create GTF with invalid chromosome
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr99\tHAVANA\tgene\t10\t30\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr99\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            
            # Gene should be skipped
            assert len(transcriptome.genes) == 0
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_invalid_strand_skipped(self, simple_genome_file):
        """Test that genes with invalid strand are skipped."""
        # Create GTF with invalid strand
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t30\t.\t.\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t.\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            
            # Gene should be skipped
            assert len(transcriptome.genes) == 0
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_multiple_genes(self, simple_genome_file):
        """Test processing multiple genes."""
        # Create GTF with multiple genes
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "GENE001"; gene_name "Gene1";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\tgene\t30\t40\t.\t+\t.\tgene_id "GENE002"; gene_name "Gene2";\n')
            f.write('chr1\tHAVANA\texon\t30\t40\t.\t+\t.\tgene_id "GENE002"; transcript_id "TRANS002";\n')
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            
            assert len(transcriptome.genes) == 2
            assert transcriptome.get_gene("Gene1") is not None
            assert transcriptome.get_gene("Gene2") is not None
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_multiple_transcripts_per_gene(self, simple_genome_file):
        """Test gene with multiple transcripts."""
        # Create GTF with one gene, two transcripts
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t40\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            f.write('chr1\tHAVANA\ttranscript\t10\t30\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";\n')
            f.write('chr1\tHAVANA\ttranscript\t15\t40\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS002";\n')
            f.write('chr1\tHAVANA\texon\t15\t25\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS002";\n')
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            
            gene = transcriptome.get_gene("TestGene")
            assert len(gene.transcripts) == 2
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_string_representation(self, simple_genome_file, simple_gtf_file):
        """Test that transcriptome has proper string representation."""
        transcriptome = generate_transcriptome_object(simple_gtf_file, simple_genome_file)
        
        trans_str = str(transcriptome)
        assert "Transcriptome" in trans_str
        assert "genes=" in trans_str


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_gtf_file(self, simple_genome_file):
        """Test with empty GTF file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write("# Empty GTF file\n")
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            assert len(transcriptome.genes) == 0
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_gene_without_transcripts(self, simple_genome_file):
        """Test gene entry without transcript entries."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\tHAVANA\tgene\t10\t30\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";\n')
            # No transcript or exon entries
            gtf_path = f.name
        
        try:
            transcriptome = generate_transcriptome_object(gtf_path, simple_genome_file)
            # Gene might be created but should have no transcripts
            if "TestGene" in transcriptome.genes:
                gene = transcriptome.get_gene("TestGene")
                assert len(gene.transcripts) == 0
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
