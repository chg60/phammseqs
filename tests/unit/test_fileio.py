"""Test the code used for FASTA reading and writing, and Genbank
flatfile reading."""

import pathlib
import unittest

from phamerate.fileio import *

file_dir = pathlib.Path(__file__).parent.parent.joinpath("files/fileio")

# Genome files to test
charlie_gff3 = file_dir.joinpath("Charlie.gff3")
mmb_no_wrap = file_dir.joinpath("MMB_no_wrap.faa")
mmb_wrap = file_dir.joinpath("MMB_wrap.faa")
phrann_gb = file_dir.joinpath("Phrann.gb")
phrann_fasta = file_dir.joinpath("Phrann.faa")


class TestSniffFormat(unittest.TestCase):
    def test_sniff_wrapped_fasta(self):
        """Test that sniff_format() returns "fasta" if wrapped.
        """
        self.assertEqual(sniff_format(mmb_wrap), "fasta")

    def test_sniff_unwrapped_fasta(self):
        """Test that sniff_format() returns "fasta" if not wrapped.
        """
        self.assertEqual(sniff_format(mmb_no_wrap), "fasta")

    def test_sniff_genbank(self):
        """Test that sniff_format() returns "genbank" for flatfile.
        """
        self.assertEqual(sniff_format(phrann_gb), "genbank")

    def test_no_sniff_gff(self):
        """Test that sniff_format() returns None for GFF3 file.
        """
        self.assertIsNone(sniff_format(charlie_gff3))


class TestFastaRead(unittest.TestCase):
    def test_read_wrapped_fasta(self):
        """Test that fasta.read_fasta() can parse a line-wrapped FASTA
        file.
        """
        # Check that 70 sequences are parsed
        wrap_seqs = read_fasta(mmb_wrap)
        self.assertEqual(len([x for x in wrap_seqs]), 70)

    def test_read_unwrapped_fasta(self):
        """Test that fasta.read_fasta() can parse non-wrapped FASTA.
        """
        # Check that 70 sequences are parsed
        unwrap_seqs = read_fasta(mmb_no_wrap)
        self.assertEqual(len([x for x in unwrap_seqs]), 70)

    def test_wrapped_unwrapped_equal(self):
        """Test that fasta.read_fasta() gets the same contents from FASTA
        files that differ only in whether they are line-wrapped.
        """
        # The two files differ ONLY in that one has sequences line-wrapped
        wrap_seqs = read_fasta(mmb_wrap)
        unwrap_seqs = read_fasta(mmb_no_wrap)

        for wrap_seq, unwrap_seq in zip(wrap_seqs, unwrap_seqs):
            # Check that headers are equal
            self.assertEqual(wrap_seq[0], unwrap_seq[0])
            # Check that translations are equal
            self.assertEqual(wrap_seq[1], unwrap_seq[1])

    def test_no_read_genbank(self):
        """Test that fasta.read_fasta() does not parse anything from a
        Genbank flatfile."""
        # Check that 0 of the 69 sequences are parsed
        gbk_seqs = read_fasta(phrann_gb)
        self.assertEqual(len([x for x in gbk_seqs]), 0)


class TestFastaReadWrite(unittest.TestCase):
    def setUp(self):
        self.sequences = read_genbank(phrann_gb)

    def test_unwrapped_round_trip(self):
        write_fasta(self.sequences, phrann_fasta, wrap=0)
        unwrap_seqs = read_fasta(phrann_fasta)

        for seq, unwrap_seq in zip(self.sequences, unwrap_seqs):
            # Check that headers are equal
            self.assertEqual(seq[0], unwrap_seq[0])
            # Check that translations are equal
            self.assertEqual(seq[1], unwrap_seq[1])

    def test_wrapped_round_trip(self):
        write_fasta(self.sequences, phrann_fasta)
        wrap_seqs = read_fasta(phrann_fasta)

        for seq, wrap_seq in zip(self.sequences, wrap_seqs):
            # Check that headers are equal
            self.assertEqual(seq[0], wrap_seq[0])
            # Check that translations are equal
            self.assertEqual(seq[1], wrap_seq[1])


class TestGenbankRead(unittest.TestCase):
    def test_read_flatfile(self):
        # Check that 67 sequences are parsed
        gbk_seqs = read_genbank(phrann_gb)
        self.assertEqual(len([x for x in gbk_seqs]), 67)

    def test_no_read_fasta(self):
        unwrap_seqs = read_genbank(mmb_no_wrap)
        self.assertEqual(len([x for x in unwrap_seqs]), 0)


if __name__ == '__main__':
    unittest.main()
