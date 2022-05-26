"""Test the code used for generating pham MSAs with Clustal Omega.

Inherently linked to the success of unit tests for sub-processing."""

import pathlib
import unittest

from phamerate.clustalo import *

clustalo_version = get_clustalo_version()
if not clustalo_version:
    print("ClustalO not found - skipping unittests in 'test_clustalo'")
    exit()

file_dir = pathlib.Path(__file__).parent.parent.joinpath("files/clustalo")

large_fasta = file_dir.joinpath("large_pham.faa")
small_fasta = file_dir.joinpath("small_pham.faa")
orpham_fasta = file_dir.joinpath("orpham.faa")
not_a_file = file_dir.joinpath("not_a_file.faa")
outfile = file_dir.joinpath("test_clustalo.aln")


class TestClustalOmega(unittest.TestCase):
    def test_clustalo_large(self):
        """Large phams have many sequences - test that this works."""
        _, o = run_clustalo(large_fasta, outfile)
        self.assertTrue(o.is_file())
        o.unlink()

    def test_clustalo_small(self):
        """Small phams have a few sequences - test that this works."""
        _, o = run_clustalo(small_fasta, outfile)
        self.assertTrue(o.is_file())
        o.unlink()

    def test_clustalo_orpham(self):
        """An orpham is not alignable and should raise an error."""
        with self.assertRaises(ClustalOmegaError):
            run_clustalo(orpham_fasta, outfile)

    def test_clustalo_fake_file(self):
        """A fake file should raise an error."""
        with self.assertRaises(ClustalOmegaError):
            run_clustalo(not_a_file, outfile)


if __name__ == '__main__':
    unittest.main()
