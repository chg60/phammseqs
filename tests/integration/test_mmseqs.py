"""Test the code used to perform pham assembly with MMseqs2.

Inherently linked to the success of unit tests for sub-processing."""

import pathlib
import shutil
import tempfile
import unittest

from phamerate.mmseqs import *
from phamerate.multiprocess import CPUS

mmseqs_version = get_mmseqs_version()
if not mmseqs_version:
    print("MMseqs2 not found - skipping unittests in 'test_mmseqs'")
    exit()

file_dir = pathlib.Path(__file__).parent.parent.joinpath("files/mmseqs")

genes_fasta = file_dir.joinpath("genes.faa")
not_a_file = file_dir.joinpath("not_a_file.faa")


class TestMMseqs2(unittest.TestCase):
    def setUp(self):
        self.tempdir = pathlib.Path(tempfile.mkdtemp(prefix="test_mmseqs2__"))

    def test_createdb_valid_fasta(self):
        """Test that a FASTA file can be used to create an MMseqs2
        sequence database."""
        seq_db = self.tempdir.joinpath("sequenceDB")
        createdb(genes_fasta, seq_db)

        # Check for three of the sequenceDB files
        self.assertTrue(seq_db.is_file())
        self.assertTrue(seq_db.with_suffix(".dbtype").is_file())
        self.assertTrue(seq_db.with_suffix(".index").is_file())

    def test_createdb_invalid_fasta(self):
        """Test that a non-existent FASTA raises an MMseqs2Error."""
        seq_db = self.tempdir.joinpath("sequenceDB")
        with self.assertRaises(MMseqs2Error):
            createdb(not_a_file, seq_db)

        # Check for three of the sequenceDB files
        self.assertFalse(seq_db.is_file())
        self.assertFalse(seq_db.with_suffix(".dbtype").is_file())
        self.assertFalse(seq_db.with_suffix(".index").is_file())

    def test_cluster_good_params_1(self):
        """Test that clustering with valid parameters works.

        MMseqs2 expects identity and coverage values to be in [0, 1]
        so 90% identity and 90% coverage should both be auto-scaled
        to 0.9, and no error raised."""
        seq_db = self.tempdir.joinpath("sequenceDB")
        clu_db = self.tempdir.joinpath("clusterDB")
        createdb(genes_fasta, seq_db)
        cluster(seq_db, clu_db, self.tempdir, identity=90, coverage=90,
                e_value=0.001, sens=1, cpus=CPUS)

        # Check for two of the clusterDB files
        self.assertTrue(clu_db.with_suffix(".dbtype").is_file())
        self.assertTrue(clu_db.with_suffix(".index").is_file())

    def test_cluster_good_params_2(self):
        """Test that clustering with valid parameters works.

        MMseqs2 expects identity and coverage values to be in [0, 1]
        so 0.9 identity and 0.9 coverage should both be ok, and no
        error raised."""
        seq_db = self.tempdir.joinpath("sequenceDB")
        clu_db = self.tempdir.joinpath("clusterDB")
        createdb(genes_fasta, seq_db)
        cluster(seq_db, clu_db, self.tempdir, identity=0.9, coverage=0.9,
                e_value=0.001, sens=1, cpus=CPUS)

        # Check for two of the clusterDB files
        self.assertTrue(clu_db.with_suffix(".dbtype").is_file())
        self.assertTrue(clu_db.with_suffix(".index").is_file())

    def test_cluster_bad_params(self):
        """Test that clustering with invalid parameters raises an
        MMseqs2Error.

        MMseqs2 expects identity and coverage values to be in [0, 1];
        900% identity and 900% coverage are out of range, so will both
        be auto-scaled down to 9 - still out of bounds, so an error
        should be raised."""
        seq_db = self.tempdir.joinpath("sequenceDB")
        clu_db = self.tempdir.joinpath("clusterDB")
        createdb(genes_fasta, seq_db)
        with self.assertRaises(MMseqs2Error):
            cluster(seq_db, clu_db, self.tempdir, identity=900,
                    coverage=900, e_value=0.001, sens=1, cpus=CPUS)

        # Check for two of the clusterDB files
        self.assertFalse(clu_db.with_suffix(".dbtype").is_file())
        self.assertFalse(clu_db.with_suffix(".index").is_file())

    def test_createseqfiledb(self):
        """"""
        seq_db = self.tempdir.joinpath("sequenceDB")
        clu_db = self.tempdir.joinpath("clusterDB")
        sf_db = self.tempdir.joinpath("seqfileDB")

        createdb(genes_fasta, seq_db)
        cluster(seq_db, clu_db, self.tempdir, identity=90, coverage=90,
                e_value=0.001, sens=1, cpus=CPUS)
        createseqfiledb(seq_db, clu_db, sf_db)

        # Check for three of the seqfileDB files
        self.assertTrue(sf_db.is_file())
        self.assertTrue(sf_db.with_suffix(".dbtype").is_file())
        self.assertTrue(sf_db.with_suffix(".index").is_file())

    def test_result2flat(self):
        seq_db = self.tempdir.joinpath("sequenceDB")
        clu_db = self.tempdir.joinpath("clusterDB")
        sf_db = self.tempdir.joinpath("seqfileDB")
        seq_out = self.tempdir.joinpath("seqPhams.txt")

        createdb(genes_fasta, seq_db)
        cluster(seq_db, clu_db, self.tempdir, identity=90, coverage=90,
                e_value=0.001, sens=1, cpus=CPUS)
        createseqfiledb(seq_db, clu_db, sf_db)
        result2flat(seq_db, seq_db, sf_db, seq_out)

        self.assertTrue(seq_out.is_file())

    def tearDown(self):
        shutil.rmtree(self.tempdir)


if __name__ == '__main__':
    unittest.main()
