"""Assort phage protein sequences into phamilies of homologs using MMseqs2."""

import argparse
import datetime
import pathlib

from phamerate.multiprocess import CPUS

DATE = datetime.datetime.now().strftime("%d_%b_%Y")
OUT_DIR = pathlib.Path().cwd().joinpath(f"phamerate__{DATE}")

# Sequence-sequence
SEQ_MIN_SEQ_ID = 35.0
SEQ_COVERAGE = 80.0
SEQ_EVALUE = 0.001
SEQ_SENSITIVITY = 7
SEQ_CLUSTER_MODE = 0
SEQ_CLUSTER_STEPS = 1

HMM_MIN_SEQ_ID = 15.0
HMM_COVERAGE = 70.0
HMM_EVALUE = 0.001
HMM_SENSITIVITY = 7
HMM_CLUSTER_MODE = 0
HMM_CLUSTER_STEPS = 3

EPILOG = """
Steinegger M. and Söding J. MMseqs2 enables sensitive protein
sequence searching for the analysis of massive data sets. Nature
Biotechnology, 2017. doi: 10.1038/nbt.3988"""


def parse_args():
    """Parse command line arguments."""
    p = argparse.ArgumentParser(description=__doc__, prog="phamerate",
                                epilog=EPILOG)

    p.add_argument("infile", nargs="+", type=pathlib.Path,
                   help="path to input file(s) in FASTA or Genbank flatfile "
                        "format")

    s = p.add_argument_group("MMseqs2 sequence-sequence clustering arguments")
    s.add_argument("--identity",
                   type=float, default=SEQ_MIN_SEQ_ID, metavar='',
                   help=f"percent identity for sequence-sequence clustering "
                        f"[default: {SEQ_MIN_SEQ_ID}%%]")
    s.add_argument("--coverage",
                   type=float, default=SEQ_COVERAGE, metavar='',
                   help=f"percent coverage for sequence-sequence clustering "
                        f"[default: {SEQ_COVERAGE}%%]")
    s.add_argument("--evalue",
                   type=float, default=SEQ_EVALUE, metavar='',
                   help=f"E-value threshold for sequence-sequence clustering "
                        f"[default: {SEQ_EVALUE}]")
    s.add_argument("--sensitivity",
                   type=float, default=SEQ_SENSITIVITY, metavar='',
                   help=f"sensitivity: 1 favors speed, 7 favors "
                        f"sensitivity [default: {SEQ_SENSITIVITY}]")
    s.add_argument("--cluster-mode",
                   type=int, default=SEQ_CLUSTER_MODE, metavar='',
                   help=f"clustering algorithm [default: {SEQ_CLUSTER_MODE}]")
    s.add_argument("--cluster-steps",
                   type=int, default=SEQ_CLUSTER_STEPS, metavar='',
                   help=f"number of steps for sequence-sequence clustering "
                        f"to proceed in [default: {SEQ_CLUSTER_STEPS}]")

    h = p.add_argument_group("MMseqs2 profile-sequence clustering arguments")
    h.add_argument("--hmm-identity",
                   type=float, default=HMM_MIN_SEQ_ID, metavar='',
                   help=f"percent identity for profile-consensus clustering "
                        f"[default: {HMM_MIN_SEQ_ID}%%]")
    h.add_argument("--hmm-coverage",
                   type=float, default=HMM_COVERAGE, metavar='',
                   help=f"percent coverage for profile-consensus clustering "
                        f"[default: {HMM_COVERAGE}%%]")
    h.add_argument("--hmm-evalue",
                   type=float, default=HMM_EVALUE, metavar='',
                   help=f"E-value threshold for profile-consensus clustering "
                        f"[default: {HMM_EVALUE}]")
    h.add_argument("--hmm-sensitivity",
                   type=float, default=HMM_SENSITIVITY, metavar='',
                   help=f"sensitivity: 1 favors speed, 7 favors "
                        f"sensitivity [default: {HMM_SENSITIVITY}]")
    h.add_argument("--hmm-cluster-mode",
                   type=int, default=HMM_CLUSTER_MODE, metavar='',
                   help=f"clustering algorithm [default: {HMM_CLUSTER_MODE}]")
    h.add_argument("--hmm-cluster-steps",
                   type=int, default=HMM_CLUSTER_STEPS, metavar='',
                   help=f"number of steps for profile-consensus clustering "
                        f"to proceed in [default: {HMM_CLUSTER_STEPS}]")
    h.add_argument("--skip-hmm", action="store_true",
                   help="do not perform profile-consensus clustering")

    p.add_argument("-v", "--verbose", action="store_true",
                   help="print progress messages to the console")
    p.add_argument("-d", "--debug", action="store_true",
                   help="run in debug mode")
    p.add_argument("-a", "--align-phams", action="store_true",
                   help="use Clustal Omega to align phams (could take "
                        "awhile...)")
    p.add_argument("-p", "--pangenome", action="store_true",
                   help="pangenome analysis à la Roary (only meaningful if "
                        "given one input file per genome)")
    p.add_argument("-c", "--cpus", type=int, default=CPUS,
                   help=f"number of threads to use [default: {CPUS}]")
    p.add_argument("-o", "--outdir", default=OUT_DIR, type=pathlib.Path,
                   help=f"path to directory where output files should go "
                        f"[default: {OUT_DIR}]")

    return p.parse_args()
