"""Functions for running Clustal Omega"""

from phamerate import subprocess

VERSION = "1.2.4"


class ClustalOmegaError(Exception):
    """Error caused by running Clustal Omega."""
    pass


def get_clustalo_version():
    """Return the Clustal Omega version string.

    :return: version
    """
    command = "clustalo --version"
    try:
        stdout, stderr = subprocess.run(command)
        if stderr:
            print(stderr)
            raise ClustalOmegaError(f"command failed: {command}")
        return stdout.rstrip()
    except FileNotFoundError:
        return None


def run_clustalo(infile, outfile, threads=1, verbose=False):
    """Use clustalo to generate an MSA from the sequences in `infile`
    and store the results in `outfile`.

    :param infile: path to a FASTA multiple sequence file
    :type infile: pathlib.Path
    :param outfile: path to a FASTA multiple sequence alignment file
    :type outfile: pathlib.Path
    :param threads: number of threads to use
    :type threads: int
    :param verbose: return stdout
    :type verbose: bool
    :return: infile, outfile
    """
    command = f"clustalo -i {infile} --infmt=fasta -o {outfile} " \
              f"--outfmt=fasta --output-order=tree-order " \
              f"--threads={threads} --seqtype=protein --force"

    stdout, stderr = subprocess.run(command)
    if verbose and stdout:
        print(stdout)
    if stderr and not stderr.startswith("WARNING: "):
        print(stderr)
        raise ClustalOmegaError(f"command failed: {command}")

    return infile, outfile
