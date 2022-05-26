from phamerate.fileio.fasta import read_fasta, write_fasta
from phamerate.fileio.genbank import read_genbank


def sniff_format(filepath):
    """Sniff the first line of the file to determine the file format.

    Return options are: "fasta", "genbank", or None

    :param filepath: the path to the file we need to sniff
    :type filepath: pathlib.Path
    :return: fmt
    """
    with open(filepath, "r") as file_sniffer:
        line = file_sniffer.readline()

    if line.startswith(">"):
        return "fasta"
    elif line.startswith("LOCUS"):
        return "genbank"
    else:
        return None


__all__ = ["read_fasta", "read_genbank", "sniff_format", "write_fasta"]
