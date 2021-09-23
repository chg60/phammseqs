from phamerate.cmd import run_command


def clustalo_align(fasta, align):
    """
    Align a pham using Clustal Omega

    :param fasta: path to a FASTA file to align
    :type fasta: pathlib.Path
    :param align: path where the alignment should go
    :type align: pathlib.Path
    """
    command = f"clustalo -i {fasta} --infmt=fasta -o {align} " \
              f"--outfmt=fasta --output-order=tree-order --threads=1 " \
              f"--seqtype=protein --force"
    run_command(command)
