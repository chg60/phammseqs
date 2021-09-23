from phamerate.cmd import run_command


def clustalo_align(fasta, out_dir):
    """
    Align a pham using Clustal Omega

    :param fasta:
    :type fasta: pathlib.Path
    :param out_dir: path to a FASTA file
    :type out_dir: pathlib.Path
    """
    align = out_dir.joinpath(f"{fasta.stem}.aln")
    command = f"clustalo -i {fasta} --infmt=fasta -o {align} " \
              f"--outfmt=fasta --output-order=tree-order --threads=1 " \
              f"--seqtype=protein --force"
    run_command(command)
