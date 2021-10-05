from phamerate.cmd import run_command
from phamerate.fileio import parse_fasta, write_fasta


def clustalo_align(fasta, align):
    """
    Use clustalo to align the sequences in `fasta` and write the
    alignment to `align`.

    :param fasta: path to a FASTA multiple sequence file
    :type fasta: pathlib.Path
    :param align: path to a FASTA multiple sequence alignment file
    :type align: pathlib.Path
    """
    command = f"clustalo -i {fasta} --infmt=fasta -o {align} " \
              f"--outfmt=fasta --output-order=tree-order --threads=1 " \
              f"--seqtype=protein --force"
    run_command(command)


def align_phams(phams, fasta_paths, align_paths):
    """
    Align a pham using Clustal Omega

    :param phams: the phams to align
    :type phams: list of classes.Pham
    :param fasta_paths: paths to FASTA files to align
    :type fasta_paths: list of pathlib.Path
    :param align_paths: paths where the alignments should go
    :type align_paths: list of pathlib.Path
    """
    for pham, fasta, align in zip(phams, fasta_paths, align_paths):
        # Only align phams with > 1 non-redundant sequence
        if len([x for x in pham.nr_iterator]) != 1:
            clustalo_align(fasta, align)

            # Make sure all this pham's headers are in the alignment file
            _, aligned_translations = parse_fasta(align)
            for aligned_translation in aligned_translations:
                pham.add_aligned_translation(aligned_translation)

            headers, sequences = list(), list()
            for _, geneid, translation in pham.alignment_iterator:
                headers.append(geneid), sequences.append(translation)
            write_fasta(headers, sequences, align)

        # Make sure all this pham's headers are in the FASTA file
        headers, sequences = list(), list()
        for _, geneid, translation in pham:
            headers.append(geneid), sequences.append(translation)
        write_fasta(headers, sequences, fasta)
