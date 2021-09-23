"""Functions for writing and parsing FASTA files."""


def write_fasta(headers, sequences, filepath, width=80):
    """
    Write the given headers and sequences to filepath, wrapping long
    sequence lines at the indicated width.

    :param headers: sequence labels
    :type headers: str or list of str
    :param sequences: nucleotide or amino acid sequences
    :type sequences: str or list of str
    :param filepath: a FASTA file to write to
    :type filepath: pathlib.Path or str
    :param width: maximum number of characters per sequence line
    :type width: int
    :return: filepath
    """
    if not isinstance(headers, list) or not isinstance(sequences, list):
        raise TypeError(f"headers and sequences should be lists of strings")

    # Open the file for writing
    fasta_writer = open(filepath, "w")

    for header, sequence in zip(headers, sequences):
        fasta_writer.write(f">{header}\n")
        # Use string slicing to split the sequence to satisfy the width param
        for i in range(0, len(sequence), width):
            fasta_writer.write(f"{sequence[i:i + width]}\n")

    # Close the file
    fasta_writer.close()

    return filepath


def parse_fasta(filepath):
    """
    Parse a FASTA file and return the headers and sequences

    :param filepath: a FASTA file to parse
    :type filepath: pathlib.Path or str
    :return: headers, sequences
    """
    headers, sequences = list(), list()

    # Open the file for reading
    fasta_reader = open(filepath, "r")

    cache = list()
    for line in fasta_reader:
        # If the line is a header line, flush the cache and store new header
        if line.startswith(">"):
            sequences.append("".join(cache))
            cache = list()
            headers.append(line.lstrip(">").rstrip())
        # Otherwise, append to the cache
        else:
            cache.append(line.rstrip())
    # Flush the last sequence out of the cache, and pop the empty sequence
    sequences.append("".join(cache))
    sequences = sequences[1:]

    # Close the file
    fasta_reader.close()

    return headers, sequences
