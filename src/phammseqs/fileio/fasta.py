"""Functions for reading and writing FASTA files."""


def read_fasta(filepath):
    """Read a FASTA file and return the headers and sequences.

    :param filepath: a FASTA file to parse
    :type filepath: pathlib.Path or str
    :return: sequences
    """
    headers, sequences = list(), list()

    cache = list()
    with open(filepath, "r") as fasta_reader:
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

    return zip(headers, sequences)


def write_fasta(sequences, filepath, wrap=80):
    """Write the sequences to the indicated file in FASTA format,
    wrapping long sequence lines at the indicated width.

    Use wrap=None or wrap=0 or wrap=False to turn off line-wrapping.

    :param sequences: the headers and translations to write_fasta
    :type sequences: list[tuple[str, str]]
    :param filepath: a FASTA file to write_fasta to
    :type filepath: pathlib.Path or str
    :param wrap: maximum number of characters per sequence line
    :type wrap: int
    :return: filepath
    """
    with open(filepath, "w") as fasta_writer:
        for header, sequence in sequences:
            fasta_writer.write(f">{header}\n")

            # If no wrapping, write_fasta the whole sequence
            if not wrap:
                fasta_writer.write(f"{sequence}\n")
            # Else, use string slicing to split the sequence to satisfy wrap
            else:
                for i in range(0, len(sequence), wrap):
                    fasta_writer.write(f"{sequence[i:i + wrap]}\n")

    return filepath
