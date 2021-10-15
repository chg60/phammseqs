from Bio import SeqIO


def sniff(filepath):
    """
    Sniffs the first line of the file to determine the file format.

    Options are: "fasta", "genbank", or None

    :param filepath: the path to the file we need to sniff
    :type filepath: pathlib.Path
    :return: fmt
    """
    fmt = None
    with open(filepath, "r") as file_sniffer:
        line = file_sniffer.readline()

    if line.startswith(">"):
        fmt = "fasta"
    elif line.startswith("LOCUS"):
        fmt = "genbank"

    return fmt


def parse_genbank(filepath):
    """
    Parse CDS features from records in a Genbank flatfile.

    :param filepath:
    :return: headers, sequences
    """
    headers, sequences = list(), list()
    for record in SeqIO.parse(filepath, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                except KeyError:
                    locus_tag = ""
                try:
                    protein_id = feature.qualifiers["protein_id"][0]
                except KeyError:
                    protein_id = ""
                try:
                    product = feature.qualifiers["product"][0]
                except KeyError:
                    product = "hypothetical protein"
                try:
                    translation = feature.qualifiers["translation"][0]
                except KeyError:
                    translation = feature.translate(record.seq)

                if not locus_tag and not protein_id:
                    print("Please submit a GitHub issue about this gene:")
                    print(feature)
                    print("No 'protein_id' or 'locus_tag' qualifiers...")
                    continue

                if not locus_tag:
                    header = f"{protein_id} {product}"
                else:
                    header = f"{locus_tag} {product}"

                headers.append(header)
                sequences.append(translation)

    return headers, sequences


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
