from Bio import SeqIO


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
                headers.append(feature.qualifiers["locus_tag"][0])
                sequences.append(feature.qualifiers["translation"][0])

    return headers, sequences
