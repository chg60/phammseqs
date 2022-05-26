"""Function for reading Genbank flatfiles."""

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def read_genbank(filepath):
    """Read a Genbank flatfile and return headers and sequences
    for all CDS features.

    :param filepath: the path to a Genbank flatfile
    :type filepath: pathlib.Path or str
    :return: headers, sequences
    """
    headers, sequences = list(), list()

    for record in SeqIO.parse(filepath, "genbank"):
        features = [x for x in record.features if x.type == "CDS"]

        for i, cds in enumerate(features):
            header = f"lcl|{record.id}_"

            # try to get the "protein_id" qualifier, otherwise skip
            protein_id = cds.qualifiers.get("protein_id", [None])[0]
            if protein_id:
                header += f"prot_{protein_id}_"

            # append CDS index + 1
            header += f"{i + 1}"

            # try to get the "gene" qualifier, otherwise skip
            gene = cds.qualifiers.get("gene", [None])[0]
            if gene:
                header += f" [gene={gene}]"

            # try to get the "locus_tag" qualifier, otherwise skip
            locus_tag = cds.qualifiers.get("locus_tag", [None])[0]
            if locus_tag:
                header += f" [locus_tag={locus_tag}]"

            # try to get the "product" qualifier, otherwise skip
            product = cds.qualifiers.get("product", [None])[0]
            if product:
                header += f" [protein={product}]"
            else:
                header += " [protein=hypothetical protein]"

            # append the protein_id
            if protein_id:
                header += f" [protein_id={protein_id}]"

            # append the location
            start, end = cds.location.start, cds.location.end
            if cds.location.strand == 1:
                header += f" [location={start}..{end}]"
            else:
                header += f" [location=complement({start}..{end})]"

            # append the genbank key
            header += f" [gbkey=CDS]"

            # try to get the translation or generate it - skip if cannot
            try:
                translation = str(cds.qualifiers.get("translation", [None])[0])
                if not translation:
                    translation = str(cds.translate(record.seq))
            except TranslationError:
                # Translations encoded by frameshifts or those that use
                # suppressor tRNAs may fail to translate - skip these
                print(f"Unable to translate {header} - skipping it")
                continue

            headers.append(header)
            sequences.append(translation)

    return zip(headers, sequences)
