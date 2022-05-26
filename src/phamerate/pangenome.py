"""Classes and functions to perform basic pangenome analysis."""

from phamerate.SequenceDB import SequenceDB, Pham
from phamerate import fileio

summary_template = """
Core genes (99% <= strains <= 100%):    {}
Soft core genes (95% <= strains < 99%): {}
Shell genes (15% <= strains < 95%):     {}
Cloud genes (0% <= strains < 15%):      {}
Total genes (0% <= strains <= 100%):    {}
"""

ROARY_HEADER = ["Gene", "Non unique Gene name", "No. isolates", "No. sequences",
                "Avg sequences per isolate", "Genome Fragment",
                "Order within Fragment", "Accessory Fragment",
                "Accessory Order with Fragment", "QC", "Min group size nuc",
                "Max group size nuc", "Avg group size nuc"]


class Genome(SequenceDB):
    """A class derived from SequenceDB that is better suited for the
    kinds of operations one might wish to do with smaller groups of
    sequences found in the same genome.
    """
    def __init__(self, sequences=None):
        super().__init__(sequences)
        self._pham_geneids = dict()

    def add_phams(self, phams):
        geneids = set(self.geneids)
        for i, pham in enumerate(phams):
            phamid = f"pham_{i+1}"
            shared = geneids.intersection(pham.geneids)
            if len(shared) > 0:
                self._pham_geneids[phamid] = list(shared)

    def __contains__(self, item):
        """Query a genome for whether it contains a pham."""
        return item in self._pham_geneids.keys()

    @property
    def pham_genes(self):
        for phamid, geneids in self._pham_geneids.items():
            yield phamid, geneids


def build_genomes(phams, infiles):
    """Parse input files into genomes, and match them to the phams
    their genes belong to.

    :param phams:
    :type phams: list[Pham]
    :param infiles:
    :type infiles: list[pathlib.Path]
    :return: temp_genomes
    """
    temp_genomes = dict()
    for infile in infiles:
        name = infile.stem

        fmt = fileio.sniff_format(infile)
        if fmt == "fasta":
            sequences = fileio.read_fasta(infile)
        elif fmt == "genbank":
            sequences = fileio.read_genbank(infile)
        else:
            continue

        genome = Genome(sequences)
        genome.add_phams(phams)

        temp_genomes[name] = genome

    return temp_genomes


def analyze_pangenome(phams, infiles, outdir):
    """

    :param phams:
    :type phams: list[Pham]
    :param infiles:
    :type infiles: list[pathlib.Path]
    :param outdir:
    :type outdir: pathlib.Path
    :return:
    """
    # Build genomes to analyze
    genomes = build_genomes(phams, infiles)

    # Map phams to genomes to geneids
    pham_data = dict()
    for name, genome in genomes.items():
        for phamid, geneids in genome.pham_genes:
            if phamid in pham_data:
                pham_data[phamid][name] = geneids
            else:
                pham_data[phamid] = {name: geneids}

    # Analyze pangenome
    core, soft, shell, cloud = list(), list(), list(), list()
    strain_rows, roary_rows = list(), list()

    for i, pham in enumerate(phams):
        phamid = f"pham_{i+1}"
        pham_genomes = pham_data[phamid]

        # Begin building this pham's Roary row
        seq_per_iso = float(len(pham)) / len(pham_genomes)
        roary_row = [phamid, "", str(len(pham_genomes)), str(len(pham)),
                     str(seq_per_iso), "", "", "", "", "",
                     str(3 * pham.minimum_length + 3),
                     str(3 * pham.maximum_length + 3),
                     str(3 * pham.average_length + 3)]

        # Deal with Roary rows and strain genes
        for name in genomes.keys():
            if name in pham_genomes:
                geneids = pham_genomes[name]
                roary_row.append(";".join([x.split()[0] for x in geneids]))
                for _ in geneids:
                    strain_rows.append((name, phamid))
            else:
                roary_row.append("")

        roary_rows.append(roary_row)

        fraction = float(len(pham_genomes)) / len(genomes)
        if fraction >= 0.99:
            core.append(phamid)
        elif fraction >= 0.95:
            soft.append(phamid)
        elif fraction >= 0.15:
            shell.append(phamid)
        else:
            cloud.append(phamid)

    summary_file = outdir.joinpath("summary_statistics.txt")
    with open(summary_file, "w") as summary_writer:
        summary = summary_template.format(len(core), len(soft), len(shell),
                                          len(cloud), len(phams))
        summary_writer.write(summary)

    roary_file = outdir.joinpath("gene_presence_absence.csv")
    with open(roary_file, "w") as roary_writer:
        header = ROARY_HEADER
        header.extend(genomes.keys())
        roary_writer.write(",".join(header) + "\n")
        for roary_row in roary_rows:
            roary_writer.write(",".join(roary_row) + "\n")

    strain_genes_file = outdir.joinpath("strain_genes.tsv")
    with open(strain_genes_file, "w") as strain_genes_writer:
        for name, phamid in strain_rows:
            strain_genes_writer.write(f"{name}\t{phamid}\n")
