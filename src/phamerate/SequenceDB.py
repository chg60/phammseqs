"""Classes for efficient storage and utilization of large numbers of
sequences."""

from phamerate.fileio import fasta


class SequenceDB:
    """A convenient interface for wrangling sequences and labels from
    different sources."""
    def __init__(self, sequences=None):
        """Initialize an instance of SequenceDB.

        :param sequences: geneids to add to the database
        :type sequences: list[tuple[str, str]]
        """
        self._translation_geneids: dict[str, list[str]] = dict()
        self._geneid_translations: dict[str, str] = dict()

        if sequences:
            for geneid, translation in sequences:
                self.add_gene(geneid, translation)

    def add_gene(self, geneid, translation):
        """Add a gene to the database.

        :param geneid: the gene's identifier
        :type geneid: str
        :param translation: the gene's protein sequence
        :type translation: str
        """
        if geneid in self._geneid_translations:
            raise ValueError(f"geneid '{geneid}' already in database")

        self._geneid_translations[geneid] = translation

        if translation in self._translation_geneids:
            self._translation_geneids[translation].append(geneid)
        else:
            self._translation_geneids[translation] = [geneid]

    @property
    def geneids(self):
        """Provide access to the geneids in this database.

        :return: geneids
        """
        return self._geneid_translations.keys()

    @property
    def translations(self):
        """Provide access to the translations in this database.

        :return: translations
        """
        return self._translation_geneids.keys()

    @property
    def nr_genes(self):
        """Alternative to __iter__() that allows iterating over
        non-redundant sequences rather than all sequences."""
        for translation, geneids in self._translation_geneids.items():
            yield geneids[0], translation

    def get_translation_geneids(self, translation):
        """Lookup the indicated translation and return its associated
        geneids.

        :param translation: the translation to lookup
        :type translation: str
        :return: geneids
        """
        return self._translation_geneids[translation]

    def get_geneid_translation(self, geneid):
        """Lookup the indicated geneid and return its translation.

        :param geneid: the geneid to lookup
        :type geneid: str
        :return: translation
        :rtype: str
        :raise: KeyError
        """
        return self._geneid_translations[geneid]

    def load(self, filepath, reset=False):
        """Load sequences from a FASTA file.

        :param filepath: the file to load sequences from
        :type filepath: pathlib.Path
        :param reset: clear existing sequences first
        :type reset: bool
        """
        if reset:
            self.__init__()
        for geneid, translation in fasta.read_fasta(filepath):
            self.add_gene(geneid, translation)

    def write(self, filepath, nr_only=False):
        """Dump sequences to a FASTA file.

        :param filepath: path to the file to write_fasta sequences to
        :type filepath: pathlib.Path
        :param nr_only: only write_fasta non-redundant sequences to the file
        :type nr_only: bool
        :return: filepath
        """
        sequences = list()
        if nr_only:
            for geneid, translation in self.nr_genes:
                sequences.append((geneid, translation))
        else:
            for geneid, translation in self:
                sequences.append((geneid, translation))

        return fasta.write_fasta(sequences, filepath)

    def __iter__(self):
        for geneid, translation in self._geneid_translations.items():
            yield geneid, translation

    def __len__(self):
        return len(self._geneid_translations)


class Pham(SequenceDB):
    """A class derived from SequenceDB that is better suited for the
    kinds of operations one might wish to do with smaller groups of
    related sequences.
    """
    def __init__(self, sequences=None):
        super().__init__(sequences)

    def merge(self, other):
        """Merge phams from multiple rounds of clustering.

        :param other: another pham to be merged into this one
        :type other: Pham
        """
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        for geneid, translation in other:
            try:
                self.add_gene(geneid, translation)
            except ValueError:
                pass    # Tried adding a gene that was already in the pham

    @property
    def minimum_length(self):
        """Return the length of the shortest sequence in the phamily."""
        lengths = [len(x) for x in self.translations]
        return min(lengths)

    @property
    def average_length(self):
        """Return the length of the average sequence in the phamily."""
        lengths = [len(x) for x in self.translations]
        return float(sum(lengths))/len(lengths)

    @property
    def maximum_length(self):
        """Return the length of the longest sequence in the phamily."""
        lengths = [len(x) for x in self.translations]
        return max(lengths)

    @property
    def n_eff(self):
        """Return the effective size of the pham (number of nr genes)."""
        return len(self.translations)

    @property
    def is_orpham(self):
        """Check whether this pham is an orpham (only 1 gene)."""
        return len(self) == 1

    def update_msa(self, filepath):
        """Introduce duplicate sequences to the MSA in FASTA format.

        :param filepath: path to the MSA file to update
        :type filepath: pathlib.Path
        """
        nr_sequences = fasta.read_fasta(filepath)

        sequences = list()
        for nr_geneid, align in nr_sequences:
            geneids = self.get_translation_geneids(align.replace("-", ""))
            for geneid in geneids:
                sequences.append((geneid, align))

        return fasta.write_fasta(sequences, filepath)

    def __lt__(self, other):
        if not isinstance(other, Pham):
            raise TypeError(f"cannot compare Pham with {type(other)}")

        return len(self) < len(other)
