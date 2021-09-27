"""Some classes that serve useful functions..."""


class SequenceDB:
    """
    Provide a convenient interface for wrangling translations and
    gene identifiers from different sources
    """
    def __init__(self):
        """
        Initialize an instance of SequenceDB.
        """
        self._translations_to_geneids = dict()
        self._geneids_to_translations = dict()
        self._genomes_to_geneids = dict()
        self._geneids_to_genomes = dict()

    @property
    def translations(self):
        """
        Provide read-only iterable access to the translations
        represented in this database.

        :return: translations
        """
        return self._translations_to_geneids.keys()

    @property
    def genomes(self):
        """
        Provide read-only iterable access to the genomes represented in
        this database.

        :return: genomes
        """
        return self._genomes_to_geneids.keys()

    @property
    def nr_iterator(self):
        """
        Alternative to __iter__() that allows iterating over non-redundant
        sequences rather than all sequences.

        :rtype: tuple[str, str]
        """
        for translation, geneids in self._translations_to_geneids.items():
            yield geneids[0], translation

    def add_gene(self, geneid, genome, translation):
        """
        Add a gene to the database.

        :param geneid: the gene's identifier
        :type geneid: str
        :param genome: the genome the gene comes from
        :type genome: str
        :param translation: the gene's protein sequence
        :type translation: str
        """
        if geneid in self._geneids_to_translations:
            raise ValueError(f"geneid {geneid} already in database")

        self._geneids_to_translations[geneid] = translation

        if translation in self._translations_to_geneids:
            self._translations_to_geneids[translation].append(geneid)
        else:
            self._translations_to_geneids[translation] = [geneid]

        self._geneids_to_genomes[geneid] = genome

        if genome in self._genomes_to_geneids:
            self._genomes_to_geneids[genome].append(geneid)
        else:
            self._genomes_to_geneids[genome] = [geneid]

    def get_translation_geneids(self, translation):
        """
        Lookup the indicated translation and return its associated
        geneids.

        :param translation: the translation to lookup
        :type translation: str
        :return: geneids
        :rtype: list[str]
        :raise: KeyError
        """
        return self._translations_to_geneids[translation]

    def get_geneid_translation(self, geneid):
        """
        Lookup the indicated geneid and return its translation.

        :param geneid: the geneid to lookup
        :type geneid: str
        :return: translation
        :rtype: str
        :raise: KeyError
        """
        return self._geneids_to_translations[geneid]

    def get_genome_geneids(self, genome):
        """
        Lookup the indicated genome and return its associated geneids.

        :param genome: the genome to lookup
        :type genome: str
        :return: geneids
        :rtype: list[str]
        :raise: KeyError
        """
        return self._genomes_to_geneids[genome]

    def get_geneid_genome(self, geneid):
        """
        Lookup the indicated geneid and return its genome.

        :param geneid: the geneid to lookup
        :type geneid: str
        :return: genome
        :rtype: str
        :raise: KeyError
        """
        return self._geneids_to_genomes[geneid]

    def get_genome_translations(self, genome):
        """
        Lookup the indicated genome and return its associated
        translations.

        :param genome: the genome to lookup
        :type genome: str
        :return: translations
        :rtype: list[str]
        :raise: KeyError
        """
        return [self.get_geneid_translation(x) for x in
                self.get_genome_geneids(genome)]

    def __iter__(self):
        for geneid, translation in self._geneids_to_translations.items():
            genome = self.get_geneid_genome(geneid)
            yield genome, geneid, translation

    def __len__(self):
        return len(self._geneids_to_translations)


class Pham(SequenceDB):
    """
    A type of SequenceDB, with additional convenience methods for
    performing actions that are useful when considering relatively
    small numbers of related sequences, and not so much for whole
    meta- or pan-genomes.
    """
    def __init__(self):
        super().__init__()
        self._geneids_to_aligned_translations = dict()

    @property
    def alignment_iterator(self):
        """
        Alternative to __iter__(), which allows iteration over aligned
        sequences rather than the unaligned sequences.
        """
        alignments = self._geneids_to_aligned_translations.items()
        for geneid, translation in alignments:
            genome = self.get_geneid_genome(geneid)
            yield genome, geneid, translation

    def merge(self, other):
        """
        Merge phams from multiple rounds of clustering.

        :param other: another pham to be merged into this one
        :type other: Pham
        """
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        for genome, geneid, translation in other:
            try:
                self.add_gene(geneid, genome, translation)
            except ValueError:
                pass    # Tried adding a gene that was already in the pham

    def add_aligned_translation(self, aligned_translation):
        """
        Add an aligned translation for this database to keep track of.

        :param aligned_translation: an aligned translation
        :type aligned_translation: str
        :raise: KeyError
        """
        translation = aligned_translation.replace("-", "")
        geneids = self.get_translation_geneids(translation)
        for geneid in geneids:
            self._geneids_to_aligned_translations[geneid] = aligned_translation

    def in_genome(self, genome):
        """
        Check whether this phamily is present in the indicated genome.

        :param genome: the genome to be checked against
        :type genome: str
        return: in_genome
        """
        return genome in self.genomes

    def minimum_length(self, kind="nucleotide"):
        """
        Return the length of the shortest sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: min_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        min_length = min(lengths)

        if kind == "nucleotide":
            return min_length * 3 + 3

        return min_length

    def average_length(self, kind="nucleotide"):
        """
        Return the length of the average sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: avg_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        avg_length = float(sum(lengths))/len(lengths)

        if kind == "nucleotide":
            return avg_length * 3 + 3

        return avg_length

    def maximum_length(self, kind="nucleotide"):
        """
        Return the length of the longest sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: max_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        max_length = max(lengths)

        if kind == "nucleotide":
            return max_length * 3 + 3

        return max_length

    def __lt__(self, other):
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        return len(self) < len(other)
