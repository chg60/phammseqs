"""A class for efficient storage and utilization of large numbers of
sequences."""


class SequenceDB:
    """A convenient interface for wrangling translations and
    gene identifiers from different sources.
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
        """Provide read-only iterable access to the translations
        represented in this database.

        :return: translations
        """
        return self._translations_to_geneids.keys()

    @property
    def genomes(self):
        """Provide read-only iterable access to the genomes represented
        in this database.

        :return: genomes
        """
        return self._genomes_to_geneids.keys()

    @property
    def nr_iterator(self):
        """Alternative to __iter__() that allows iterating over
        non-redundant sequences rather than all sequences.

        :rtype: tuple[str, str]
        """
        for translation, geneids in self._translations_to_geneids.items():
            yield geneids[0], translation

    def add_gene(self, geneid, genome, translation):
        """Add a gene to the database.

        :param geneid: the gene's identifier
        :type geneid: str
        :param genome: the genome the gene comes from
        :type genome: str
        :param translation: the gene's protein sequence
        :type translation: str
        """
        if geneid in self._geneids_to_translations:
            raise ValueError(f"geneid '{geneid}' already in database")

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
        """Lookup the indicated translation and return its associated
        geneids.

        :param translation: the translation to lookup
        :type translation: str
        :return: geneids
        :rtype: list[str]
        :raise: KeyError
        """
        return self._translations_to_geneids[translation]

    def get_geneid_translation(self, geneid):
        """Lookup the indicated geneid and return its translation.

        :param geneid: the geneid to lookup
        :type geneid: str
        :return: translation
        :rtype: str
        :raise: KeyError
        """
        return self._geneids_to_translations[geneid]

    def get_genome_geneids(self, genome):
        """Lookup the indicated genome and return its associated geneids.

        :param genome: the genome to lookup
        :type genome: str
        :return: geneids
        :rtype: list[str]
        :raise: KeyError
        """
        return self._genomes_to_geneids[genome]

    def get_geneid_genome(self, geneid):
        """Lookup the indicated geneid and return its genome.

        :param geneid: the geneid to lookup
        :type geneid: str
        :return: genome
        :rtype: str
        :raise: KeyError
        """
        return self._geneids_to_genomes[geneid]

    def get_genome_translations(self, genome):
        """Lookup the indicated genome and return its associated
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
