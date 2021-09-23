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
        self._translations = dict()

    @property
    def translations(self):
        return self._translations.keys()

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
        genomes = self._translations.get(translation, dict())
        if genome in genomes:
            genomes[genome].add(geneid)
        else:
            genomes[genome] = {geneid}
        self._translations[translation] = genomes

    def get_data(self, translation):
        data = list()
        for genome, geneids in self._translations[translation].items():
            for geneid in geneids:
                data.append((genome, geneid))
        return data

    def __iter__(self):
        for translation in self.translations:
            data = self.get_data(translation)
            yield translation, data

    def __len__(self):
        return sum([len(self.get_data(x)) for x in self.translations])


class Pham(SequenceDB):
    def __init__(self):
        super().__init__()

    def merge(self, other):
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        for translation, data in other:
            for genome, geneid in data:
                self.add_gene(geneid, genome, translation)

    def get_genomes(self):
        genomes = [x.keys() for x in self._translations.values()]
        return set.union(*genomes)

    def __contains__(self, item):
        if not isinstance(item, str):
            raise TypeError("Pham only contains data of type 'str'")

        return item in self._translations

    def __lt__(self, other):
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        return len(self) < len(other)

    def __str__(self):
        s = ""
        for translation, data in self:
            for genome, geneid in data:
                s += f">{geneid}\n{translation}\n"
        return s
